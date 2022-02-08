#!/usr/bin/env python
# coding: utf-8

# # Constellation Design
# Open an STK scenario with the moon as the central body and then run this code. It will find the open instance of STK and then create our desired constellation.

# STK library imports
from agi.stk12.stkdesktop import STKDesktop
from agi.stk12.stkobjects import *
from agi.stk12.stkutil import *
from agi.stk12.vgt import *
# if using astrogator uncomment the below
# from agi.stk12.stkobjects.astrogator
# if using aviator uncomment the below
# from agi.stk12.stkobjects.aviator

# Python helper library imports
import os
import time

# Get reference to the current instance of STK
STK_PID = os.getenv('STK_PID')
stk = STKDesktop.AttachToApplication(pid=STK_PID)

# Grab a handle on the STK application root.
root = stk.Root 
type(root)
# will print the scene root


# 1. Define a scenario object
scenario = root.CurrentScenario
# 2. Set the analytical time period.
scenario.SetTimePeriod('Today','+24hr')
# 3. Reset the animation time to the newly established start time.
root.Rewind()


# Create constellation object
constellation = scenario.Children.New(AgESTKObjectType.eConstellation, "SatConstellation")


# Insert the constellation of Satellites
numOrbitPlanes = 4
numSatsPerPlane = 4


# Create satellite objects in the constellation itself
root.BeginUpdate()
for orbitPlaneNum, RAAN in enumerate(range(0,180,180//numOrbitPlanes),1): #RAAN in degrees

    for satNum, trueAnomaly in enumerate(range(0,360,360//numSatsPerPlane), 1): #trueAnomaly in degrees
        
        # Insert satellite
        satellite = scenario.Children.New(AgESTKObjectType.eSatellite, f"Sat{orbitPlaneNum}{satNum}")
                
        # Select Propagator
        satellite.SetPropagatorType(AgEVePropagatorType.ePropagatorTwoBody)
        
        # Set initial state
        twoBodyPropagator = satellite.Propagator
        keplarian = twoBodyPropagator.InitialState.Representation.ConvertTo(AgEOrbitStateType.eOrbitStateClassical.eOrbitStateClassical)
        
        # Set semi major axis and eccentricity
        keplarian.SizeShapeType = AgEClassicalSizeShape.eSizeShapeSemimajorAxis
        keplarian.SizeShape.SemiMajorAxis = 10000 #km
        keplarian.SizeShape.Eccentricity = 0
        
        # set inclination, arg of perigee, and RAAN
        keplarian.Orientation.Inclination = 55 #degrees
        keplarian.Orientation.ArgOfPerigee = 0 #degrees
        keplarian.Orientation.AscNodeType = AgEOrientationAscNode.eAscNodeRAAN
        keplarian.Orientation.AscNode.Value = RAAN  #degrees
        
        # set true Anomaly (with phasing)
        keplarian.LocationType = AgEClassicalLocation.eLocationTrueAnomaly
        keplarian.Location.Value = trueAnomaly + (360//numSatsPerPlane/2)*(orbitPlaneNum%2)  #Stagger true anomalies (degrees) for every other orbital plane       
        
        # Propagate
        satellite.Propagator.InitialState.Representation.Assign(keplarian)
        satellite.Propagator.Propagate()
        
        # Add to constellation object
        constellation.Objects.AddObject(satellite)
root.EndUpdate()


# Create coverage definition
coverageDefinition = scenario.Children.New(AgESTKObjectType.eCoverageDefinition, "CoverageDefinition")

# Set grid bounds type
grid = coverageDefinition.Grid
grid.BoundsType = AgECvBounds.eBoundsGlobal

# Set resolution
grid.ResolutionType = AgECvResolution.eResolutionDistance
resolution = grid.Resolution
resolution.Distance = 100

# Add constellation as asset
coverageDefinition.AssetList.Add("Constellation/SatConstellation")
coverageDefinition.ComputeAccesses()


# Create figure of merit
figureOfMerit = coverageDefinition.Children.New(AgESTKObjectType.eFigureOfMerit, "SimpleCoverage")

# Set the definition and compute type
# this is what we may need to look more into when computing coverage, number of sats in view of the surface, etc.
figureOfMerit.SetDefinitionType(AgEFmDefinitionType.eFmSimpleCoverage)
definition = figureOfMerit.Definition
definition.Satisfaction

fomDataProvider = figureOfMerit.DataProviders.GetDataPrvFixedFromPath("Overall Value")
fomResults = fomDataProvider.Exec()

percentCoverage = fomResults.DataSets.GetDataSetByName("Minimum").GetValues()[0]

# Print data to console
print('The percent of the surface that is covered is: ', percentCoverage * 100, '%')

