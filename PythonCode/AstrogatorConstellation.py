import os
import time
import random as r

from agi.stk12.stkdesktop import STKDesktop
from agi.stk12.stkobjects import *
from agi.stk12.stkobjects.astrogator import *
from agi.stk12.utilities.colors import *
from math import *
import numpy as np
from pprint import pprint

#MOON VERSION
class orbitalprop:
    def __init__(self,stknew):
        """
        Parameters
        ----------
        stknew : Boolean, False if STK is running.
    
        Returns
        -------
        None.
        """
        if stknew:
            #Launches a new instance of STK.
            self.stk = STKDesktop.StartApplication(visible=True, userControl=True)
        else:
            #Finds STK and connects to it if it is already running.
            STK_PID = os.getenv('STK_PID')
            self.stk = STKDesktop.AttachToApplication(pid=STK_PID)
        
    def CreateScenario(self,scenarionew,timestart,duration):
        """
        Definition used to generate the scenario. Also can be used to generate new scenarios.
        Parameters
        ----------
        scenarionew : boolean
            True if we want a new scenario. False if we don't.
        timestart : string
            A string saying when to start the simulation time in STK format.
        duration : string
            A string saying how long to run the scenario in STK format.

        Returns
        -------
        None.

        """
        #Grab a handle on the STK application root.
        self.stkRoot = self.stk.Root
        if scenarionew:
            self.stkRoot.NewScenario('Constellation_Scenario')
            self.scenario = self.stkRoot.CurrentScenario
            self.scenario.SetTimePeriod(timestart, duration) #Usually 'Today' and '+48 hrs'

        else:
            self.scenario = self.stkRoot.CurrentScenario
            self.scenario.SetTimePeriod(timestart, duration)

    def generateConstellation(self,numOrbitPlanes,numSatsPerPlane,Orbits):
        self.constellation = self.scenario.Children.New(AgESTKObjectType.eConstellation, "SatConstellation")
        self.stkRoot.BeginUpdate()
        for orbitPlaneNum, RAAN in enumerate(range(0,180,180//numOrbitPlanes + 1),1): #RAAN in degrees

            for satNum, trueAnomaly in enumerate(range(0,360,360//numSatsPerPlane), 1): #trueAnomaly in degrees
                
                # Create Satellite
                satellite = self.scenario.Children.New(AgESTKObjectType.eSatellite, f"Sat{orbitPlaneNum}{satNum}")
                
                # Set Propagator
                satellite.SetPropagatorType(AgEVePropagatorType.ePropagatorAstrogator)
                self._driver = satellite.Propagator
                self._driver.MainSequence.RemoveAll()
                self._driver.Options.DrawTrajectoryIn3D = True

                # Define the initial state
                initState = self._driver.MainSequence.Insert(AgEVASegmentType.eVASegmentTypeInitialState, "Initial State", "-")
                initState.CoordSystemName = "CentralBody/Moon Inertial"
                
                # Modify parameters of initial state (keplerian)
                initState.SetElementType(AgEVAElementType.eVAElementTypeKeplerian)
                modKep = initState.Element

                modKep.PeriapsisRadiusSize = Orbits[orbitPlaneNum-1][1]   # Will change to for loop of array [2500:15000]
                modKep.ArgOfPeriapsis = Orbits[orbitPlaneNum-1][4]
                modKep.Eccentricity = Orbits[orbitPlaneNum-1][0]
                modKep.Inclination = Orbits[orbitPlaneNum-1][2]   # Will change to for loop of array [0:90]
                modKep.RAAN = Orbits[orbitPlaneNum-1][3]
                modKep.TrueAnomaly = trueAnomaly + (360//numSatsPerPlane/2)*(orbitPlaneNum%2) # Stagger true anomaly
                
                # Propagate with moon high precision propagator
                propagate = self._driver.MainSequence.Insert(AgEVASegmentType.eVASegmentTypePropagate, "Propagate", "-")
                propagate.Properties.DisplayCoordinateSystem = 'CentralBody/Moon Inertial'
                propagate.PropagatorName = "Moon HPOP Default v10"
                propagate.StoppingConditions["Duration"].Properties.Trip = 432000 #432000
                # Run MCS
                self._driver.RunMCS()
                # Add to constellation object
                self.constellation.Objects.AddObject(satellite)

        self.stkRoot.EndUpdate()
                
    def currentPositionKep(sat):
        #Give a sat you are interested in and return its Keplerian position.
        return [r.random(),r.randint(10000,100000),r.randint(1,10),r.randint(1,10),r.randint(1,10),r.randint(1,359)]
    
    def currentPositionCart(sat):
        #Give a sat you are interested in and return its Cartesian position.
        return [r.randint(10000,100000),r.randint(10000,100000),r.randint(10000,100000)]
    
    def lineOfSight(sat,locationofinterest):
        #Give a sat and location of interest and say whether there is a line of sight and give the distance if so.
        yesorno = r.randint(1,2)
        if yesorno==True:
            return True,r.randint(10000,100000) #Says that yes there is line of sight and returns the distance.
        else:
            return False,0 #Says no and returns a NaN value.

    def percentCoverage(self):
        """
        Calculates the percentage of the central body covered by the constellation.
        Returns
        -------
        Float
            The percentage of the surface covered by the constellation.
        """
        # Create coverage definition
        coverageDefinition = self.scenario.Children.New(AgESTKObjectType.eCoverageDefinition, "CoverageDefinition")
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
        return percentCoverage*100

    def clearScene(self, orbitPlaneNum, numSatsPerPlane):
        satellites = self.scenario.Children
        self.scenario.Children.Unload(AgESTKObjectType.eConstellation, "SatConstellation")
        for i in range(1, orbitPlaneNum+1):
            for j in range(1, numSatsPerPlane+1):
                satellites.Unload(AgESTKObjectType.eSatellite, f"sat{i}{j}")

orbitalelements=[[0,10000,75,20,0],[0,10000,75,80,0],[0,10000,75,140,0],[0,10000,75,200,0]]
propa = orbitalprop(False)
propa.CreateScenario(False,'Today','+120 hrs')
propa.generateConstellation(4,4,orbitalelements)
propa.clearScene(4, 4)
# print(propa.percentCoverage())
##stkRoot.EndUpdate()



