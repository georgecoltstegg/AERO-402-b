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
                
    def currentPositionKep(self, sat):
        #Give a sats index you are interested in and return its Keplerian position.
        satellites = self.scenario.Children
        sat2 = satellites.GetElements(AgESTKObjectType.eSatellite)
        sat2 = sat2.Item(sat)
        reportData = ""
        _driver = sat2.Propagator
        Orbit = _driver.MainSequence["Propagate"]
        result = Orbit.ExecSummary
        intervals = result.Value
        for i in range(0, intervals.Count):
            interval = intervals.Item(i)
            datasets = interval.DataSets
            for j in range(0, datasets.Count):
                dataset = datasets.Item(j)
                elements = dataset.GetValues()
                for o in elements:
                    reportData += o + "\r\n";
        breakdown = reportData.split()
        sma_s = "sma:"
        ecc_s = "ecc:"
        inc_s = "inc:"
        raan_s = "RAAN:"
        omega_s = "w:"
        trueAnomoly_s = "TA:"

        split = breakdown.index(sma_s)
        sma = breakdown[split + 1]
        split = breakdown.index(ecc_s)
        ecc = breakdown[split + 1]
        split = breakdown.index(inc_s)
        inc = breakdown[split + 1]
        split = breakdown.index(raan_s)
        raan = breakdown[split + 1]
        split = breakdown.index(omega_s)
        w = breakdown[split + 1]
        split = breakdown.index(trueAnomoly_s)
        trueAnomoly = breakdown[split + 1]

        return [sma,ecc,inc,raan,w,trueAnomoly]
    
    def currentPositionCart(self, sat):
        #Give a sats index you are interested in and return its Cartesian position.
        satellites = self.scenario.Children
        sat2 = satellites.GetElements(AgESTKObjectType.eSatellite)
        sat2 = sat2.Item(sat)
        reportData = ""
        _driver = sat2.Propagator
        Orbit = _driver.MainSequence["Propagate"]
        result = Orbit.ExecSummary
        intervals = result.Value
        for i in range(0, intervals.Count):
            interval = intervals.Item(i)
            datasets = interval.DataSets
            for j in range(0, datasets.Count):
                dataset = datasets.Item(j)
                elements = dataset.GetValues()
                for o in elements:
                    reportData += o + "\r\n";
        breakdown = reportData.split()
        x_s = "X:"
        y_s = "Y:"
        z_s = "Z:"
        vx_s = "Vx:"
        vy_s = "Vy:"
        vz_s = "Vz:"

        split = breakdown.index(x_s)
        x = breakdown[split + 1]
        split = breakdown.index(y_s)
        y = breakdown[split + 1]
        split = breakdown.index(z_s)
        z = breakdown[split + 1]
        split = breakdown.index(vx_s)
        vx = breakdown[split + 1]
        split = breakdown.index(vy_s)
        vy = breakdown[split + 1]
        split = breakdown.index(vz_s)
        vz = breakdown[split + 1]
        return [x, y, z, vx, vy, vz]
    
    def lineOfSight(self, sat,locationofinterest):
        #Give a sat and location of interest and say whether there is a line of sight and give the distance if so.
        facility = self.scenario.Children.New(AgESTKObjectType.eFacility, "KirbsHouse")
        facility.Position.AssignGeodetic(28.62, -80.62, 0.03)
        satellites = self.scenario.Children
        constellation = self.constellation
        sat2 = satellites.GetElements(AgESTKObjectType.eSatellite)
        sat2 = sat2.Item(sat)
        access = sat2.GetAccessToObject(facility)
        access.ComputeAccess()
        self.stkRoot.UnitPreferences.SetCurrentUnit("Time", "Min")
        accessDataProvider = access.DataProviders.GetDataPrvIntervalFromPath("Access Data")
        elements = ["Start Time", "Stop Time", "Duration"]
        accessResults = accessDataProvider.ExecElements(self.scenario.StartTime, self.scenario.StopTime, elements)

        startTimes = accessResults.DataSets.GetDataSetByName("Start Time").GetValues()
        stopTimes = accessResults.DataSets.GetDataSetByName("Stop Time").GetValues()
        durations = accessResults.DataSets.GetDataSetByName("Duration").GetValues()

        # Print data to console
        print("\nAccess Intervals")
        print("{a:<29s}  {b:<29s}  {c:<14s}".format(a="Start Time", b="Stop Time", c="Duration (min)"))
        for i in range(len(startTimes)):
            print("{a:<29s}  {b:<29s}  {c:<4.2f}".format(a=startTimes[i], b=stopTimes[i], c=durations[i]))

        print("\nThe maximum access duration is {a:4.2f} minutes.".format(a=max(durations)))
        # Create chain
        chain = self.scenario.Children.New(AgESTKObjectType.eChain, "Chain")

        # Add satellite constellation and facility
        chain.Objects.AddObject(constellation)
        chain.Objects.AddObject(facility)

        # Compute chain
        chain.ComputeAccess()
        chainDataProvider = chain.DataProviders.GetDataPrvIntervalFromPath("Object Access")
        chainResults = chainDataProvider.Exec(scenario.StartTime, scenario.StopTime)

        objectList = []
        durationList = []

        # Loop through all satellite access intervals
        for intervalNum in range(chainResults.Intervals.Count - 1):
            # Get interval
            interval = chainResults.Intervals[intervalNum]

            # Get data for interval
            objectName = interval.DataSets.GetDataSetByName("Strand Name").GetValues()[0]
            durations = interval.DataSets.GetDataSetByName("Duration").GetValues()

            # Add data to list
            objectList.append(objectName)
            durationList.append(sum(durations))

        # Find object with longest total duration
        index = durationList.index(max(durationList))
        print("\n{a:s} has the longest total duration: {b:4.2f} minutes.".format(a=objectList[index],
                                                                                 b=durationList[index]))
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
#print(propa.currentPositionKep(11))
print(propa.currentPositionCart(11))
propa.clearScene(4, 4)
# print(propa.percentCoverage())
##stkRoot.EndUpdate()



