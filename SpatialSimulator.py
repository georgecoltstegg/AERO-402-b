# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 07:08:19 2022

@author: mtthwtrotter
"""

import os
import time
import random as r
from datetime import datetime

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
            STK_PID =13636 #os.getenv('STK_PID')19516
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
            self.stkRoot.Rewind()

    def generateConstellation(self,numOrbitPlanes,numSatsPerPlane,Orbits):
        self.constellation = self.scenario.Children.New(AgESTKObjectType.eConstellation, "SatConstellation")
        #self.stkRoot.BeginUpdate()
        for orbitPlaneNum, RAAN in enumerate(range(0,180,180//numOrbitPlanes + 1),1): #RAAN in degrees

            for satNum, trueAnomaly in enumerate(range(0,360,360//numSatsPerPlane), 1): #trueAnomaly in degrees
                
                # Create Satellite
                satellite = self.scenario.Children.New(AgESTKObjectType.eSatellite, f"Sat{orbitPlaneNum}{satNum}")
                
                # Set Propagator
                satellite.SetPropagatorType(AgEVePropagatorType.ePropagatorAstrogator)
                self._driver = satellite.Propagator
                self._driver.MainSequence.RemoveAll()
                self._driver.Options.DrawTrajectoryIn3D = False
                # Define the initial state
                initState = self._driver.MainSequence.Insert(AgEVASegmentType.eVASegmentTypeInitialState, "Initial State", "-")
                initState.CoordSystemName = "CentralBody/Moon Inertial"
                
                # Modify parameters of initial state (keplerian)
                initState.SetElementType(AgEVAElementType.eVAElementTypeKeplerian)
                initState.OrbitEpoch = self.scenario.StartTime
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
                propagate.StoppingConditions["Duration"].Properties.Trip = 54000 #86000 one day, 2,360,580 one lunar orbit around Earth, A circular orbit at 20km takes 71 hours, a 15k takes 45 hours.
                # Run MCS
                self._driver.RunMCS()
                # Add to constellation object
                self.constellation.Objects.AddObject(satellite)

        #self.stkRoot.EndUpdate()
                
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
    
    def satPositionsCart(self, sat):
        #Give a sats index you are interested in and return its Cartesian position.
        satellites = self.scenario.Children
        sat2 = satellites.GetElements(AgESTKObjectType.eSatellite)
        sat2 = sat2.Item(sat)
        times = sat2.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(self.scenario.StartTime,self.scenario.StopTime,60).DataSets.GetDataSetByName("Time").GetValues()
        xs = sat2.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(self.scenario.StartTime,self.scenario.StopTime,60).DataSets.GetDataSetByName("x").GetValues()
        ys = sat2.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(self.scenario.StartTime,self.scenario.StopTime,60).DataSets.GetDataSetByName("y").GetValues()
        zs = sat2.DataProviders.Item('Cartesian Position').Group.Item('ICRF').Exec(self.scenario.StartTime,self.scenario.StopTime,60).DataSets.GetDataSetByName("z").GetValues()
        return times,xs,ys,zs
    
    def facilityPositionsCart(self, fac):
        #Give a facilities index you are interested in and return its Cartesian position.
        facilities = self.scenario.Children
        facility = facilities.GetElements(AgESTKObjectType.eFacility)
        facility = facility.Item(fac)
        self.stkRoot.UnitPreferences.SetCurrentUnit("Time", "Min")
        facChooseDP = facility.DataProviders.Item('Points Choose System')
        dataProvCenter = facChooseDP.Group.Item('Center')
        # Choose the referense system you want to report the Center point in
        dataProvCenter.PreData = 'CentralBody/Moon Inertial'
        rptElems = [['Time'], ['x'], ['y'], ['z']]
        results = dataProvCenter.ExecElements(self.scenario.StartTime, self.scenario.StopTime, 60, rptElems)
        datasets = results.DataSets
        Time = datasets.GetDataSetByName('Time').GetValues()
        facTODx = datasets.GetDataSetByName('x').GetValues()
        facTODy = datasets.GetDataSetByName('y').GetValues()
        facTODz = datasets.GetDataSetByName('z').GetValues()
        """
        #times = sat2.DataProviders.Item('Cartesian Position').Exec().DataSets.GetDataSetByName("Time").GetValues()
        xs = sat2.DataProviders.Item('Cartesian Position').Exec().DataSets.GetDataSetByName("x").GetValues()
        ys = sat2.DataProviders.Item('Cartesian Position').Exec().DataSets.GetDataSetByName("y").GetValues()
        zs = sat2.DataProviders.Item('Cartesian Position').Exec().DataSets.GetDataSetByName("z").GetValues()
        """
        return Time,facTODx,facTODy,facTODz

    def generateFacility(self, locationofinterest, numname):
        facility = self.scenario.Children.New(AgESTKObjectType.eFacility, "L"+str(numname))
        facility.Position.AssignGeodetic(locationofinterest[0], locationofinterest[1], locationofinterest[2])
    def gatherScenarioElements(self):
        self.satellites_and_locations = self.scenario.Children
        self.sat2 = self.satellites_and_locations.GetElements(AgESTKObjectType.eSatellite)
        self.location2 = self.satellites_and_locations.GetElements(AgESTKObjectType.eFacility)
    def lineOfSight2(self,sat_num,location_num):
        sat = self.sat2.Item(sat_num)
        location = self.location2.Item(location_num)
        access = sat.GetAccessToObject(location)
        access.ComputeAccess()
        self.stkRoot.UnitPreferences.SetCurrentUnit("Time", "Min")
        accessDataProvider = access.DataProviders.GetDataPrvIntervalFromPath("Access Data")
        elements = ["Start Time", "Stop Time", "Duration"]
        accessResults = accessDataProvider.ExecElements(self.scenario.StartTime, self.scenario.StopTime, elements)
        startTimes = accessResults.DataSets.GetDataSetByName("Start Time").GetValues()
        stopTimes = accessResults.DataSets.GetDataSetByName("Stop Time").GetValues()
        durations = accessResults.DataSets.GetDataSetByName("Duration").GetValues()
        return startTimes[0],stopTimes[0],durations[0]
    def lineOfSight(self, sat,locationofinterest):
        #Give a sat and location of interest and say whether there is a line of sight and give the distance if so.
        facility = self.scenario.Children.New(AgESTKObjectType.eFacility, "KirbsHouse")
        facility.Position.AssignGeodetic(locationofinterest[0], locationofinterest[1], locationofinterest[2])
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
        return startTimes,stopTimes,durations

    def clearScene(self, orbitPlaneNum, numSatsPerPlane):
        satellites = self.scenario.Children
        self.scenario.Children.Unload(AgESTKObjectType.eConstellation, "SatConstellation")
        #self.scenario.Children.Unload(AgESTKObjectType.eCoverageDefinition, "CoverageDefinition")
        #self.scenario.Children.Unload(AgESTKObjectType.eFacility, "KirbsHouse")
        #self.scenario.Children.Unload(AgESTKObjectType.eChain, "Chain")
        for i in range(1, orbitPlaneNum+1):
            for j in range(1, numSatsPerPlane+1):
                satellites.Unload(AgESTKObjectType.eSatellite, f"sat{i}{j}")
    def STKDateParser(self, STKdate):
        STKdate=STKdate.split(' ')
        STKtime=STKdate[3].split(':')
        STKmicroseconds=STKtime[2].split('.')
        if STKdate[1] == 'Jan': STKdate[1]=1
        elif STKdate[1] == 'Feb': STKdate[1]=2
        elif STKdate[1] == 'Mar': STKdate[1]=3
        elif STKdate[1] == 'Apr': STKdate[1]=4
        return datetime(int(STKdate[2]),STKdate[1],int(STKdate[0]),int(STKtime[0]),int(STKtime[1]),int(STKmicroseconds[0]),0)

numlats=4
numlongs=4

FILE_NUM = 2 ######FILL THIS OUT FOR FILES

#GOOD ARCH:
#numplanes=[5]
#numsats=[4]
#orbitalelements=[[[.004,8025.9,39.53,-120,270],[.004,8025.9,39.53,-60,270],[.004,8025.9,39.53,0,270],[.004,8025.9,39.53,60,270],[.004,8025.9,39.53,120,270]]]
numplanes=[]
numsats=[]
orbitalelements=[]
with open('Architecture_Set'+str(FILE_NUM)+'.txt') as f:
    for k in range(59):
        print('Architecture loaded:'+str(k))
        number_of_planes = int(f.readline())
        number_of_sats = int(f.readline())
        string_orbits=[]
        for n in range(number_of_planes):
            string_orbits.append(f.readline().split(','))
        float_orbits=[]
        for n in range(number_of_planes):
            float_orbits.append([])
            for item in string_orbits[n]:
                float_orbits[-1].append(float(item))
        numplanes.append(number_of_planes)
        numsats.append(number_of_sats)
        orbitalelements.append(float_orbits)

"""
for n in range(len(numplanes)):
    print(numplanes[n])
    print(numsats[n])
    print(orbitalelements[n])
"""
propa = orbitalprop(False)
propa.CreateScenario(False,'1 Jan 2025 00:00:00.000','+15 hrs')
counter=0
for n in range(1,numlats):
    for i in range(1,numlongs):
        propa.generateFacility((-90+180*(n/numlats), -180+360*(i/numlongs),0), counter)
        counter+=1

for n in range(len(orbitalelements)):
    t0=time.time()
    print('Architecture '+str(n))
    propa.generateConstellation(numplanes[n],numsats[n],orbitalelements[n])
    #Find how the orbit has been perturbed.
    keps=[]
    for i in range(numplanes[n]*numsats[n]):
        currentkep = propa.currentPositionKep(i)
        keps.append(currentkep)
    #Gather the cartesian coordinates at different times for all of the satellites
    #Has coordinates (Satellite,Time)
    sats_times,sats_xs,sats_ys,sats_zs = [],[],[],[]
    for i in range(numplanes[n]*numsats[n]):
        sat_times,sat_xs,sat_ys,sat_zs = propa.satPositionsCart(i)
        sats_times.append(sat_times)
        sats_xs.append(sat_xs)
        sats_ys.append(sat_ys)  
        sats_zs.append(sat_zs)
    facs_times,facs_xs,facs_ys,facs_zs = [],[],[],[]
    for i in range((numlats-1)*(numlongs-1)):
        fac_times,fac_xs,fac_ys,fac_zs = propa.facilityPositionsCart(i)
        facs_times.append(fac_times)
        facs_xs.append(fac_xs)
        facs_ys.append(fac_ys)  
        facs_zs.append(fac_zs)
    #Write down when and different spots can see satellites
    sights = []
    #Load all the elements ahead of time for speed
    propa.gatherScenarioElements()
    #Check each point
    #This makes coordinates (Point,Satellite,(First sight, last sight, duration))
    for i in range((numlats-1)*(numlongs-1)):
        sights.append([])
        #Check each satellite
        for j in range(numplanes[n]*numsats[n]):
            try:
                sighted = propa.lineOfSight2(j,i)
                sights[-1].append([sighted[0],sighted[1],sighted[2]])
            except:
                sights[-1].append([])
            #print('Point '+str(i)+' can see satellite '+str(j)+' during these times:')
            #print(sighted)
            
    #We've gathered positions and times that are sighted, so now we check to see how many satellites are spotted and where their at at each time step.
    #Check to see if each satellite can see 4 satellites at all times
    print(sats_times[0][0])
    #Find the distance between points at time steps where they can be seen
    #The results are stored as (Time,Point,(PointCoordinate,SatelliteCoordinate1,SatelliteCoordinate2,...,SatelliteCoordinate3))
    all_locations = []
    #Check each time
    for k in range(len(sats_times[0])):
        #Check each point
        all_locations.append([])
        for i in range((numlats-1)*(numlongs-1)):
            all_locations[-1].append([])
            #Check each satellite
            all_locations[-1][-1].append([facs_xs[i][k],facs_ys[i][k],facs_zs[i][k]])
            for j in range(numplanes[n]*numsats[n]):
                print(k,i,j)
                if len(sights[i][j])==3:
                    isAfterOrEqualTo = (sights[i][j][0] < sats_times[0][k] or sights[i][j][0] == sats_times[0][k])
                    isBeforeOrEqualTo = (sights[i][j][1] > sats_times[0][k] or sights[i][j][1] == sats_times[0][k])
                    if isAfterOrEqualTo and isBeforeOrEqualTo:
                        all_locations[-1][-1].append([sats_xs[j][k],sats_ys[j][k],sats_zs[j][k]])
                    else:
                        pass
                else:
                    pass
    #Print initial architecture orbital elements
    f=open('ArchitectureResults'+str(FILE_NUM)+'-'+str(n)+'.txt','a')
    f.write(str(numplanes[n])+','+str(numsats[n])+'\n')
    for k in range(numplanes[n]):
        for i in range(5):
            f.write(str(orbitalelements[n][k][i])+',')
        f.write('\n')
    
    #Print final architecture for satellites
    for k in range(len(keps)):
        for i in range(6):
            f.write(str(keps[k][i])+',')
        f.write('\n')
    
    #Print positions of satellites relative to time and location
    for k in range(len(sats_times[0])):
        for i in range((numlats-1)*(numlongs-1)):
            for j in range(len(all_locations[k][i])):
                f.write(str(all_locations[k][i][j][0])+','+str(all_locations[k][i][j][1])+','+str(all_locations[k][i][j][2])+'\n')
            f.write('\n')
    f.close()
    propa.clearScene(numplanes[n],numsats[n])
    print('That arch took ',time.time()-t0)

