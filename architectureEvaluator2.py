# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 11:59:08 2022

@author: mtthwtrotter
"""
import linkbudgetcalc as lbc
import AstrogatorConstellation as ac

class archEvaluator():
    def __init__(self,arch):
        self.arch = arch
    
    def printit(self):
        """
        Parameters
        ----------
        arch : List or tuple
            Contains all the parameters for the architecture.
            The tuple will have the following parameters (ellipses indicate arbitrary lengths)
            [Number of planes, Number of satellite per plane, ...Antenna types on satellites..., ...Antenna specifications
             on satellites..., Uplink frequency, Downlink frequency, GPS frequency, Max power for COMMs, Max power for GPS,
             ...Orbital information for each plane..., GPS Planes, COMM Planes]

        Returns
        -------
        int
            0 means unsuccessful architecture generation and 1 means successful power generation.
        """
        if isinstance(self.arch,list) or isinstance(self.arch,tuple):
            #Make sure everything is in order
            if isinstance(self.arch[0],int):
                pass
            else:
                print('Number of planes is not an integer.')
                print('No check will be made.')
                return 0
            if isinstance(self.arch[1],int):
                pass
            else:
                print('Number of satellites per plane is not an integer.')
                print('No check will be made.')
                return 0
            for i in range(2,self.arch[0]*self.arch[1]+2):
                if isinstance(self.arch[i],str):
                    if self.arch[i]=='parabolic' or self.arch[i]=='helix' or self.arch[i]=='patch':
                        pass
                    else:
                        print(str(self.arch[i])+' is not a valid antenna type. Please choose parabolic, helix, or patch.')
                        print('No check will be made.')
                        return 0
                else:
                    print('Antenna type is not a string.')
                    print('No check will be made.')
                    return 0
            for i in range(self.arch[0]*self.arch[1]+2,2*self.arch[0]*self.arch[1]+2):
                if self.arch[i-self.arch[0]*self.arch[1]] == 'parabolic':
                    if isinstance(self.arch[i],float) or isinstance(self.arch[i],int):
                        pass
                    else:
                        print('A parabolic antenna has a non-float or non-integer aperture size.')
                        print('No check will be made.')
                        return 0
                elif self.arch[i-self.arch[0]*self.arch[1]] == 'helix':
                    if isinstance(self.arch[i],list) or isinstance(self.arch[i],list):
                        if len(self.arch[i])==3:
                            for j in range(3):
                                if isinstance(self.arch[i][j],float) or isinstance(self.arch[i][j],int):
                                    pass
                                else:
                                    print('A helix antenna has a non-integer or non-float parameter specified.')
                                    print('No check will be made.')
                                    return 0
                        else:
                            print('A helix antenna has more or less than 3 parameters.')
                            print('No check will be made.')
                            return 0
                    else:
                        print('A helix antenna was not given a list or tuple of parameters.')
                        print('No check will be made.')
                        return 0
                else:
                    if isinstance(self.arch[i],list) or isinstance(self.arch[i],list):
                        if len(self.arch[i])==2:
                            for j in range(2):
                                if isinstance(self.arch[i][j],float) or isinstance(self.arch[i][j],float):
                                    pass
                                else:
                                    print('A patch antenna has a non-integer or non-float parameter specified.')
                                    print('No check will be made.')
                                    return 0
                        else:
                            print('A patch antenna has more or less than 2 parameters.')
                            print('No check will be made.')
                            return 0
                    else:
                        print('A patch antenna was not given a list or tuple of parameters.')
                        print('No check will be made.')
                        return 0
            if isinstance(self.arch[2*self.arch[0]*self.arch[1]+2],float) or isinstance(self.arch[2*self.arch[0]*self.arch[1]+2],int):
                pass
            else:
                print('Uplink frequency is not a float or an integer.')
                print('No check will be made.')
                return 0
            if isinstance(self.arch[2*self.arch[0]*self.arch[1]+3],float) or isinstance(self.arch[2*self.arch[0]*self.arch[1]+3],int):
                pass
            else:
                print('Downlink frequency is not a float or an integer.')
                print('No check will be made.')
                return 0
            if isinstance(self.arch[2*self.arch[0]*self.arch[1]+4],float) or isinstance(self.arch[2*self.arch[0]*self.arch[1]+4],int):
                pass
            else:
                print('GPS frequency is not a float or an integer.')
                print('No check will be made.')
                return 0
            if isinstance(self.arch[2*self.arch[0]*self.arch[1]+5],float) or isinstance(self.arch[2*self.arch[0]*self.arch[1]+5],int):
                pass
            else:
                print('Max power for COMMs is not a float or an integer.')
                print('No check will be made.')
                return 0
            if isinstance(self.arch[2*self.arch[0]*self.arch[1]+6],float) or isinstance(self.arch[2*self.arch[0]*self.arch[1]+6],int):
                pass
            else:
                print('Max power for GPS is not a float or an integer.')
                print('No check will be made.')
                return 0
            for i in range(2*self.arch[0]*self.arch[1]+7,2*self.arch[0]*self.arch[1]+7+self.arch[0]):
                if isinstance(self.arch[i],list) or isinstance(self.arch[i],tuple):
                    if len(self.arch[i])==5:
                        if isinstance(self.arch[i][0],float) or isinstance(self.arch[i][0],int):
                            if self.arch[i][0]>0 or self.arch[i][0]<1:
                                pass
                            else:
                                print('Eccentricity is not between 0 and 1.')
                                print('No check will be made.')
                                return 0
                        else:
                            print('Eccentricity is not a float or an integer.')
                            print('No check will be made.')
                            return 0
                        if isinstance(self.arch[i][1],float) or isinstance(self.arch[i][1],int):
                            if self.arch[i][1] < 1737400/1000:
                                print('The plane is within the moons radius.')
                                print('No check will be made.')
                                return 0
                            else:
                                pass
                        else:
                            print('Semi-major axis is not a float or an int.')
                            print('No check will be made.')
                            return 0
                        if isinstance(self.arch[i][2],float) or isinstance(self.arch[i][2],int):
                            pass
                        else:
                            print('Inclination is not a float or an int.')
                            print('No check will be made.')
                            return 0
                        if isinstance(self.arch[i][3],float) or isinstance(self.arch[i][3],int):
                            pass
                        else:
                            print('Longitude of the ascending node is not a float or an int.')
                            print('No check will be made.')
                            return 0
                        if isinstance(self.arch[i][4],float) or isinstance(self.arch[i][4],int):
                            pass
                        else:
                            print('Argument of periapsis is not a float or an int.')
                            print('No check will be made.')
                            return 0
                    else:
                        print('More or less than 5 orbital parameters were given for a plane.')
                        print('No check will be made.')
                        return 0
                else:
                    print('Orbital parameters were not given in a list or a tuple.')
                    print('No check will be made.')
                    return 0
            if isinstance(self.arch[2*self.arch[0]*self.arch[1]+7+self.arch[0]],list) or isinstance(self.arch[2*self.arch[0]*self.arch[1]+7+self.arch[0]],tuple):
                if len(self.arch[2*self.arch[0]*self.arch[1]+7+self.arch[0]])<=self.arch[0]:
                    for k in range(len(self.arch[2*self.arch[0]*self.arch[1]+7+self.arch[0]])):
                        if isinstance(self.arch[2*self.arch[0]*self.arch[1]+7+self.arch[0]][k],int):
                            pass
                        else:
                            print('A plane given for GPS is not an integer.')
                            print('No check will be made.')
                            return 0
                        if -1<self.arch[2*self.arch[0]*self.arch[1]+7+self.arch[0]][k]<=self.arch[0]:
                            pass
                        else:
                            print('A plane given for GPS is less than 0 or outside the bounds of the number of planes.')
                            print('No check will be made.')
                            return 0
                else:
                    print('More planes given than planes that exist for GPS.')
                    print('No check will be made.')
                    return 0
            else:
                print('A list is not given for the planes that will do GPS.')
                print('No check will be made.')
                return 0
            if isinstance(self.arch[2*self.arch[0]*self.arch[1]+8+self.arch[0]],list) or isinstance(self.arch[2*self.arch[0]*self.arch[1]+8+self.arch[0]],tuple):
                if len(self.arch[2*self.arch[0]*self.arch[1]+8+self.arch[0]])<=self.arch[0]:
                    for k in range(len(self.arch[2*self.arch[0]*self.arch[1]+8+self.arch[0]])):
                        if isinstance(self.arch[2*self.arch[0]*self.arch[1]+8+self.arch[0]][k],int):
                            pass
                        else:
                            print('A plane given for GPS is not an integer.')
                            print('No check will be made.')
                            return 0
                        if -1<self.arch[2*self.arch[0]*self.arch[1]+8+self.arch[0]][k]<=self.arch[0]:
                            pass
                        else:
                            print('A plane given for COMMs is less than 0 or outside the bounds of the number of planes.')
                            print('No check will be made.')
                            return 0
                else:
                    print('More planes given than planes that exist for COMMs.')
                    print('No check will be made.')
                    return 0
            else:
                print('A list is not given for the planes that will do COMMs.')
                print('No check will be made.')
                return 0
            print('The parameters pass the check.')
            print('This architecture has:')
            print(str(self.arch[0])+' planes.')
            print(str(self.arch[1])+' satellites per plane.')
            for i in range(self.arch[0]):
                for j in range(self.arch[1]):
                    print('Satellite '+str(j)+' in plane '+str(i)+' has a '+self.arch[2+self.arch[1]*i+j]+' type antenna which has parameters '+str(self.arch[2+self.arch[0]*self.arch[1]+self.arch[1]*i+j]))
            print(self.arch[2*self.arch[0]*self.arch[1]+2],' GHz Uplink Frequency')
            print(self.arch[2*self.arch[0]*self.arch[1]+3],' GHz Downlink Frequency')
            print(self.arch[2*self.arch[0]*self.arch[1]+4],' GHz GPS Frequency')
            print(self.arch[2*self.arch[0]*self.arch[1]+5],' W Max COMM power for each satllite')
            print(self.arch[2*self.arch[0]*self.arch[1]+6],' W Max GPS power for each satllite')
            for i in range(self.arch[0]):
                print('Plane '+str(i)+' has orbital parameters '+str(self.arch[2*self.arch[0]*self.arch[1]+7+i]))
            print('Planes '+str(self.arch[2*self.arch[0]*self.arch[1]+7+self.arch[0]])+' have GPS capabilities.')
            print('Planes '+str(self.arch[2*self.arch[0]*self.arch[1]+8+self.arch[0]])+' have COMM capabilities.')
            return 1
        else:
            print("Architecture decision checker's printit was passed a variable that was not a list or a tuple.")
            print('No check will be made.')
            return 0
        
    def generateComObjects(self):
        sats=[]
        self.all_powers=[]
        for n in range(self.arch[0]):
            power=0
            if n in self.arch[-2]:
                power += self.arch[2*self.arch[0]*self.arch[1]+6]
            if n in self.arch[-1]:
                power += self.arch[2*self.arch[0]*self.arch[1]+5]
            for i in range(self.arch[1]):
                sats.append(lbc.comobject(self.arch[2+self.arch[1]*n+i],self.arch[2+self.arch[0]*self.arch[1]+self.arch[1]*n+i],[0,0,0],power,1/12)) #1/12 is the assumed pointing error.
            self.all_powers.append(power)
        return sats
    
    def totalPowerofConstellation(self):
        return self.all_powers
    
    def generateOrbitalPropInstance(self):
        propa = ac.orbitalprop(True)
        propa.CreateScenario(True,'Today','+120 hrs')
        orbitalelements=[]
        for i in range(self.arch[0]):
            orbitalelements.append(self.arch[2*self.arch[0]*self.arch[1]+7+i])
        propa.generateConstellation(self.arch[0],self.arch[1],orbitalelements)
        return propa
    
    def NetworkCostTuple(self):
        return [500000, 1.61, 63, 9.53, 57.5, 0.8, 1.65, 0.01, 0.1, 3, 150, 2500]
    
#currentarch = archEvaluator([4,1,'parabolic','parabolic','helix','parabolic',3.2,3.4,[1,2,3],2.8,10,20,30,300,200,[0,10000,75,20,0],[0,10000,75,80,0],[0,10000,75,140,0],[0,10000,75,200,0],[0,1,2],[0,1,2,3]])
#currentarch.printit()
#sats_com_objects = currentarch.generateComObjects()
#propa = currentarch.generateOrbitalPropInstance()
#print(currentarch.totalPowerofConstellation())