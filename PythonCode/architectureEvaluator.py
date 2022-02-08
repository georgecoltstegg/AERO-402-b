# -*- coding: utf-8 -*-

# George Colts-Tegg

import time
import numpy as np
#import orbitalPropogator


class architecture:
    def __init__(self, antennaType, antennaSize, orbitalElems, maxPower, massAntennaType, massSat):
        self.antennaType = antennaType
        self.antennaSize = antennaSize
        self.orbitalElems = orbitalElems
        self.maxPower = maxPower
        self.massAntenna = massAntennaType
        self.massSat = massSat

    def orbitEval(self):
        # initialize the orbital elements
        inclination = self.orbitalElems[0]
        altitude = self.orbitalElems[1]
        numSats = self.orbitalElems[2]
        numPlanes = self.orbitalElems[3]
        eccentricity = self.orbitalElems[4]
        RAAN = self.orbitalElems[5]
        trueAnomaly = self.orbitalElems[6]
        argOfPerigee = self.orbitalElems[7]

        # insert more code/decisions after learning more about STK


    def antennaEval(self):
        if self.antennaType == 'patch':
            # Evaluate if the mass is within calculated parameters
            mass_antenna = 4*self.massAntenna[0]
            total_mass = mass_antenna + self.massSat

            power = 250 # random value
            if power > self.maxPower:
                print('The power required exceeds our max power')
                return 0
            else:
                print('The power required does not exceed max power for comms')

        elif self.antennaType == 'helix':
            mass_antenna = 4 * self.massAntenna[1]
            total_mass = mass_antenna + self.massSat

            power = 250 # random value
            if power > self.maxPower:
                print('The power required exceeds our max power')
                return 0
            else:
                print('The power required does not exceed max power for comms')

        elif self.antennaType == 'parabolic':
            mass_antenna = 4 * self.massAntenna[2]
            total_mass = mass_antenna + self.massSat

            power = 250
            if power > self.maxPower:
                print('The power required exceeds our max power')
            else:
                print('The power required does not exceed max power for comms')

        else:
            print('Invalid antennaType')


# antennaEval will go more in depth and orbitEval will as well once we learn more about STK, etc.
