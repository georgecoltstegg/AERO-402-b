# Kevin Ott

import random
import numpy as np
# Number of planes
# Number of satellites per plane
# Antennas on satellites (Look at the code in Github, the properties you have to define are based on the antenna you choose. For example, if you choose parabolic you define aperture only, for helical you define length, diameter, number of turns, etc.)
# Uplink frequency
# Downlink frequency
# GPS frequency
# Max power on satellites
# Orbital information for each plane. (Eccentricity, Semimajor axis, Inclination, Longitude of ascending node, Argument of periapsis)

def LHS(nV, nS):
    # Initialization:
    k = 0
    # Creation of a list dictionary
    x = {}
    x[k] = []

    # Loop elements (part1)
    for i in range(nV):
        x1 = []

        for j in range(nS):
            a = ((float(j)) / nS)
            b = ((float(j) + 1) / nS)
            listesample = random.uniform(a, b)
            x1.append(listesample)

        # Select a random number nP times between each Sample and for each Var
        for k in range(nS):
            listechoice = random.choice(x1)
            x.setdefault(k, []).append(listechoice)
            x1.remove(listechoice)

    return list(x.values())

def mapLHS(birds, bounds, nS, nV):
    s = (nS, nV)
    pos = np.zeros(s)
    for mu in range(nS):
        for nu in range(nV):
            #Code for when first or second design variable is discrete
            if nu == 0:
                test = birds[mu][nu]
                if test > 0.5:
                    pos[mu][nu] = bounds[nu][1]
                else:
                    pos[mu][nu] = bounds[nu][0]
            elif nu == 1 or nu == 2:
                pos[mu][nu] = int(round(birds[mu][nu] * (bounds[nu][1]-bounds[nu][0]) + bounds[nu][0]))
            else:
                pos[mu][nu] = birds[mu][nu] * (bounds[nu][1]-bounds[nu][0]) + bounds[nu][0]
    return pos

#Binary DVs
argper = [90, 270]

#multi dimensional variables
parabolic = [0.1, 2.5] # diamater in meters
patch = [[0.04, 2.2], [0.04, 2.2]] #length, width both in meter
helical = [[0.1, 1.6], [2,64], [0.01, 1.28]] #diamter, number of turns, seperation
antennaType = [parabolic, patch, helical]
# Latin hyper square does not account for different antennas
# we are currently assuming customer and sat to have parabolic antennas in the bounds of 0.1 to 2.5 meters

#### begin Temp code
satAntenna = [0.1, 2.5]
groundAntenna = [0.1, 2.5]
#### end Temp code

sBand = [2000, 4000] #MHz
xBand = [8000, 12000] #MHz
kaBand = [27000, 40000] #MHz
upFreq = [sBand, xBand, kaBand]
downFreq = [sBand, xBand, kaBand]


#Discrete Variables
numPlanes = [1, 8]
numSatsPerPlane = [1, 8]

#Continous variables
gpsFreq = [500, 3000] #MHz
maxPower = [50, 400] #Watts
alt = [700, 30000] #km
ecc = [0, 0.05]
sma = [2000, 20000] #km
inc = [0, 180] #deg
raan = [0, 30] #deg

# creating different sets of bounds for LHS Mapping
designs = []
numV = 14 # number of variables
numS = 20 # number of seeds for lhs

# Calling LHS function
constellations = LHS(numV, numS)
count = 0
for i in range(3):
    for j in range(3):
        bounds = [argper, numPlanes, numSatsPerPlane, upFreq[i], downFreq[j], gpsFreq, maxPower, alt, ecc, sma, inc, raan, satAntenna, groundAntenna]
        designs.append(mapLHS(constellations, bounds, numS, numV))

# designs[0] -> up = sband, down = sband; designs[1]->up = sband, down = xband, designs[2]->up = sband, down=kaBand.. etc

# ECHO PRINT
print("Order of design variables in matrix:")
for i in [0,1,2,3,4,5,6,7,8,9,10,11,12,13]:
    text = ["Argument of perigee","Number of Planes", "Number of SatsPerPlane", "Uplink Frequency", "Downlink Frequency", "GNSS Freq", "Max Power", "Altitude", "Eccentricity", "Semi Major Axis", "Inclination Angle", "RAAN", "Satellite Antenna Diameter", "Ground Antenna Diameter"]
    print("\n", text[i])
    print(str(designs[0][0][i]) + ", " + str(designs[0][1][i]))

#DATA DUMP
print(designs[0])



