#2/1/2022
#Link Budget Calculator Portion for the AERO 402 CNCC Project
#This program estimates power necessary to maintain a certain bit error rate (BER)
#using the link budget equation. It also allows the user to create communication objects
#that store properties necessary for the link budget equation and to calculate dilution
#of precision, which this program also does.

import math as m
import numpy as np
class comobject():
    def __init__(self, antenna_type, antenna_size, location, max_power,pointing_error):
        self.antenna_type = antenna_type
        self.antenna_size = antenna_size
        self.location = location
        self.max_power = max_power
        self.pointing_error = pointing_error

    def calcGain(self,wl):
        #Returns the gain of the antenna in decibels. Needs the wavelength.
        #Parabolic needs a single float or integer representing diameter.
        #Patch needs a list or tuple with two elements representing the two sides of a rectangle.
        #Helix needs a list or tuple of three elements representing circumference, number of turns,
        #and vertical seperation in that order.
        if self.antenna_type == 'patch':
            if isinstance(self.antenna_size,list) or isinstance(self.antenna_size,tuple):
                if len(self.antenna_size)==2:
                    return 10*m.log10(.55*4*m.pi*self.antenna_size[0]*self.antenna_size[1]/wl**2)
                else:
                    print('Incorrect number of elements in list/tuple for antenna type patch.')
                    return 0
            else:
                print('Incorrect parameters loaded for antenna type patch.')
                return 0
        elif self.antenna_type == 'helix':
            if isinstance(self.antenna_size,list) or isinstance(self.antenna_size,tuple):
                if len(self.antenna_size)==3:
                    return 10*m.log10(6.2*self.antenna_size[0]**2*self.antenna_size[1]*self.antenna_size[2]/wl**3)
                else:
                    print('Incorrect number of elements in list/tuple for antenna type helix.')
                    return 0
            else:
                print('Incorrect parameters loaded for antenna type helix.')
                return 0
        elif self.antenna_type == 'parabolic':
            if isinstance(self.antenna_size,float) or isinstance(self.antenna_size,int):
                return 10*m.log10(((m.pi*self.antenna_size)/wl)**2*.55)
            else:
                print('Incorrect parameters loaded for antenna type parabolic.')
                return 0
        else:
            print('Invalid antenna type.')
            return 0

class dil_of_pre():
    def __init__(self):
        pass

    def calcR(self,satellite,interest):
        return m.sqrt((satellite.location[0]-interest[0])**2+(satellite.location[1]-interest[1])**2+(satellite.location[2]-interest[2])**2)

    def pseudorangeA(self,satellites,interest):
        arrays=[]
        for n in range(len(satellites)):
            array=[]
            for i in range(3):
                array.append((satellites[n].location[i]-interest[i])/self.calcR(satellites[n],interest))
            array.append(-1)
            arrays.append(array)
        return np.matrix(arrays)

    def covarianceQ(self,A):
        return np.linalg.inv(A.T*A)

    def calcDOP(self,satellitelist,interest):
        Q = self.covarianceQ(self.pseudorangeA(satellitelist,interest))
        return m.sqrt(Q[0,0]+Q[1,1]+Q[2,2])

class LinkBudg():
    def __init__(self):
        self.ebn0 = 16 #dB, for 10^-5 BER with BPSK encoding.
        self.lineloss = -1 #dB, assumed from SMAD
        self.c = 3 * 10 ** 8  # speed of light (m/s)

    def calcPower(self, transmitter, receiver, wl, data_rate):
        f = self.c / wl
        k = 1.38 * 10 ** -23  # boltzmann constant
        S=m.sqrt((transmitter.location[0]-receiver.location[0])**2+(transmitter.location[1]-receiver.location[1])**2+(transmitter.location[2]-receiver.location[2])**2)
        trans_point_loss = -12 * (transmitter.pointing_error) ** 2
        rec_point_loss = -12 * (receiver.pointing_error) ** 2
        space_loss = 20*m.log10(self.c)-20*m.log10(4*m.pi)-20*m.log10(S)-20*m.log10(f)
        L = -trans_point_loss * rec_point_loss * space_loss * self.lineloss  # losses
        Gt = transmitter.calcGain(wl)
        Gr = receiver.calcGain(wl)
        Rb = 20*m.log10(data_rate)
        Ts = 20*m.log10(135)  #Needs to be fixed later but is a valid approximation as long as we stay around .5 GHz to 10 GHz
        print(S)
        print(trans_point_loss,rec_point_loss,space_loss)
        print(self.ebn0,k,Ts,Rb,Gt,Gr,L)
        Pt = (self.ebn0 + k + Ts + Rb) / (Gt + Gr + L)
        print(Pt)
        return 10**(Pt/10) #Converts from dBW to W

satellite1 = comobject('parabolic',3.2,[8,2,1],100,1/6)
satellite2 = comobject('parabolic',3.2,[500000,1000000,60000000],100,1/6)
satellite3 = comobject('parabolic',3.2,[10,23,11],100,1/6)
dilscenario1=dil_of_pre()
#print(dilscenario1.calcDOP([satellite1,satellite2,satellite3],[100000,11000,3000000]))
LinkBudg1 = LinkBudg()
print(LinkBudg1.calcPower(satellite1,satellite2,0.1125,100000))


