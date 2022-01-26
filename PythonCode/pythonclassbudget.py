# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 04:07:44 2022

@author: Eric Comstock IV
"""

import numpy as np

class budget:
    def __init__(self, numsats = 1):
        self.numsats = numsats
        self.ctotal = 0#total cost
        self.csub = []#cost breakdown of system
        self.msat = 0#mass of 1 satellite
        self.msub = []#mass breakdown 1 satellite
        self.legend = ['comms', 'adcs', 'thermal', 'power', 'structure', 'propulsion', 'personeel']
    def populate_budget(self, costkg, mc, pc, madcs, padcs, theff, mth, pmpw, strfrac, deltavdeorbit, life, deltavskyr, ispfuel, numengr, salary, tdes, admin):
        self.msub.append(mc)#mass of comms
        self.msat += mc
        cc = 10 ** 3 * (20 + 230 * mc ** 0.59) + mc * costkg#cost of comms
        self.csub.append(self.numsats * cc)
        self.ctotal += self.numsats * cc
        
        self.msub.append(madcs)#mass of adcs
        self.msat += madcs
        cadcs = 10 ** 3 * np.min([madcs ** 0.73 * 186 - 364, 560.9]) + madcs * costkg#cost of adcs
        self.csub.append(self.numsats * cadcs)
        self.ctotal += self.numsats * cadcs
        
        self.msub.append(mth)#mass of thermal
        self.msat += mth
        cth = 10 ** 3 * (2640 + 416 * mth ** 0.66) + mth * costkg#cost of thermal
        self.csub.append(self.numsats * cth)
        self.ctotal += self.numsats * cth
        
        p = (pc + padcs) / theff#power budget calculated via thermal efficincy
        mp = p * pmpw#power mass per watt and power used to find power system mass
        self.msub.append(mp)#mass of power
        self.msat += mp
        cp = (183 * p ** 0.29) * 10 ** 3 + mp * costkg#cost of power
        self.csub.append(self.numsats * cp)
        self.ctotal += self.numsats * cp
        
        mstr = strfrac * (mc + madcs + mp + mth)#strfrac is structure to everything else ratio (not including propulsion - that gets its own structure)
        self.msub.append(mstr)#mass of power
        self.msat += mstr
        cstr = 10 ** 6 * (1.47 + np.log(mstr) * 0.07 * mstr) + mstr * costkg#cost of structure
        self.csub.append(self.numsats * cstr)
        self.ctotal += self.numsats * cstr
        
        mprop = (mc + madcs + mp + mth + mstr) * (np.exp((deltavdeorbit + life * deltavskyr) * 1.1 / ispfuel) - 1)
        self.msub.append(mprop)#mass of propulsion
        self.msat += mprop
        cprop = (0.704 + 0.0235 * (mprop / 6) ** 1.261) * 10 ** 6 + mprop * costkg#cost of propulsion
        self.csub.append(self.numsats * cprop)
        self.ctotal += self.numsats * cprop
        
        clabor = numengr * tdes * salary * admin#number of engineers times design time in years times engineer salary times admin. overhead
        self.csub.append(clabor + 0.0)
        self.ctotal += clabor
        self.msub.append(0)
    def print_budget(self):
        print('Total system cost (dollars):', self.ctotal)
        print('Number of satellites: ', self.numsats)
        print('Mass of each satellite (kg): ', self.msat)
        print()
        print('Cost Breakdown')
        print('    Component   Cost($)         Mass(kg)')
        for i in range(len(self.legend)):
            print('    ' + str(self.legend[i]).ljust(12, ' ') + str(round(self.csub[i], 2)).ljust(16, ' ') + str(round(self.msub[i], 3)).ljust(12, ' '))

sys = budget(16)
sys.populate_budget(500000, 1.61, 63, 6.53, 57.5, 1.2, 1.65, 0.01, 0.1, 300, 3, 150, 2500, 7, 150000, 1, 6)
sys.print_budget()