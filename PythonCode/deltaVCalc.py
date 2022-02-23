from math import sqrt, sin, cos, atan, asin

def deltaVCalculation(orbitalelements, currentPositionKep,numSatsPerPlane,orbitPlaneNum):
    #satellite 1 desired orbital element
    a1, e1, i1, right_ascension1, f1 = [orbitalelements[0][1], 1.*pow(10.,-15.), orbitalelements[0][2], orbitalelements[0][3], (360//numSatsPerPlane/2)*(orbitPlaneNum%2)]
    G = 6.67408* pow(10.,-11.) #m3 kg-1 s-2
    M_moon = 7.34767309 *pow(10.,22.) #kg
    r1 = a1
    h1 = sqrt(2*G*M_moon) * sqrt((r1*r1)/(r1+r1)) #circular orbit thus ra = rb = r

    #pulling satellites 1 current oribtal elements after 1 day
    a2,e2,i2,right_ascension2,h2,f2 = currentPositionKep

    i = i2                              # inclination
    right_ascension = right_ascension2  # right ascension
    e = e2                              # eccentricity
    a = a2                              # semi-major axis
    f = f2                              # true anomaly
    p = a*(1-pow(e,2.))                 # semilatus rectum
    h = h2                              # orbit angular momentum
    r = a2                              # scalar orbit radius
    ascending_node = 0.                 # ascending node angle
    orbit_normal_impulse = 0.  # magnitude of orbit normal impulsive delta V, used to adjust inclination angle and
    theta = right_ascension + f         # true latitude angle
    nu = sqrt(1. - pow(e, 2.))          # another convenient form of eccentricity measurement
    omega = 0.                          # argument of perigee
    E1 = atan((sqrt(1.-pow(e1,2.))*sin(f1))/(e1+cos(f1)))
    E2 = atan((sqrt(1.-pow(e2,2.))*sin(f2))/(e2+cos(f2)))
    M1 = E1 - e*sin(E1)                 # mean anomaly
    M2 = E2 - e*sin(E2)                 # mean anomaly

    delta_right_ascension = right_ascension2 - right_ascension1  # desired change to right ascension measure
    delta_i = i2-i1     # desired change to inclination angle
    delta_omega = 0.    # change in argument of perigee
    delta_M = M2-M1     # change in mean anomaly

    na = h / p * nu
    theta_c = atan(delta_right_ascension * sin(i) / delta_i)  # critical true latitude angle at which to perform orbit normal

    #calculate delta V for each orbital element correction
    magnitude_orbit_normal_impulse = h / r * sqrt(pow(delta_i, 2.) + pow(delta_right_ascension, 2.) * sin(pow(i, 2.)))
    magnitude_delta_radial_impulse_perigee = -(na / 4) * (
                (pow((1. + e), 2.) / nu) * (delta_omega + delta_right_ascension * cos(i)) + delta_M)
    magnitude_delta_radial_impulse_apogee = (na / 4) * (
                (pow((1. - e), 2.) / nu) * (delta_omega + delta_right_ascension * cos(i)) + delta_M)

    delta_a = a2-a1     # desired change in apogee
    delta_e = e2-e1     # desired change in eccentricity
    na = h / (a * nu)
    magnitude_delta_tangential_impulse_perigee = (na * nu / 4) * (delta_a / a + delta_e / (1 + e))
    magnitude_delta_tangential_impulse_apogee = (na * nu / 4) * (delta_a / a - delta_e / (1 - e))

    magnitude_delta = ((abs(magnitude_orbit_normal_impulse) + abs(magnitude_delta_radial_impulse_perigee)
         +abs(magnitude_delta_radial_impulse_apogee) + abs(magnitude_delta_tangential_impulse_perigee) + 
         abs(magnitude_delta_tangential_impulse_apogee))*365.*3.*1000.) #m/s
    print("3 Year delta V: " +  str(magnitude_delta) + " m/s")

#TODO: implement burns at certain times during a full 3 year simulation
# adjusted best when the spacecraft passes through either the polar or the equatorial regions
def ascendingNode_inclinationAngle_Change(a, e, h, r, i, omega, f, delta_i, delta_right_ascension, delta_a, delta_e,
                                          J2):
    nu = sqrt(1. - pow(e, 2.))
    magnitude_orbit_normal_impulse = h / r * sqrt(pow(delta_i, 2.) + pow(delta_right_ascension, 2.) * sin(pow(i, 2.)))
    theta_c = atan(delta_right_ascension * sin(i) / delta_i)
    print("Delta V normal: " + magnitude_orbit_normal_impulse)
    print("Critical True Latitude angle: " + theta_c)

    # changes in inclination angle affect semi-major axis and eccentricity thus must be accounted for
    re = 1.  # radius of Earth
    delta_a += (-3. / 2. * J2 * sin(2. * i) * pow(re, 2.) / a * (
                1 / pow(nu, 3.) + pow((a / r), 3.) * (cos(2. * omega + 2. * f) - 1.))) * delta_i
    delta_e += (-J2 / 4. * pow(re, 2.) / pow(a, 2.) * sin(2. * i) / pow(nu, 2.) *
                (e / 4. * (11. + (80. * pow(cos(i), 2.)) / (1. - 5. * pow(cos(i), 2.)) + (200. * pow(cos(i), 4.)) / pow(
                    (1. - 5. * pow(cos(i), 2.)), 2.)) * cos(2. * omega)
                 + 3. / (e * pow(nu, 4.)) * ((pow(a / r, 3.) - 1. / pow(nu, 3.)) + (pow(a / r, 3.) - 1. / pow(nu, 4.)))
                 - 3. * cos(2. * omega + f) - cos(2. * omega + 3. * f))) * delta_i


def argumentOfPerigee_meanAnomaly_Change(h, p, nu, delta_omega, delta_M, delta_right_ascension):
    na = h / p * nu
    magnitude_delta_radial_impulse_perigee = -(na / 4) * (
            (pow((1. + e), 2.) / nu) * (delta_omega + delta_right_ascension * cos(i)) + delta_M)
    magnitude_delta_radial_impulse_apogee = (na / 4) * (
            (pow((1. - e), 2.) / nu) * (delta_omega + delta_right_ascension * cos(i)) + delta_M)
    print("Delta V radial:")
    print("Delta Vr perigee: " + magnitude_delta_radial_impulse_perigee)
    print("Delta Vr apogee: " + magnitude_delta_radial_impulse_apogee)


def semiMajorAxis_eccentricity_Change(h, a, e, delta_a, delta_e):
    nu = sqrt(1. - pow(e, 2.))
    na = h / (a * nu)
    magnitude_delta_tangential_impulse_perigee = (na * nu / 4.) * (delta_a / a + delta_e / (1. + e))
    magnitude_delta_tangential_impulse_apogee = (na * nu / 4.) * (delta_a / a - delta_e / (1. - e))

orbitalelements=[[0,10000,75,20,0],[0,10000,75,80,0],[0,10000,75,140,0],[0,10000,75,200,0]]
currentPositionKep = [10009.3357330703711341, 0.0019435337176272, 74.65442301229011, 19.97660004412605, 71.5193249662165, 318.4903359246459]
numSatsPerPlane = 4
orbitPlaneNum = 4
deltaVCalculation(orbitalelements,currentPositionKep, numSatsPerPlane, orbitPlaneNum)


