##########################################################
##
## SNvisMap_int

##
## sky map of Milky Way supernova probability per solid angle

##
## calculated for a grid of points in (l,b) Galactic coordinates
## and plotted as a countour map

##
## method of calculation:

##   given supernova rate density q(R,z) in Galactocentric cylindrical coords

##   total rate is:  R_SN = int q dV

##   rate per solid angle is:  dR/dOmega = int q dV/dOmega = int q r^2 dr

##   and so probability density = probability per solid ange is:

##      dP/dOmega = dR/dOmega / R_SN

##
## Inputs:

##   SN_type:  core collapse, Type Ia

##   SN_dist:  supernova rate density distribution

##   zoom:  show smaller longitude field

##   N_l1q, N_b1q:  number of latitude (l) and longitude (b) points;
##                  this sets the resolution of the map and the runtime

##
## Outputs:

##   contour map of probability per solid angle

############################################################


import numpy as np

import matplotlib.pyplot as plt

import matplotlib.colorbar as colorbar

import scipy.integrate as integrate

import matplotlib.ticker as ticker

from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter

from math import gamma

import time


##### start the clock
t0 = time.time()


### welcome messages
print ("Sky map of MW SN probability")

#print ("imposes reflection symmetry in Galactic latitude and longitude")
#print( "i.e., calculates for one quadrant and copies to the others")


## model parameters

## supernova density distribution model

##
## supernova type

PlotType_Ia = 'Ia_fiducial'
PlotType_CC = 'CC_fiducial'
PlotType_CC26Al = 'CC26Al'
PlotType_Ia26Al = 'Ia26Al'
PlotType_CCbigR = 'CC_bigR'
PlotType_IabigR = 'Ia_bigR'

PlotType = PlotType_CCbigR
PlotType = PlotType_IabigR
PlotType = PlotType_Ia26Al
PlotType = PlotType_CC
PlotType = PlotType_Ia
PlotType = PlotType_CC26Al

#SN_dist = 'spherical'
#SN_dist = 'thindisk'
#SN_dist = 'Adams'
#SN_dist = 'Green'
#SN_dist = '60Fe'
#SN_dist = '26Al'


if (PlotType == PlotType_CC):
    SN_type = 'CC'
    SN_dist = 'CC'
    SN_label = "Core-Collapse"
elif (PlotType == PlotType_Ia):
    SN_type = 'Ia'
    SN_dist = 'Ia'
    SN_label = "Type Ia"
elif (PlotType == PlotType_CC26Al):
    SN_type = 'CC'
    SN_dist = '26Al'
    SN_label = "Core-Collapse"
elif (PlotType == PlotType_Ia26Al):
    SN_type = 'Ia'
    SN_dist = '26Al'
    SN_label = "Ia"    
else:
    print ("Bad plot type: %s" % (PlotType))



## plot region is zoomed in longitidue

#ZoomLong = True

ZoomLong = False

full_lat = False
full_lat = True


## resolution: longitude and latitude points in one quadrant

N_l = 180+1
N_b = 15+1
N_l = 2*360+1
N_b = 2*30+1
N_b = 2*30+1
N_l = 2*45+1
N_l = 2*45+1
N_b = 2*90+1
N_l = 1*45+1
N_b = 1*30+1
N_l = 4*45+1
N_b = 4*90+1



def dP_drdldb(r,bee,ell):
    # integrand for probability integral

    global R_sun, z_sun

    r_upper = 6.*R_sun
    V_tot = 4.*np.pi*(r_upper**3)/3.

    z = r * np.sin(bee) + z_sun

    r_par = r * np.cos(bee)

    R = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(ell) + R_sun**2)


    integrand = np.cos(bee) * r**2 * q_SN(R,z)
    #integrand_test = np.cos(bee) * r**2 / V_tot
    
    return integrand


def dPvis_drdldb(r,bee,ell):
    # integrand for probability integral

    global R_sun, z_sun

    r_upper = 6.*R_sun
    V_tot = 4.*np.pi*(r_upper**3)/3.

    z = r * np.sin(bee) + z_sun

    r_par = r * np.cos(bee)

    R = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(ell) + R_sun**2)


    integrand = np.cos(bee) * r**2 * q_SN(R,z)
    #integrand_test = np.cos(bee) * r**2 / V_tot
    
    return integrand




def dPsn_drdOmega(radius):

    # integrand of supernova probability per solid angle on sky

    # dN_sn/dr dOmega = r^2 q_sn(R,z) dr

    # where R = R(r,l,b) and z = z(r,b)

    global l_rad, b_rad

    global R_sun, z_sun

    z = radius * np.sin(b_rad) + z_sun

    r_par = radius * np.cos(b_rad)

    R = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(l_rad) + R_sun**2)


    dP_drdOmega = radius**2 * q_SN(R,z)


    return dP_drdOmega



def dPsn_drdl(radius):

    

    # integrand of supernova probability per solid angle on sky

    # dN_sn/dr dOmega = r^2 q_sn(R,z) dr

    # where R = R(r,l,b) and z = z(r,b)

    global l_rad

    global R_sun, z_sun

    z = radius * np.sin(b_rad) + z_sun

    r_par = radius * np.cos(b_rad)

    R = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(l_rad) + R_sun**2)


    dP_drdOmega = radius**2 * q_SN(R,z)


    return dP_drdOmega




def q_SN(R_gc,z_gc):

    # supernova density in Galactocentric cylindrical coordinates

    global R_disk, h_disk
    global R_thin, h_thin
    global R_thick, h_thin
    global R_26Al, h_26Al
    global R_60Fe, h_60Fe

    global SN_dist


    if (SN_dist == 'Adams'):

        q_0 = 1./(4.*np.pi*R_disk**2*h_disk)

        q = q_0 * np.exp(-R_gc/R_disk) * np.exp(-np.abs(z_gc)/h_disk)

    elif (SN_dist == 'spherical'):

        q_0 = 1./(8.*np.pi*R_disk**3)

        r_gc = np.sqrt(R_gc**2 + z_gc**2)

        q = np.exp(-r_gc/R_disk)

    elif (SN_dist == 'thindisk'):

        q_0 = 1./(8.*np.pi*R_disk*3)

        r_gc = np.sqrt(R_gc**2 + z_gc**2)

        q = np.exp(-r_gc/R_disk)

    elif (SN_dist == 'Green'):
        alpha = 1.09
        beta = 3.87
        R_sun = 8.5 # kpc
        R_sn = 295
        R_0 = R_sun / beta # kpc
        h = 0.095 # kpc
        r_gc = np.sqrt(R_gc**2 + z_gc**2)
        q = (R_sn/(4*np.pi*gamma(alpha+2)*beta**alpha*R_0**2*h))*((abs(r_gc)/R_0)**alpha)*np.exp(-abs(r_gc)/R_0)*np.exp(-abs(z_gc)/h)


    elif (SN_dist == 'Ia'):

        h_1 = h_thin
        R_1 = R_thin
        p_1 = 0.5

        h_2 = h_thick
        R_2 = R_thick
        p_2 = 1. - p_1

        #Q_int = 8.*np.pi * (0.5 * h_1 * R_1**2 + 0.5 * h_2 * R_2**2)
        #norm = 1/Q_int

        norm_1 = 1./(4.*np.pi * h_1 * R_1**2) 
        q_1 = norm_1 * np.exp( -R_gc/R_1 ) * np.exp( -np.abs(z_gc)/h_1)

        norm_2 = 1./(4.*np.pi * h_2 * R_2**2) 
        q_2 = norm_2 * np.exp( -R_gc/R_2 ) * np.exp( -np.abs(z_gc)/h_2)

        q = p_1*q_1 + p_2*q_2


    elif (SN_dist == 'CC'):

        h_1 = h_thin
        R_1 = R_thin

        q = (1./(4.*np.pi * R_1**2 * h_1)) * np.exp( -R_gc/R_1 ) * np.exp( -np.abs(z_gc)/h_1)


    elif (SN_dist == '26Al'):

        h_1 = h_26Al
        R_1 = R_26Al

        q = (1./(4.*np.pi * R_1**2 * h_1)) * np.exp( -R_gc/R_1 ) * np.exp( -np.abs(z_gc)/h_1)


    elif (SN_dist == '60Fe'):

        h_1 = h_60Fe
        R_1 = R_60Fe

        q = (1./(4.*np.pi * R_1**2 * h_1)) * np.exp( -R_gc/R_1 ) * np.exp( -np.abs(z_gc)/h_1)


    else:

        q = 0


    return q



def Psn_int(l_min,l_max,b_min,b_max):

    global R_sun

    r_min = 0.
    r_max = 6.*R_sun

    b_max_LSST = (10./90.)*np.pi/2.

    def b_inf(l):
        return b_min
    def b_sup(l):
        return b_max
    
    def r_inf(l,b):
        return r_min
    def r_sup(l,b):
        return FindDist(m_lim,M_SN,l,b)
    
    #P_1q, err = integrate.tplquad(dP_drdldb,l_min,l_max,b_inf,b_sup,r_inf,r_sup)
    #P_sn = 4. * P_1q
    P_sn, err = integrate.tplquad(dPvis_drdldb,l_min,l_max,b_inf,b_sup,r_inf,r_sup)

    print ("P and err", P_sn, err)

    return P_sn



def rho_dust(R,z):

    #"Dust Density Function"

    #R_thin = 2.9  # kpc

    #h_thin = 0.095 # kpc


    rho = np.exp(-R/R_thin)*np.exp(-np.abs(z)/h_thin)

    rho = rho / (R_thin * (1. - np.exp(-R_sun/R_thin)))

    return rho


def dAv_dr(radius, l_rad, b_rad):
#"Extinction Rate due to Dust"
    global R_sun, z_sun

    z_cyl = z_sun + radius * np.sin(b_rad) # solarcentric radius component in plane

    r_par = radius * np.cos(b_rad) # Galactocentric cylindrical radius
    
    R_cyl = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(l_rad) + R_sun**2)



    Av_gc = 30.0

    dAv_dr = Av_gc * rho_dust(R_cyl,z_cyl)


    return dAv_dr


def A_ext(l_rad,b_rad,radi):

#"Total dimming due to dust"

    Sigfunct = lambda r: dAv_dr(r, l_rad, b_rad)
    #so dAv_dr will intgrate correctly
    
    Sigma, err = integrate.quad(Sigfunct, 0. ,radi)
    #get magnitude loss due to extinction

    return Sigma



def m_app(r0,l0,b0,M_abs):
    d_0 = 0.010 # kpc
    mu = 5.*np.log10(r0/d_0)
    A_los = A_ext(l0,b0,r0)
    m_los = M_abs + mu + A_los
    return m_los


def dm_dr(r0,l0,b0):
    dmu_dr = 5./(np.log(10.)*r0)
    dA_dr = dAv_dr(r0,l0,b0)
    dmdr = dmu_dr + dA_dr
    return dmdr
    

def FindDist(m_lim,M_SN,ll,bb):
    # find distance to supernova
    global R_sun
    global N_obscure

    r_max = 6.*R_sun
    m_max = m_app(r_max,ll,bb,M_SN)

    r_old = r_max
    m_old = m_max

    eps = 0.01

    count = 0

    if (m_old > m_lim):

        N_obscure = N_obscure+1

        while (np.abs(m_old-m_lim) > eps):
            # use Newton's method
            count = count+1
            r_new = r_old - (m_app(r_old,ll,bb,M_SN)-m_lim)/dm_dr(r_old,ll,bb)
            if (r_new < 0.):
                r_new = r_old/2.
            m_new = m_app(r_new,ll,bb,M_SN)
            #print "%i: (l,b)=(%.2f,%.2f)deg, and (%.2f kpc, %.2f mag) -> (%.2f kpc, %.2f mag)" \
            #    % (count, ll*180./np.pi, bb*180./np.pi, r_old, m_old, r_new, m_new)
            r_old = r_new
            m_old = m_new
            
    return r_old
        

########################################
# geometry parameters

R_sun = 8.7  # kpc
z_sun = 0.000 # kpc
z_sun = 0.020 # kpc

R_cc = 2.9  # kpc
R_ia = 2.4 # kpc
h_cc = 0.05 # kpc
h_ia = 0.8 # kpc

# TRILEGAL
h_dust = 0.110 # kpc

h_thin = 0.095 # kpc
R_thin = 2.9 # kpc

h_thick = 0.800 # kpc
R_thick = 2.4 # kpc


# Gamma Lines
# from radioactive gamma lines:
# Wang+ 2020 https://ui.adsabs.harvard.edu/abs/2020ApJ...889..169W/abstract
# double exponential fit to 26Al and 60Fe gamma line maps
## 26Al is more reliable
R_26Al = 7.0 # kpc   error:  +1.5 -1.0
h_26Al = 0.8 # kpc   error:  +0.3 -0.2
## 60Fe less reliable
R_60Fe = 3.5 # kpc   error:   +2.0 -1.5
h_60Fe = 0.3 # kpc   error:   +2.0 -0.2

########################################


if (ZoomLong):

    l_max_deg = 90.

else:

    l_max_deg = 180.



if (SN_type == 'Ia'):

    h_disk = h_ia

    R_disk = R_ia
 # kpc

    b_max_deg = 15.
    if (full_lat):
        b_max_deg = 90.

    labtext = "Type Ia Supernovae"

    if (SN_dist == 'thindisk'):

        h_disk = 0.350 # kpc
elif (SN_type == 'CC'):

    h_disk = h_cc

    R_disk = R_cc

    b_max_deg = 15.
    if (full_lat):
        b_max_deg = 90.

    labtext = "Core-Collapse Supernovae"

    if (SN_dist == 'thindisk'):

        h_disk = 0.350 # kpc

else:

    print ("bad SN type!")




##### calcuation begins
basename = "SNProb_"
zsunlabel = "zsun%.0f" % (z_sun*1.e3)
thindisklabel = "Rthin%.1f_hthin%.0f" % (R_thin,h_thin*1.e3)
thickdisklabel = "Rthick%.1f_hthick%.0f" % (R_thick,h_thick*1.e3)
reslabel = "res%ix%i" % (N_l,N_b)
figbasename = basename + SN_type + "_" + SN_dist + "_" + zsunlabel + "_" + thindisklabel + "_" + thickdisklabel + "_" + reslabel


dataname_summary = figbasename + ".dat"
DataFile = open(dataname_summary,"w")

dataname_points = figbasename + "_prob.dat"
DataFile_points = open(dataname_points,"w")


print ("here goes P_SN")

CalcTotalProb = False

if (CalcTotalProb):

    #l_up = np.pi
    #b_up = np.pi/2.
    b_0 = (15./90.) * np.pi/2.
    l_0 = np.pi/2.
    l_lo = 0.
    l_up = l_0
    b_lo = -b_0
    b_up = b_0
    P_1q = Psn_int(l_lo,l_up,b_lo,b_up)
    l_lo = l_0
    l_up = np.pi
    b_lo = -b_0
    b_up = b_0
    P_2q = Psn_int(l_lo,l_up,b_lo,b_up)
    l_lo = -np.pi
    l_up = -l_0
    b_lo = -b_0
    b_up = b_0
    P_3q = Psn_int(l_lo,l_up,b_lo,b_up)
    l_lo = -l_0
    l_up = 0
    b_lo = -b_0
    b_up = b_0
    P_4q = Psn_int(l_lo,l_up,b_lo,b_up)
    
    P_tot = P_1q + P_2q + P_3q + P_4q
    print ("P_1q, P_2q, P_3q, P_41, P_tot = %.4f, %.4f, %.4f, %.4f, %.4f" \
           % (P_1q, P_2q, P_3q, P_4q,P_tot))



dP_dOmega = np.zeros((N_l,N_b))

lat = np.zeros((N_l,N_b))

long = np.zeros((N_l,N_b))


P_sum = 0.
P_int_sum = 0.
P_ext_sum = 0.

#b_lim = np.zeros(N_l1q)
#l_lim = np.zeros(N_l1q)
b_lim = np.zeros(N_l)
l_lim = np.zeros(N_l)


#for i in range(0,N_l1q):
for i in range(0,N_l):

    l_deg = l_max_deg * (1. - 2.*np.float(i)/np.float(N_l-1))

    l_rad = l_deg*np.pi/180.

    b_lim[i] = 10. * (1 - l_deg / 90.)
    l_lim[i] = l_deg

    if (i%10 == 0):
        print ("%i " % i, end='', flush=True)
    if (i == N_l-1):
        print ("")

    for j in range(0,N_b):

        b_deg = b_max_deg * (2.*np.float(j)/np.float(N_b-1)-1.)

        b_rad = b_deg*np.pi/180.

        r_lim = 6.*R_sun

        dPdOmega, err = integrate.quad(dPsn_drdOmega,0.,r_lim)


        P_sum = P_sum + dPdOmega * np.cos(b_rad)
        if (b_deg < 10.* (1 - l_deg/ 90.)):
            P_int_sum += dPdOmega * np.cos(b_rad)
        else:
            P_ext_sum += dPdOmega * np.cos(b_rad)


        dP_dOmega[i,j] = dPdOmega

        long[i,j] = l_deg

        lat[i,j] = b_deg



dl_deg = lat[0,1] - lat[0,0]
db_deb = long[0,0] - long[1,0]
dOmega = dl_deg * db_deb * (np.pi/180.)**2
#print ("dOmega, P_sum, P_int_sum, P_ext_sum: %.3e %.3e %.3e %.3e" % (dOmega, P_sum, P_int_sum, P_ext_sum))
P_sum = P_sum * dOmega
P_int_sum = P_int_sum * dOmega
P_ext_sum = P_ext_sum * dOmega


print( "SN type: %s; SN distrubtion: %s" % (SN_type,SN_dist))
DataFile.write( "SN type: %s\n" % SN_label)


P_max = np.max(dP_dOmega)
maxlocs = np.where(dP_dOmega==P_max)
maxlats = lat[(dP_dOmega==P_max)]
maxlongs = long[(dP_dOmega==P_max)]
maxPs = dP_dOmega[(dP_dOmega==P_max)]
print ("there is %i maximum:" % np.size(maxlats))
print ("l [deg]\tb [deg] \tdP/dOmega_max [sr^-1]")
for j in range(np.size(maxlats)):
    print ("%6.2f \t%6.2f \t%.2e" % (maxlats[j],maxlongs[j],maxPs[j]))
    
DataFile.write("there are %i maxima:\n" % np.size(maxlats))
DataFile.write("l [deg]\tb [deg]\tdP/dOmega_max [sr^-1]\n")
for j in range(np.size(maxlats)):
    DataFile.write("%6.2f \t%7.2f\t%.2e\n" % (maxlats[j],maxlongs[j],maxPs[j]))


P_min = np.min(dP_dOmega)
minlocs = np.where(dP_dOmega==P_min)
minlats = lat[(dP_dOmega==P_min)]
minlongs = long[(dP_dOmega==P_min)]
minPs = dP_dOmega[(dP_dOmega==P_min)]
print ("there are %i minima:" % np.size(minlats))
print ("l [deg]\tb [deg]\tdP/dOmega_min [sr^-1]")
for j in range(np.size(minlats)):
    print ("%6.2f \t%7.2f\t%.2e" % (minlats[j],minlongs[j],minPs[j]))

DataFile.write("there are %i minima:\n" % np.size(minlats))
DataFile.write("l [deg]\tb [deg]\tdP/dOmega_min\n")
for j in range(np.size(minlats)):
    DataFile.write("%7.2f\t%6.2f \t%.2e\n" % (minlats[j],minlongs[j],minPs[j]))

P_max_deg2 = P_max * (np.pi/180.)**2
print ("max probability density:  %.2e deg^-2" % P_max_deg2)
DataFile.write("max probability density:  %.2e deg^-2\n" % P_max_deg2)


#P_tot = 4.*P_sum*(b_max_deg/(N_b1q-1.))*(l_max_deg/(N_l1q-1.))*(np.pi/180.)**2
#print ("estimated P_tot = ",P_tot)

f_int = P_int_sum/P_sum
f_ext = P_ext_sum/P_sum 
print ("interior probability:  P_int = %.4f, f_int = %.4f" % (P_int_sum,f_int))
print ("exterior probability:  P_ext = %.4f, f_ext = %.4f" % (P_ext_sum,f_ext))
print ("total probability:  P_tot = %.4f" % (P_sum))

DataFile.write("interior probabilty:  P_int = %.4f, f_int = %.4f\n" % (P_int_sum,f_int))
DataFile.write("exterior probabilty:  P_ext = %.4f, f_ext = %.4f\n" % (P_ext_sum,f_ext))
DataFile.write("total probabilty:  P_tot = %.4f\n" % (P_sum))


t1 = time.time()
print ("time to calculate:  %.2f sec" % (t1-t0))

####


fig = plt.figure(figsize=(15.,8))

if (PlotType == PlotType_CC26Al):
    plt.title(r"%s Supernova Probability Density: ${}^{\mathbf{26}}{\mathbf{Al}}$ Model" % (SN_label),weight='bold',fontsize=20)
else:
    plt.title("%s Supernova Probability Density" % (SN_label),weight='bold',fontsize=20)

ScaleType = 2 # [P/P_max]
ScaleType = 3 # [P/P_max] Aitoff
ScaleType = 0 # [sr^-1] linear
ScaleType = 1 # [sr^-1] log


if (ScaleType != 3):
    ax1 = plt.subplot()
    ax1.set_xlim(180.,-180.)
    ax1.set_ylim(-15.,15.)

# cs = ax1.contour(long,lat,dP_dOmega/P_max,levs)
#cs = ax1.contourf(long,lat,dP_dOmega/P_max,levs, cmap=plt.cm.jet)
#cs = ax1.contourf(long,lat,dP_dOmega,levs, cmap=plt.cm.jet, locator=ticker.LogLocator())


Punits = "deg2"


if (ScaleType == 0):
    ############################################################
    ## scale units [sr^-1]
    ############################################################
    if (Punits == "sr"):
        levs =  P_max * np.linspace(0.00, 1.00, 301)
        #cs = ax1.contourf(long,lat,dP_dOmega, levs, cmap=plt.cm.nipy_spectral)
        #cs = ax1.contourf(long,lat,dP_dOmega, levs, cmap=plt.cm.PuRd)
        cs = ax1.contourf(long,lat,dP_dOmega, levs, cmap=plt.cm.jet)
        #cs = ax1.contourf(long,lat,dP_dOmega, levs, cmap=plt.cm.jet, \
            #                  vmin=0.00, vmax=50.0)
        #plt.clim(0.,2.)
        ScaleLabel = r'Probability $dP/d\Omega \ [\rm sr^{-1}]$'
        scale_flabel = "sr"
    elif (Punits == "deg2"):
        levs =  P_max * np.linspace(0.00, 1.00, 301) * (np.pi/180.)**2
        cs = ax1.contourf(long, lat, dP_dOmega*(np.pi/180.)**2, levs, cmap=plt.cm.jet)
        ScaleLabel = r'Probability $dP/d\Omega \ [\rm deg^{-2}]$'
        scale_flabel = "deg2"
elif (ScaleType == 1):
    ############################################################
    ## log scale, units [sr^-2]
    ############################################################
    if (Punits == "sr"):
        levs = P_max*np.logspace(-5.2, 0., 301)
        cs = ax1.contourf(long,lat,dP_dOmega,levs, cmap=plt.cm.jet, norm=LogNorm(), extend='both')
        ScaleLabel = r'Probability $dP/d\Omega \ [\rm sr^{-1}]$'
        scale_flabel = "sr"
    elif (Punits == "deg2"):
    ############################################################
    ## log scale, units [deg^-2]
    ############################################################
        levs = (np.pi/180.)**2 * P_max*np.logspace(-4., 0., 301)
        #cs = ax1.contourf(long,lat,(np.pi/180.)**2*dP_dOmega,levs, cmap=plt.cm.jet, norm=LogNorm())
        cs = ax1.contourf(long,lat,(np.pi/180.)**2*dP_dOmega,levs, cmap=plt.cm.seismic, norm=LogNorm())
        ScaleLabel = r'Probability $dP/d\Omega \ [\rm deg^{-2}]$'
        scale_flabel = "deg2"
elif (ScaleType == 2):    
    ############################################################
    ## P/P_max:  dimensionless
    ############################################################
    levs = [0.001,0.003,0.01,0.03,0.1,0.3]
    levs = [0.0625,0.125,0.25,0.5,0.95]
    levs = [0.01,0.03,0.1,0.3,0.99]
    #levs = np.linspace(0.00, 0.99, 301)
    levs = np.linspace(0.00, 1.0, 301)
    cs = ax1.contourf(long,lat,dP_dOmega/P_max,levs, cmap=plt.cm.jet)
    ScaleLabel = r'Probability $P/P_{\rm max}$'
    scale_flabel = "prel"
elif (ScaleType == 3):
    ############################################################
    ## full Aitoff projection
    ############################################################
    ax1 = plt.subplot(111,projection="aitoff")
    ax1.grid(True)
    levs = [0.1,0.5]
    #ax1.plot(long/(2.*np.pi),lat/(2.*np.pi),dP_dOmega/P_max,'r.')
    #cs = ax1.contour(long,lat,dP_dOmega/P_max,levs)
    l_radian = long * np.pi/180.
    b_radian = lat * np.pi/180.
    cs = ax1.contour(l_radian,b_radian,dP_dOmega/P_max,levs)
    #cs = ax1.imshow
    ScaleLabel = r'Probability $P/P_{\rm max}$'
    scale_flabel = "Ait"
else:
    print ("Bad ScaleType")
    
print ("contour plotted")


ax1.set_xlabel(r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold")

ax1.set_ylabel(r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold")

#ax1.text(-0.9*l_max_deg,+2.*b_max_deg/3.,labtext,fontsize=20,weight='bold', color = 'white')


if (SN_type=="Ia"):
    plt.scatter([327.6-360., 4.5, 120.1], [+14.6, +6.8, +1.4], facecolors='w', edgecolors='y', zorder=10, marker=(5,1), s=200)
    plt.annotate('SN1006', (327.6-360.,+14.6), xytext=(-37.,14.2), color='w', zorder=10, fontsize=15)
    plt.annotate('SN1604 (Kepler)',(4.5, 6.8), xytext=(0.0,6.4), color = 'w',zorder=10, fontsize=15)
    plt.annotate('SN1572 (Tycho)', (120.1, 1.4), xytext=(116.,1.0), color = 'w',zorder=10, fontsize=15)
elif (SN_type=="CC"):
    plt.scatter([-175.4, 130.7], [-5.8, +3.1], facecolors='blue', marker=(5,1), zorder=10, s=200, edgecolors='y')
    plt.annotate('SN1054 (Crab)', (-175.4, -5.8), xytext=(-110,-6.2), color='k',zorder=10, fontsize=15)
    plt.annotate('SN1181',(130.7,3.1), xytext=(125,2.8), color='k', zorder=10, fontsize=15)
    maybe_CC_color = 'azure'
    maybe_CC_color = 'lightskyblue'
    maybe_CC_color = 'deepskyblue'
    plt.scatter([315.4-360., 11.2, 111.7], [-2.3, -0.3, -2.1], facecolors=maybe_CC_color, marker=(5,1), zorder=10, s=200, edgecolors='royalblue')
    plt.annotate('SN185', (-45.,-4.), xytext=(-45.,-4.), color=maybe_CC_color, zorder=10, fontsize=15)
    plt.annotate('SN386', (25.,-1.75), xytext=(25.,-1.75), color=maybe_CC_color, zorder=10, fontsize=15)
    plt.annotate('Cas A', (125.,-4.), xytext=(125.,-4.), color=maybe_CC_color, zorder=10, fontsize=15)

    



draw_LSST = True
draw_LSST = False

if (draw_LSST):
    labtext = "region excluded from baseline WFD"
    ax1.text(0.,-b_max_deg/3.,labtext,fontsize=20,weight='bold', color = 'red', horizontalalignment='center')
    plt.plot(l_lim, b_lim, 'r-', linewidth = 3)
    plt.plot(l_lim,-1 * b_lim, 'r-', linewidth = 3)
    plt.plot(-1*l_lim, b_lim, 'r-', linewidth = 3)
    plt.plot(-1*l_lim, -1*b_lim, 'r-', linewidth = 3)


############################################################
## add celestial equator

def AddCelEq():

    dec_NGP_deg = 27.12835323 
    delta_NGP = dec_NGP_deg * np.pi/180.
    
    RA_NGP_deg = 192.85949646
    alpha_NGP = RA_NGP_deg * np.pi/180.
    
    alpha_cel = np.linspace(0.,2.*np.pi,200)
    d_alpha = alpha_cel - alpha_NGP
    
    b_ceq_rad = np.arcsin(np.cos(delta_NGP)*np.cos(d_alpha))
    b_ceq = b_ceq_rad * 180./np.pi
    dl_ceq_rad = np.arctan2(np.sin(d_alpha)/np.cos(b_ceq_rad),-np.sin(delta_NGP)*np.cos(d_alpha)/np.cos(b_ceq_rad))
    dl_ceq_rad = np.arctan2(np.sin(d_alpha),-np.sin(delta_NGP)*np.cos(d_alpha)) 
    dl_ceq = dl_ceq_rad * 180./np.pi

    l_NCP = 122.932
    
    l_ceq = l_NCP - dl_ceq

    #      1234567 1234567
    print (" l_ceq   b_ceq ")
    for jj in range(200):
        if (l_ceq[jj] > 180.):
            l_ceq[jj] = l_ceq[jj] - 360.
            #print ("%3i %7.2f %7.2f" % (jj, l_ceq[jj], b_ceq[jj]))
            
            #print ("l_ceq = ",l_ceq)
            #print ("b_ceq = ",b_ceq)

            sort_inds = np.argsort(l_ceq)

            l_ceq_rad = np.zeros(200)
            b_ceq_rad = np.zeros(200)
            
            for jj in range(200):
                l_ceq_rad[jj] = l_ceq[sort_inds[jj]] * np.pi/180.
                b_ceq_rad[jj] = b_ceq[sort_inds[jj]] * np.pi/180.
                
                
                #l_ceq_rad = l_ceq * np.pi/180.
                #b_ceq_rad = b_ceq * np.pi/180.
                plt.plot(l_ceq_rad, b_ceq_rad, 'y-', linewidth = 3)


PlotCelEq = False                
                
if (PlotCelEq):
    AddCelEq()

                
############################################################

# ax1.set_aspect('equal')

if (ScaleType == 1):
    formatter = LogFormatter(10, labelOnlyBase=False)
    #cbar = fig.colorbar(cs, pad=0.01, ticks=[1.e-5, 1.e-4, 1.e-3, 1.e-2], format=formatter)
    cbar = fig.colorbar(cs, pad=0.01, ticks=[1.e-7,1.e-6,1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.,10.,100.])
    cbar.ax.tick_params(labelsize=16)
    #cbar = fig.colorbar(cs, pad=0.01)
    cbar.set_label(ScaleLabel, fontsize = 20, weight = 'bold')
#else:
elif (ScaleType != 3):
    cbar = fig.colorbar(cs, format='%.2e', pad=0.01)
    cbar.set_label(ScaleLabel, fontsize = 20, weight = 'bold')
    
#cbar.set_label(r'Probability $P/P_{\rm max}$', fontsize = 20, weight = 'bold')
#cbar.set_label(ScaleLabel, fontsize = 20, weight = 'bold')

#    cbar.set_ticks(cbticks)
#    cbar.set_ticklabels(cbticks)


plt.tick_params(axis='both', which='major', labelsize=20)



basename = "SNProb_"
zsunlabel = "zsun%.0f" % (z_sun*1.e3)
thindisklabel = "Rthin%.1f_hthin%.0f" % (R_thin,h_thin*1.e3)
thickdisklabel = "Rthick%.1f_hthick%.0f" % (R_thick,h_thick*1.e3)
reslabel = "res%ix%i" % (N_l,N_b)
figbasename = basename + SN_type + "_" + SN_dist + "_" + zsunlabel + "_" + thindisklabel + "_" + thickdisklabel + "_" + reslabel + "_" + scale_flabel

if ZoomLong:
    plotname = figbasename + "_zoom"
else:
    plotname = figbasename 

#plt.show()

figname_png = figbasename+".png"
fig.savefig(figname_png)
print ("plot written to: %s" % figname_png)
fig.savefig(figbasename+".eps")
fig.savefig(figbasename+".pdf")

t2 = time.time()

print ("total time:  %.2f sec" % (t2-t0))


PrintPoints = False
PrintPoints = True

if (PrintPoints):
    DataFile_points.write ("%i %i # N_l N_b\n" % (N_l, N_b))
    for i in range(N_l):
        for j in range(N_b):
            DataFile_points.write ("%.3f %.3f %.3e\n" % (lat[i,j], long[i,j], dP_dOmega[i,j]))

    print ("Data ponts written to:  %s" % (dataname_points))


DataFile.close()
DataFile_points.close()



def LongDis():
    print ("\nLong Dis")

    dataname_long = figbasename + "_longdist.dat"
    DataFileLong = open(dataname_long,"w")

    dP_dl_Simp = np.zeros(N_l)

    for jj in range(N_l):
        b_rads = lat[jj,...] * np.pi/180.
        dP_dOmega_slice = dP_dOmega[jj,...]
        dP_dl_integrand = np.cos(b_rads) * dP_dOmega_slice
        dP_dl_Simp[jj] = integrate.simps(dP_dl_integrand, b_rads)
        
        #print ("l,dP/dl = %.3f %.3e" \
        #       % (long[jj,0],dP_dl_Simp[jj]))
        DataFileLong.write ("%.3f %.3e\n" \
               % (long[jj,0],dP_dl_Simp[jj]))
        
    DataFileLong.close()
               
    return    

LongDis()


def LatDis():
    print ("\nLat Dis")

    dataname_lat = figbasename + "_latdist.dat"
    DataFileLat= open(dataname_lat,"w")

    dP_db_Simp = np.zeros(N_b)

    for jj in range(N_b):
        l_rads = long[...,jj] * np.pi/180.
        dP_dOmega_slice = dP_dOmega[...,jj]
        dP_db_integrand = dP_dOmega_slice
        ##### minus sign since long distribution goes from high to low
        dP_db_Simp[jj] = - integrate.simps(dP_db_integrand, l_rads) 
        
        #print ("b,dP/db = %.3f %.3e" \
        #       % (lat[0,jj],dP_db_Simp[jj]))
        DataFileLat.write ("%.3f %.3e\n" \
               % (lat[0,jj],dP_db_Simp[jj]))
        
    DataFileLat.close()
               
    return    

LatDis()

def LatitudeDist():

    print ("\nLatitude Distribution")
    
    dP_db = np.zeros((N_l,4))
    
    dataname_lat = figbasename + "_latdist.dat"
    DataFileLat= open(dataname_lat,"w")
    
    for jj in range(N_b):
        b_rads = b_radian[0,jj]
        b_degs = lat[0,jj]
        dPdb_integrand = lambda r, l_rads:  dP_drdldb(r,b_rads,l_rads)
        l_inf = l_radian[-1,0]
        l_sup = l_radian[0,0]
        r_inf = 0.
        
        l_mid = l_sup/2.
        
        dPdb1, err = integrate.dblquad(dPdb_integrand,0,l_mid,r_inf,lambda ll: FindDist(m_lim,M_SN,ll,b_rads))
        dPdb2, err = integrate.dblquad(dPdb_integrand,l_mid,l_sup,r_inf,lambda ll: FindDist(m_lim,M_SN,ll,b_rads))
        dPdb3, err = integrate.dblquad(dPdb_integrand,l_inf,-l_mid,r_inf,lambda ll: FindDist(m_lim,M_SN,ll,b_rads))
        dPdb4, err = integrate.dblquad(dPdb_integrand,-l_mid,0,r_inf,lambda ll: FindDist(m_lim,M_SN,ll,b_rads))
        dP_db[jj,0] = dPdb1
        dP_db[jj,1] = dPdb2
        dP_db[jj,2] = dPdb3
        dP_db[jj,3] = dPdb4
        dP_db_tot = dPdb1 + dPdb2 + dPdb3 + dPdb4
        
        #print ("b,dP/db = %.2f %.3e %.3e %.3e %.3e %.3e" \
        #       % (b_degs,dP_db_tot,dP_db[jj,0],dP_db[jj,1],dP_db[jj,2],dP_db[jj,3]))
        DataFileLat.write ("%.3f %.3e %.3e %.3e %.3e %.3e\n" % (b_degs,dP_db_tot,dP_db[jj,0],dP_db[jj,1],dP_db[jj,2],dP_db[jj,3]))
        
        DataFileLat.close()
        return
    


def LongitudeDist():

    print ("\nLongitude Distribution")
    
    dP_dl = np.zeros(N_l)
    
    dataname_long = figbasename + "_longdist.dat"
    DataFileLong = open(dataname_long,"w")
    
    for jj in range(N_l):
        l_rads = l_radian[jj,0]
        l_degs = long[jj,0]
        dPdl_integrand = lambda r, b_rads:  dP_drdldb(r,b_rads,l_rads)
        b_inf = b_radian[0,0]
        b_sup = b_radian[0,-1]
        r_inf = 0.

        dPdl, err = integrate.dblquad(dPdl_integrand,b_inf,b_sup,r_inf,lambda bb: FindDist(m_lim,M_SN,l_rads,bb))
        dP_dl[jj] = dPdl

        #print ("l,dP/dl = %.2f %.3e" % (l_degs,dP_dl[jj]))
        DataFileLong.write ("%.3f %.3e\n" % (l_degs, dP_dl[jj]))

        DataFileLong.close()
        return

    
