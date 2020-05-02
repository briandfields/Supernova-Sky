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

##   zoom:  show smaller field

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
SN_dist = 'spherical'

SN_dist = 'thindisk'

SN_dist = 'Adams'



# SN_dist = 'Green'



##
## supernova type



SN_type = 'Ia'
SN_dist = 'Ia'
SN_label = "Type Ia"

SN_type = 'CC'
SN_dist = 'CC'
SN_label = "Core-Collapse"



#SN_type = SN_dist


### passband
band = "V"

m_lim = 0.
m_lim = 80.
m_lim = 2.

N_obscure = 0

if (SN_type == 'CC'):
    M_SN = -16.
elif (SN_type == 'Ia'):
    M_SN = -17.8


## plot region is zoomed in longitidue

#zoom = True

zoom = False

## resolution: longitude and latitude points in one quadrant

#N_l = 2*N_l1q
#N_b = 2*N_b1q

#N_l = 180
#N_b = 10
N_l = 180+1
N_b = 15+1
N_l = 2*360+1
N_b = 2*30+1
N_l = 2*45+1
N_b = 2*30+1

#N_l1q = 2*90
#N_b1q = 2*10
#N_l1q = 90
#N_b1q = 10


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




def q_SN(R_gc,z_gc):

    # supernova density in Galactocentric cylindrical coordinates

    global R_disk, h_disk

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
        return r_max
    
    P_1q, err = integrate.tplquad(dP_drdldb,l_min,l_max,b_inf,b_sup,r_inf,r_sup)
    P_sn = 4. * P_1q

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
R_thin = 2.9 # kpc

h_thick = 0.800 # kpc
R_thick = 2.4 # kpc



if (zoom):

    l_max_deg = 90.

else:

    l_max_deg = 180.



if (SN_type == 'Ia'):

    h_disk = h_ia

    R_disk = R_ia
 # kpc

    b_max_deg = 15.

    labtext = "Type Ia Supernovae"

    if (SN_dist == 'thindisk'):

        h_disk = 0.350 # kpc
elif (SN_type == 'CC'):

    h_disk = h_cc

    R_disk = R_cc

    b_max_deg = 15.

    labtext = "Core-Collapse Supernovae"

    if (SN_dist == 'thindisk'):

        h_disk = 0.350 # kpc

else:

    print ("bad SN type!")





print ("here goes P_SN")

#l_up = np.pi
#b_up = np.pi/2.
b_0 = (15./90.) * np.pi/2.
l_0 = np.pi/2.
l_lo = 0.
l_up = l_0
b_lo = 0.
b_up = b_0
#P_indisk = Psn_int(l_lo,l_up,b_lo,b_up)
l_lo = l_0
l_up = np.pi
b_lo = 0.
b_up = b_0
#P_outdisk = Psn_int(l_lo,l_up,b_lo,b_up)
l_lo = 0.
l_up = l_0
b_lo = b_0
b_up = np.pi/2.
#P_inpole = Psn_int(l_lo,l_up,b_lo,b_up)
l_lo = l_0
l_up = np.pi
b_lo = b_0
b_up = np.pi/2.
#P_outpole = Psn_int(l_lo,l_up,b_lo,b_up)

#P_tot = P_indisk + P_outdisk + P_inpole + P_outpole
#print "P_indisk, P_outdis, P_inpole, P_outpole, P_tot = %.4f, %.4f, %.4f, %.4f, %.4f" \
#    % (P_indisk, P_outdisk, P_inpole, P_outpole,P_tot)



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
    #l_deg = l_max_deg * (2.*np.float(i)/np.float(N_l-1) - 1.)
    #l_deg = l_max_deg * np.float(i)/np.float(N_l1q-1)
    #l_deg = l_max_deg * i/(N_l1q-1)

    l_rad = l_deg*np.pi/180.

    b_lim[i] = 10. * (1 - l_deg / 90.)
    l_lim[i] = l_deg

    if (i%10 == 0):
        print ("%i " % i, end='', flush=True)
    if (i == N_l-1):
        print ("")

    #for j in range(0,N_b1q):
    for j in range(0,N_b):

        b_deg = b_max_deg * (2.*np.float(j)/np.float(N_b-1)-1.)
        #b_deg = b_max_deg * float(j)/np.float(N_b1q-1)
        #b_deg = b_max_deg * j/(N_b1q-1)

        b_rad = b_deg*np.pi/180.

        r_lim = FindDist(m_lim,M_SN,l_rad,b_rad)

        dPdOmega, err = integrate.quad(dPsn_drdOmega,0.,r_lim)


        P_sum = P_sum + dPdOmega * np.cos(b_rad)
        if (b_deg < 10.* (1 - l_deg/ 90.)):
            P_int_sum += dPdOmega * np.cos(b_rad)
        else:
            P_ext_sum += dPdOmega * np.cos(b_rad)


        dP_dOmega[i,j] = dPdOmega

        long[i,j] = l_deg

        lat[i,j] = b_deg


        #dP_dOmega[N_l1q-i-1,N_b1q+j] = dPdOmega

        #dP_dOmega[N_l1q+i,N_b1q+j] = dPdOmega

        #dP_dOmega[N_l1q-i-1,N_b1q-j-1] = dPdOmega

        #dP_dOmega[N_l1q+i,N_b1q-j-1] = dPdOmega


        #long[N_l1q-i-1,N_b1q-j-1] = l_deg

        #long[N_l1q-i-1,N_b1q+j] = l_deg

        #long[N_l1q+i,N_b1q-j-1] = -l_deg

        #long[N_l1q+i,N_b1q+j] = -l_deg


        #lat[N_l1q-i-1,N_b1q-j-1] = -b_deg

        #lat[N_l1q-i-1,N_b1q+j] = b_deg

        #lat[N_l1q+i,N_b1q-j-1] = -b_deg

        #lat[N_l1q+i,N_b1q+j] = b_deg



dl_deg = lat[0,1] - lat[0,0]
db_deb = long[0,0] - long[1,0]
dOmega = dl_deg * db_deb * (np.pi/180.)**2
print ("dOmega, P_sum, P_int_sum, P_ext_sum: %.3e %.3e %.3e %.3e" % (dOmega, P_sum, P_int_sum, P_ext_sum))
P_sum = P_sum * dOmega
P_int_sum = P_int_sum * dOmega
P_ext_sum = P_ext_sum * dOmega


print( "SN type ", SN_type)


P_max = np.max(dP_dOmega)
P_max_deg2 = P_max * (np.pi/180.)**2
print ("max probability density:  %.2e deg^-2" % P_max_deg2)


#P_tot = 4.*P_sum*(b_max_deg/(N_b1q-1.))*(l_max_deg/(N_l1q-1.))*(np.pi/180.)**2
#print ("estimated P_tot = ",P_tot)

f_int = P_int_sum/P_sum
f_ext = P_ext_sum/P_sum 
print ("interior propabilty:  P_int = %.4f, f_int = %.4f" % (P_int_sum,f_int))
print ("exterior propabilty:  P_ext = %.4f, f_ext = %.4f" % (P_ext_sum,f_ext))
print ("total propabilty:  P_tot = %.4f" % (P_sum))

print ("%i of %i obscured sightlines" % (N_obscure,N_l*N_b))


t1 = time.time()
print ("time to calculate:  %.2f sec" % (t1-t0))

####

levs = [0.001,0.003,0.01,0.03,0.1,0.3]

levs = [0.0625,0.125,0.25,0.5,0.95]

levs = [0.01,0.03,0.1,0.3,0.99]

#levs = np.linspace(0.00, 0.99, 301)
levs = np.linspace(0.00, 1.0, 301)

fig = plt.figure(figsize=(15.,8))

#plt.title("Sky Map: Type Ia Supernova Probability",weight='bold',fontsize=20)
plt.title("%s Supernovae in the %s Band" % (SN_label, band),weight='bold',fontsize=20)

plt.text(170.,13.,r"$%s_{\rm max} > %.1f$" % (band,m_lim), color='white', fontsize=20)


ax1 = plt.subplot()

ax1.set_xlim(180.,-180.)
ax1.set_ylim(-15.,15.)

# cs = ax1.contour(long,lat,dP_dOmega/P_max,levs)
cs = ax1.contourf(long,lat,dP_dOmega/P_max,levs, cmap=plt.cm.jet)

print ("contour plotted")

#ax1 = plt.subplot(111,projection="aitoff")

#ax1.grid(True)
#ax1.plot(long/(2.*np.pi),lat/(2.*np.pi),dP_dOmega/P_max,'r.')

#cs = ax1.contour(lat/(2.*np.pi),long/(2.*np.pi),dP_dOmega/P_max,levs)

#cs = ax1.contour(long,lat,dP_dOmega/P_max,levs)

#cs = ax1.imshow

ax1.set_xlabel(r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold")

ax1.set_ylabel(r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold")

#ax1.text(-0.9*l_max_deg,+2.*b_max_deg/3.,labtext,fontsize=20,weight='bold', color = 'white')


if (SN_type=="Ia"):
    plt.scatter([327.6-360., 4.5, 120.1], [+14.6, +6.8, +1.4], facecolors='w', edgecolors='y', zorder=10, marker=(5,1), s=200)
    plt.annotate('SN1006', (327.6-360.,+14.6), xytext=(-37.,14.2), color='w', zorder=10, fontsize=15)
    plt.annotate('SN1604 (Kepler)',(4.5, 6.8), xytext=(0.0,6.4), color = 'w',zorder=10, fontsize=15)
    plt.annotate('SN1572 (Tycho)', (120.1, 1.4), xytext=(116.,1.0), color = 'w',zorder=10, fontsize=15)
elif (SN_type=="CC"):
    plt.scatter([-175.4, 130.7], [-5.8, +3.1], facecolors='w', marker=(5,1), zorder=10, s=200, edgecolors='y')
    plt.annotate('SN1054 (Crab)', (-175.4, -5.8), xytext=(-110,-6.2), color='w',zorder=10, fontsize=15)
    plt.annotate('SN1181',(130.7,3.1), xytext=(125,2.8), color='w', zorder=10, fontsize=15)




draw_LSST = True
draw_LSST = False

if (draw_LSST):
    labtext = "region excluded from baseline WFD"
    ax1.text(0.,-b_max_deg/3.,labtext,fontsize=20,weight='bold', color = 'red', horizontalalignment='center')
    plt.plot(l_lim, b_lim, 'r-', linewidth = 3)
    plt.plot(l_lim,-1 * b_lim, 'r-', linewidth = 3)
    plt.plot(-1*l_lim, b_lim, 'r-', linewidth = 3)
    plt.plot(-1*l_lim, -1*b_lim, 'r-', linewidth = 3)



# ax1.set_aspect('equal')

cbar = fig.colorbar(cs, format='%.2f', pad=0.01)
cbar.set_label(r'Probability $P/P_{\rm max}$', fontsize = 20, weight = 'bold')

plt.tick_params(axis='both', which='major', labelsize=20)



basename = "SNVisProb_"
magname = "%s%.1f" % (band,m_lim)
zsunlabel = "zsun%.0f" % (z_sun*1.e3)
thindisklabel = "Rthin%.1f_hthin%.0f" % (R_thin,h_thin*1.e3)
thickdisklabel = "Rthick%.1f_hthick%.0f" % (R_thick,h_thick*1.e3)
figbasename = basename + SN_type + "_" + magname + "_" + zsunlabel + "_" + thindisklabel + "_" + thickdisklabel

if zoom:
    plotname = figbasename + "_zoom"
else:
    plotname = figbasename 

#plt.show()

fig.savefig(figbasename+".png")
#fig.savefig(figbasename+".eps")
#fig.savefig(figbasename+".pdf")

t2 = time.time()

print ("total time:  %.2f sec" % (t2-t0))
