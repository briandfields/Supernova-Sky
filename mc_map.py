import numpy as np

import matplotlib.pyplot as plt

#from mpl_toolkits.mplot3d
#import Axes3D

from numpy.random import rand

from scipy import integrate

import math

import time

from scipy.stats import kde

from sys import exit


def SNstats(ell,bee):
    Nquad = np.zeros(4)
    for j in range(len(ell)):
        if ( (ell[j] > 0.) and (bee[j] > 0.) ):
            Nquad[0] = Nquad[0] + 1
        elif ( (ell[j] < 0.) and (bee[j] > 0.) ):
            Nquad[1] = Nquad[1] + 1
        elif ( (ell[j] > 0.) and (bee[j] < 0.) ):
            Nquad[2] = Nquad[2] + 1
        elif ( (ell[j] < 0.) and (bee[j] < 0.) ):
            Nquad[3] = Nquad[3] + 1
        else:
            print "should not be here"

    for j in range(len(Nquad)):
        print "quadrant %i: %i events" % (j, Nquad[j])
        
    return

t0 = time.time()

# Variables, numbers and strings

R_sun = 8.5# kpc, distance from galactic center to sun

h_sun = 0.02 # kpc
h_sun = 0.00 # kpc
#kpc, height of sun above galactic plane

print ("Sun's Galactocentric coords: (R,z) = (%.2f, %.3f) kpc" % (R_sun, h_sun))

R_thin = 2.9 # kpc
# kpc, scale radius of Milky Way

h_thin = 0.095 # kpc
# kpc, scale height of Milky Way

R_thick = 2.4  # kpc
#kpc

h_thick = 0.8 #kpc

stddev = 1
#magnitude difference of standard deviation

sn_num = int(4e5)
sn_num = int(9e3)
sn_num = int(1e4)
sn_num = int(1e5)
sn_num = int(1e5)
sn_num = int(2.e3)
sn_num = int(1e6)
sn_num = int(2e4)
sn_num = int(5e3)
sn_num = int(1e5)
sn_num = int(5.e3)
sn_num = int(1.e3)
sn_num = int(1e6)


# number of random supernovae generated

band = "K"
band = "V"
#band we're looking in


SN_type = "Ia"
SN_label = "Type Ia"


SN_type = "Ia"
SN_label = "Type Ia"

SN_type = "CC"
SN_label = "Core Collapse"


#sn type we're looking at

minvis = 0.0
minvis = 2.0
# 2: easily noticeable everywhere; 6: noticeable in dark spaces; 7: human vision limit; 26: LSST


print ("Running with:")
print ("N_SN = %i" % sn_num)
print ("Type:  %s = %s" % (SN_type, SN_label))
print ("band:  %s" % band)
print ("m_lim = %.2f mag" % minvis)


if ((SN_type != "CC") and (SN_type != "Ia")):
    print "bad SN type ::%2::" % (SN_type)



if band == "V":

    A_gc = 30

    if SN_type == "CC": absmag = -16

    elif SN_type == "Ia": absmag = -18.5



elif band == "K":

        A_gc = 3.51

        if SN_type == "CC": absmag = -17

        elif SN_type == "Ia": absmag = -17.8



else:

    #exit("Invalid band, please use either V or K band")
    print "Invalid band, please use either V or K band"







extinct = absmag - minvis #total exinction along sightline needed









# STEP 1: generate random supernovae

np.random.seed(137137)

mu1 = rand(sn_num) # random values to use for r, theta, and z

mu2 = rand(sn_num)

mu3 = rand(sn_num)





def radius_dist(y,u):

    return np.log((y+1)/(1-u))



a = np.zeros(sn_num)

k = [0.5]*sn_num

o = np.zeros(sn_num)

for i in range(0,sn_num):

    while a[i] != 1:

        o[i] = radius_dist(k[i],mu1[i])

        if abs(o[i]-k[i]) <= 0.0001:

            a[i] = 1

        else:

            k[i] = o[i]





r = np.zeros_like(o)

theta = np.zeros_like(o)

z = np.zeros_like(o)

if SN_type == "CC":

    # r = R_thin*o

    # theta = 2*np.pi*mu2 # uniform angular distribution from 0 to 2pi
    
    # z = h_thin*np.log(1/1-mu3) # exponential height distribution
    
    for i in range(len(o)):

        r[i] = R_thin*o[i]

        theta[i] = 2.*np.pi*mu2[i]
        # uniform angular distribution from 0 to 2pi

        z[i] = h_thin*np.log(1./(1.-mu3[i]))
        # exponential height distribution

        if i%2 == 0:
            z[i] = - z[i]


elif SN_type == "Ia":

    for i in range(len(o)):

        if i%2 == 0:

            r[i] = R_thin*o[i]

            theta[i] = 2.*np.pi*mu2[i]
            # uniform angular distribution from 0 to 2pi

            z[i] = h_thin*np.log(1./(1.-mu3[i]))
            if i%4 == 0:
                z[i] = - z[i]
            # exponential height distribution

        else:

            r[i] = R_thick*o[i]

            theta[i] = 2*np.pi*mu2[i]
            # uniform angular distribution from 0 to 2pi
            
            z[i] = h_thick*np.log(1./(1.-mu3[i]))
            if i%4 == 3:
                z[i] = - z[i]
            # exponential height distribution
            
else:

    print ("SN type = ::%s::") % (SN_type)
    exit("Invalid SN type, make sure it is either CC or Ia")



print(sn_num, "supernovae generated...")





x = np.zeros_like(r)

y = np.zeros_like(r)



for i in range(len(mu1)):

    x[i] = r[i]*np.cos(theta[i]) # used for plotting in step 5
    
    y[i] = r[i]*np.sin(theta[i]) # same use as x array
    



# STEP 2: find relative position to Sun

### set sun at x = r_sun, y = z = 0



l = np.zeros_like(x) # galactic longitude w.r.t. Sun

b = np.zeros_like(y) # galactic latitude w.r.t. Sun



radius = np.sqrt((x-R_sun)**2 + y**2 + (z-h_sun)**2)

b = np.arcsin((z-h_sun)/radius)

for i in range(len(radius)):

    if x[i] <= R_sun:
        #arcsin only covers from -pi/2 to pi/2

        l[i] = np.arcsin(y[i] / np.sqrt((x[i]-R_sun)**2 + y[i]**2))


    else:
        # cover the rest of the sky (whole thing goes from -pi/2 to 3pi/2)
        
        #l[i] = np.arcsin(y[i] / np.sqrt((x[i]-R_sun)**2 + y[i]**2))  + np.pi
        l[i] = np.arcsin(y[i] / np.sqrt((x[i]-R_sun)**2 + y[i]**2))
        if (y[i] > 0.):
            l[i] = np.pi - l[i] 
        else:
            l[i] = -l[i] - np.pi

    #b[i] = np.arcsin((z[i]-h_sun) / radius[i])



#getting arrays in degrees and centering around 0 for plotting

l_deg = np.zeros_like(l)

b_deg = np.zeros_like(b)

for i in range(len(l)):

    l_deg[i] = math.degrees(l[i])

    b_deg[i] = math.degrees(b[i])


for j in range(len(l_deg)):

    if l_deg[j] > 180:

        l_deg[j] = l_deg[j]-360


print "for all supernovae: in (l,b)"
SNstats(l_deg,b_deg)
print "for all supernovae: in (b,b)"
SNstats(b_deg,b_deg)
print "for all supernovae: in (l,l)"
SNstats(l_deg,l_deg)
print "for all supernovae:  in (x,y)"
SNstats(x,y)
print "for all supernovae:  in (y,z)"
SNstats(y,z)
print "for all supernovae:  in (x,z)"
SNstats(x,z)



# STEP 3: get dimming for each supernova

def rho_dust(R,z):

    #"Dust Density Function"

    R_thin = 2.9  # kpc

    h_thin = 0.095 # kpc


    rho = np.exp(-R/R_thin)*np.exp(-np.abs(z)/h_thin)

    rho = rho / (R_thin * (1. - np.exp(-R_sun/R_thin)))

    return rho





def dAv_dr(radius, l_rad,b_rad):

#"Extinction Rate due to Dust"

    z_cyl = h_sun + radius * np.sin(b_rad) # solarcentric radius component in plane

    r_par = radius * np.cos(b_rad) # Galactocentric cylindrical radius
    
    R_cyl = np.sqrt(r_par**2 - 2.*R_sun*r_par*np.cos(l_rad) + R_sun**2)



    Av_gc = 30.0

    dAv_dr = Av_gc * rho_dust(R_cyl,z_cyl)


    return dAv_dr



def distmod(rad):

#"Distance Modulus, takes radii in kpc"

    return 5*np.log10(rad*1000/10)



def dustdim(l_rad,b_rad,radi):

#"Total dimming due to dust"

    Sigfunct = lambda r: dAv_dr(r, l_rad, b_rad)
    #so dAv_dr will intgrate correctly
    
    Sigma, err = integrate.quad(Sigfunct, 0. ,radi)
    #get magnitude loss due to extinction

    return Sigma



dimmed = np.zeros_like(radius)

for i in range(len(radius)):

    dimmed[i] = distmod(radius[i]) + dustdim(l[i], b[i], radius[i])



print('Extinction calculated...')




# STEP 4: sort visible and invisible supernova

mask = dimmed > -extinct

darkx = x[mask] 

darky = y[mask] # invisible supernovae

darkz = z[mask]

#title = "Estimated Probability after " + str(sn_num) + " Simulated Supernovae"

newmask = ~mask

brightx = x[newmask]

brighty = y[newmask] # visible supernovae

brightz = z[newmask]



brightl = l_deg[newmask]

brightb = b_deg[newmask]

diml = l_deg[mask]
dimb = b_deg[mask]

print('Supernovae sorted...')

print('%1.2f percent of supernovae have an apparent magnitude greater than %d' % (len(brightx)/sn_num*100, minvis))

# totl = np.append(brightl, brightl)

# totb = np.append(brightb, -brightb)


print "Stats for visible supernovae"
Nq_tot = SNstats(brightl,brightb)


t1 = time.time()
hours = (t1-t0)/3600
minutes = (hours - int(hours))*60
seconds = (minutes-int(minutes))*60

print('Time to run: %d hours %d minutes %.1f seconds' %(int(hours), int(minutes), seconds))


print('Plotting supernovae...')

fig0 = plt.figure(figsize=(15.,8.))

#plt.scatter(diml, dimb, c='k', marker='.')
#plt.scatter(brightl, brightb, c='tab:red', marker='.')
plt.plot(diml, dimb, 'k.',zorder=1)
plt.plot(brightl, brightb, 'r.', zorder=2)

plt.title('%d %s Supernovae in the %s Band' %(sn_num, SN_label, band),fontsize=20)

plt.text(-170.,13.,r"$%s_{\rm max} > %.1f$" % (band,minvis),color='blue',fontsize=20)

plt.xlabel(r'Galactic Longitude $\ell$',fontsize=20)

plt.ylabel(r'Galactic Latitude $b$',fontsize=20)

plt.xlim(-180.,180.)
plt.ylim(-10.,10.)

plt.savefig("scat.png")
                  
                  

fig = plt.figure(figsize=(15.,8.))

ax1 = plt.subplot()

nbins=1000

X = brightl
Y = brightb

k = kde.gaussian_kde([X,Y])
xi, yi = np.mgrid[X.min():X.max():nbins*1j, Y.min():Y.max():nbins*1j]

zi = k(np.vstack([xi.flatten(), yi.flatten()]))
ax1.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.jet)
#ax1.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)

plt.title('%d %s Supernovae in the %s Band' %(sn_num, SN_label, band),fontsize=20)

ax1.set_xlabel(r'Galactic Longitude $\ell$ [deg]',fontsize=20)

ax1.set_ylabel(r'Galactic Latitude $b$ [deg]',fontsize=20)

#mag_label = r"$%s_{\rm max} > %.1f$" % (band,minvis),
plt.text(+170.,13.,r"$%s_{\rm max} > %.1f$" % (band,minvis),color='white',fontsize=20)
#plt.text(-160.,8.,mag_label,color='white',fontsize=20)

ax1.set_xlim(180.,-180.)
ax1.set_ylim(-15.,15.)



if (SN_type=="Ia"):
    plt.scatter([327.6-360., 4.5, 120.1], [+14.6, +6.8, +1.4], facecolors='w', edgecolors='y', zorder=10, marker=(5,1), s=200)
    plt.annotate('SN1006', (327.6-360.,+14.6), xytext=(-37.,14.2), color='w', zorder=10, fontsize=15)
    plt.annotate('SN1604 (Kepler)',(4.5, 6.8), xytext=(2.5,6.4), color = 'w',zorder=10, fontsize=15)
    plt.annotate('SN1572 (Tycho)', (120.1, 1.4), xytext=(116.,1.0), color = 'w',zorder=10, fontsize=15)
elif (SN_type=="CC"):
    plt.scatter([-175.4, 130.7], [-5.8, +3.1], facecolors='w', marker=(5,1), zorder=10, s=200, edgecolors='y')
    plt.annotate('SN1054 (Crab)', (-175.4, -5.8), xytext=(-125,-6.0), color='w',zorder=10, fontsize=15)
    plt.annotate('SN1181',(130.7,3.1), xytext=(127,3.0), color='w', zorder=10, fontsize=15)

magname = "%s%.1f" % (band,minvis)
lgsnnum = "%.1f" % (np.log10(sn_num))
zsunlabel = "zsun%.0f" % (h_sun*1.e3)
figbasename = "SkyMap_" + SN_type + "_" + magname + "_" + lgsnnum + "_" + zsunlabel
figname_png = figbasename+".png"
figname_pdf = figbasename+".pdf"
plt.savefig(figname_png)
print ("plot written to: %s" % figname_png)
#plt.savegig(figname_pdf)
#plt.savefig("vis.pdf")



# plt.subplot(111, projection = 'aitoff')

# plt.hist2d(totl, totb, bins = (360, 180), normed = True, cmap = plt.cm.jet)

#plt.scatter(brightl, brightb, c = 'r', marker = '.')

# k = kde.gaussian_kde([X,Y])
# xi, yi = np.mgrid[X.min():X.max():nbins*1j, Y.min():Y.max():nbins*1j]
# zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
# # Make the plot
# plt.pcolormesh(xi, yi, zi.reshape(xi.shape))


# plt.title('%d %s Supernovae in the %s-Band' %(sn_num, SN_type, band))

# plt.xlabel('Longitude')

# plt.ylabel('Latitude')

# plt.grid()

# plt.ylim(-10,10)

# #plt.colorbar(label = 'Probability per Square Degree')

# plt.show()



t2 = time.time()
hours = (t2-t0)/3600
minutes = (hours - int(hours))*60
seconds = (minutes-int(minutes))*60

print('Time to run: %d hours %d minutes %.1f seconds' %(int(hours), int(minutes), seconds))
        
