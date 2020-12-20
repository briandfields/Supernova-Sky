##########################################################
##
## SNSummary



import numpy as np

import matplotlib.pyplot as plt

import matplotlib.colorbar as colorbar

import scipy.integrate as integrate

import matplotlib.ticker as ticker

import astropy as ap
from astropy.coordinates import SkyCoord

from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter

from math import gamma

from scipy.stats import kde

import time



Data_visCC = "./DATA/SNVisProb_CC_CC_V2.0_M-15.0_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"
Data_visIa = "./DATA/SNVisProb_Ia_Ia_V2.0_M-14.7_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"

ccfile = open(Data_visCC)
line1 = ccfile.readline()
line1 = line1.strip()
print ("line1 = ",line1)
N_l_CC, N_b_CC, juk1, junk2, junk3 = line1.split(" ")
ccfile.close()

N_l_CC = int(N_l_CC)
N_b_CC = int(N_b_CC)

iafile = open(Data_visIa)
line1 = iafile.readline()
line1 = line1.strip()
print ("line1 = ",line1)
N_l_Ia, N_b_Ia, j1, j2, j3 = line1.split(" ")
iafile.close()

N_l_Ia = int(N_l_Ia)
N_b_Ia = int(N_b_Ia)

lat_CC, long_CC, dP_dOmega_CC = np.loadtxt(Data_visCC, dtype='float', delimiter=' ', skiprows=1, usecols=(0,1,2), unpack=True)

lat_Ia, long_Ia, dP_dOmega_Ia = np.loadtxt(Data_visIa, dtype='float', delimiter=' ', skiprows=1, usecols=(0,1,2), unpack=True)

lat_CC = np.reshape(lat_CC, (N_l_CC, N_b_CC))
long_CC = np.reshape(long_CC, (N_l_CC, N_b_CC))
dP_dOmega_CC = np.reshape(dP_dOmega_CC, (N_l_CC, N_b_CC))

lat_Ia = np.reshape(lat_Ia, (N_l_Ia, N_b_Ia))
long_Ia = np.reshape(long_Ia, (N_l_Ia, N_b_Ia))
dP_dOmega_Ia = np.reshape(dP_dOmega_Ia, (N_l_Ia, N_b_Ia))



############################################################
## full Aitoff projection
############################################################

fig = plt.figure(figsize=(15.,8))


#ax1 = plt.subplot(111,projection="aitoff")
ax1 = plt.subplot(111,projection="mollweide")
#ax1.set_ylim(-np.pi/4.,np.pi/4.)

ax1.grid(True)

levs = [0.1,0.5,0.95]
levs = [0.5,1.,2.]
levs = [0.5,4.,32.,100.]
levs = [0.5,2.,8.,32.]
levs = [0.1,0.3,1.,3.]
levs = [0.2,0.6,2.,6.]
levs = [0.3,1.,3.,10.]
levs = [0.25,1.,3.,10.]

colors_cc = ['cornflowerblue','royalblue','mediumblue','darkblue']
colors_ia = ['lightcoral','indianred','firebrick','maroon']

colors_cc = ['blue','mediumblue','darkblue','black']
colors_ia = ['red','firebrick','maroon','black']

colors_cc = ['cornflowerblue', 'blue','mediumblue','darkblue']
colors_ia = ['lightcoral', 'red','firebrick','maroon']

#levs = [np.exp(-2.),np.exp(-0.5)]
#colors_cc = ['royalblue','mediumblue']
#colors_ia = ['indianred','firebrick']




########## Type Ia ##########

alpha_transp = 0.6

l_radian = long_Ia * np.pi/180.
b_radian = lat_Ia * np.pi/180.

P_max_Ia = np.max(dP_dOmega_Ia)
print ("P_max_Ia = %.2e" % (P_max_Ia))

cs = ax1.contourf(-l_radian,b_radian,dP_dOmega_Ia,levs,colors=colors_ia,alpha=alpha_transp)


cs = ax1.imshow

########## core collapse ##########

l_radian = long_CC * np.pi/180.
b_radian = lat_CC * np.pi/180.

P_max_CC = np.max(dP_dOmega_CC)
print ("P_max_CC = %.2e" % (P_max_CC))

cs = ax1.contourf(-l_radian,b_radian,dP_dOmega_CC,levs,colors=colors_cc,alpha=alpha_transp)

cs = ax1.imshow


ScaleLabel = r'Probability $P/P_{\rm max}$'
    
print ("contour plotted")


#ax1.set_xlabel(r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold")

#ax1.set_ylabel(r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold")

hist_Ia_l = np.array([327.6-360., 4.5, 120.1]) * np.pi/180.
hist_Ia_b = np.array([+14.6, +6.8, +1.4]) * np.pi/180.

hist_Ia_names = np.array(["SN1006", "SN1604 (Kepler)", "SN1572 (Tycho)"])
hist_Ia_l_lab = np.array([-20., 42., 122.]) * np.pi/180.
hist_Ia_b_lab = np.array([+30., 20., 7.]) * np.pi/180.
hist_Ia_tx = np.array([-20., 42., 122.]) * np.pi/180.
hist_Ia_ty = np.array([20., 14., 7.]) * np.pi/180.

##### plot historical SNIa
Ia_color = 'r'
plt.scatter(-hist_Ia_l, hist_Ia_b, facecolors=Ia_color, edgecolors='y', zorder=10, marker=(5,1), s=200)
for j in range(3):
    plt.annotate(hist_Ia_names[j], (-hist_Ia_l_lab[j],hist_Ia_b_lab[j]), \
                xytext=(-hist_Ia_tx[j],hist_Ia_ty[j]), color='r', zorder=10, fontsize=15)

##### plot historical CC    
hist_CC_l = np.array([-175.4, 130.7]) * np.pi/180.
hist_CC_b = np.array([-5.8, +3.1]) * np.pi/180.

hist_CC_names = np.array(["SN1054 (Crab)", "SN1181"])
hist_CC_l_lab = np.array([-124., 153.]) * np.pi/180.
hist_CC_b_lab = np.array([-9., 7.]) * np.pi/180.
hist_CC_tx = np.array([-124.,153.]) * np.pi/180. 
hist_CC_ty = np.array([-9., 7.]) * np.pi/180.


plt.scatter(-hist_CC_l, hist_CC_b, facecolors='b', marker=(5,1), zorder=10, s=200, edgecolors='y')
for j in range(2):
    plt.annotate(hist_CC_names[j], (-hist_CC_l_lab[j],hist_CC_b_lab[j]), \
                xytext=(-hist_CC_tx[j],hist_CC_ty[j]), color='darkblue', zorder=10, fontsize=15)

##### plot possible CC

#draw_circle = plt.Circle((-(348.-360.)*np.pi/180.,0.*np.pi/180.), radius=8.*np.pi/180., fill=False, color='deepskyblue', zorder=10)
#ax1.gcf().gca().add_artist(draw_circle)
#ax1.add_patch(draw_circle)


maybe_CC_l = np.array([315.4-360.,11.2,111.7]) * np.pi/180.
maybe_CC_b = np.array([-2.3,-0.3,-2.1]) * np.pi/180.
maybe_CC_names = np.array(["SN185","SN386","Cas A"])
maybe_CC_l_lab = np.array([-40.,15.,125.]) * np.pi/180.
maybe_CC_b_lab = np.array([-3.,-8,3.]) * np.pi/180. 
maybe_CC_tx = np.array([-35.,20.,125.]) * np.pi/180.
maybe_CC_ty = np.array([-13.,-12.,-10.]) * np.pi/180.

maybe_CC_color = 'lightskyblue'
maybe_CC_color = 'deepskyblue'
plt.scatter(-maybe_CC_l, maybe_CC_b, facecolors=maybe_CC_color, marker=(5,1), zorder=10, s=200, edgecolors='royalblue')
for j in range(3):
    plt.annotate(maybe_CC_names[j], (-maybe_CC_l_lab[j],maybe_CC_b_lab[j]), \
                xytext=(-maybe_CC_tx[j],maybe_CC_ty[j]), color=maybe_CC_color, zorder=10, fontsize=15)


############################################################
## add celestial equator

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
ax1.plot(-l_ceq_rad, b_ceq_rad, 'y-', linewidth = 3)

ax1.set_xticklabels([150,120,90,60,30,0,-30,-60,-90,-120,-150])
ax1.set_yticklabels([-75,-60,-45,-30,-15,0.,15,30,45,60])


############################################################


plt.tick_params(axis='both', which='major', labelsize=15)


plt.title("Naked-Eye Supernovae", weight='bold', fontsize=20)

band = "V"
m_lim = 2.0
plt.text(170.*np.pi/180.,45.*np.pi/180.,r"$%s_{\rm max} > %.1f$" % (band,m_lim), color='white', fontsize=20)

plt.text(20.*np.pi/180., 67.*np.pi/180., r"celestial equator", color='y', fontsize=15)

l_NCP = 122.932
b_NCP = 27.1284

l_SCP = 302.932-360.
b_SCP = -27.1284

plt.text(-l_NCP*np.pi/180., b_NCP*np.pi/180., r"NCP", color='y', fontsize=15, multialignment='center')
plt.text(-l_SCP*np.pi/180., b_SCP*np.pi/180., r"SCP", color='y', fontsize=15, multialignment='center')


figbasename = "SNvis_sky"

figname_png = figbasename+".png"
fig.savefig(figname_png)
print ("plot written to: %s" % figname_png)

figname_eps = figbasename+".eps"
fig.savefig(figname_eps)
print ("plot written to: %s" % figname_eps)


############################################################
##### naked-eye, equatorial coordinates

PlotEquatorial = False

if (PlotEquatorial):

    fig_eq = plt.figure(figsize=(15.,8))


    ax1eq = plt.subplot(111,projection="aitoff")


    ax1eq.grid(True)

    levs = [0.1,0.5,0.95]
    levs = [0.5,1.,2.]
    levs = [0.5,4.,32.,100.]
    levs = [0.5,2.,8.,32.]
    levs = [0.3,1.,3.,10.]

    colors_cc = ['cornflowerblue','royalblue','mediumblue','darkblue']
    colors_ia = ['lightcoral','indianred','firebrick','maroon']

    colors_cc = ['blue','mediumblue','darkblue','black']
    colors_ia = ['red','firebrick','maroon','black']

    colors_cc = ['cornflowerblue', 'blue','mediumblue','darkblue']
    colors_ia = ['lightcoral', 'red','firebrick','maroon']

    #levs = [np.exp(-2.),np.exp(-0.5)]
    #colors_cc = ['royalblue','mediumblue']
    #colors_ia = ['indianred','firebrick']




    ########## Type Ia ##########

    alpha_transp = 0.6

    l_radian = long_Ia * np.pi/180.
    b_radian = lat_Ia * np.pi/180.

    P_max_Ia = np.max(dP_dOmega_Ia)
    print ("P_max_Ia = %.2e" % (P_max_Ia))

    print ("shape is " ,l_radian.shape)

    levs_eq = [0.01,0.03,0.1,0.3]
    levs_eq = levs

    ra_radian = 0.*l_radian
    dec_radian = 0.*b_radian
    for j in range(N_l_Ia):
        for k in range(N_b_Ia):
            coord = SkyCoord(l_radian[j,k],b_radian[j,k], frame='galactic', unit='radian')
            ra_radian[j,k] = coord.icrs.ra.radian
            if (ra_radian[j,k] > np.pi):
                ra_radian[j,k] = ra_radian[j,k] - 2.*np.pi
                dec_radian[j,k] = coord.icrs.dec.radian
                #print ("j, k, l, b, ra, dec = %i %i %.3f %.3f %.3f %.3f" % (j,k,l_radian[j,k],b_radian[j,k],ra_radian[j,k],dec_radian[j,k]))

    cs = ax1eq.contourf(-ra_radian,dec_radian,dP_dOmega_Ia,levs_eq,colors=colors_ia,alpha=alpha_transp)


    cs = ax1eq.imshow

    ########## core collapse ##########

    l_radian = long_CC * np.pi/180.
    b_radian = lat_CC * np.pi/180.
    
    P_max_CC = np.max(dP_dOmega_CC)
    print ("P_max_CC = %.2e" % (P_max_CC))

    cs = ax1eq.tricontourf(-ra_radian,dec_radian,dP_dOmega_CC,levs_eq,colors=colors_cc,alpha=alpha_transp)

    cs = ax1eq.imshow


    ScaleLabel = r'Probability $P/P_{\rm max}$'
    
    print ("contour plotted")


    #ax1.set_xlabel(r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold")

    #ax1.set_ylabel(r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold")

    hist_Ia_l = np.array([327.6-360., 4.5, 120.1]) * np.pi/180.
    hist_Ia_b = np.array([+14.6, +6.8, +1.4]) * np.pi/180.

    hist_Ia_names = np.array(["SN1006", "SN1604 (Kepler)", "SN1572 (Tycho)"])
    hist_Ia_l_lab = np.array([-20., 42., 140.]) * np.pi/180.
    hist_Ia_b_lab = np.array([+30., 20., -7.]) * np.pi/180.
    hist_Ia_tx = np.array([-20., 42., 140.]) * np.pi/180.
    hist_Ia_ty = np.array([20., 14., -7.]) * np.pi/180.

    hist_CC_l = np.array([-175.4, 130.7]) * np.pi/180.
    hist_CC_b = np.array([-5.8, +3.1]) * np.pi/180.

    hist_CC_names = np.array(["SN1054 (Crab)", "SN1181"])
    hist_CC_l_lab = np.array([-124., 145.]) * np.pi/180.
    hist_CC_b_lab = np.array([-9., 8.]) * np.pi/180.
    hist_CC_tx = np.array([-124.,145.]) * np.pi/180. 
    hist_CC_ty = np.array([-9., 8.]) * np.pi/180.


    plt.scatter(-hist_Ia_l, hist_Ia_b, facecolors=Ia_color, edgecolors='y', zorder=10, marker=(5,1), s=200)
    for j in range(3):
        plt.annotate(hist_Ia_names[j], (-hist_Ia_l_lab[j],hist_Ia_b_lab[j]), \
                xytext=(-hist_Ia_tx[j],hist_Ia_ty[j]), color='r', zorder=10, fontsize=15)
    
    plt.scatter(-hist_CC_l, hist_CC_b, facecolors='b', marker=(5,1), zorder=10, s=200, edgecolors='y')
    for j in range(2):
        plt.annotate(hist_CC_names[j], (-hist_CC_l_lab[j],hist_CC_b_lab[j]), \
                xytext=(-hist_CC_tx[j],hist_CC_ty[j]), color='darkblue', zorder=10, fontsize=15)

    ax1eq.set_xticklabels([150,120,90,60,30,0,-30,-60,-90,-120,-150])
    ax1eq.set_yticklabels([-75,-60,-45,-30,-15,0.,15,30,45,60])

    ylimits = np.array([-45.,+45.]) * (np.pi/180.0)
    ax1eq.set_ylim(ylimits)


    ############################################################


    plt.tick_params(axis='both', which='major', labelsize=15)


    plt.title("Naked-Eye Supernovae", weight='bold', fontsize=20)

    band = "V"
    m_lim = 2.0
    plt.text(170.*np.pi/180.,45.*np.pi/180.,r"$%s_{\rm max} > %.1f$" % (band,m_lim), color='white', fontsize=20)



    figbasename = "SNvis_eqsky"

    figname_png = figbasename+".png"
    fig_eq.savefig(figname_png)
    print ("plot written to: %s" % figname_png)

    figname_eps = figbasename+".eps"
    fig_eq.savefig(figname_eps)
    print ("plot written to: %s" % figname_eps)


############################################################
##### Intrinsic Probability and SNR Locations
############################################################


Data_allCC = "./DATA/SNVisProb_CC_V65.0_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"
Data_allIa = "./DATA/SNVisProb_Ia_V65.0_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"

ccfile = open(Data_allCC)
line1 = ccfile.readline()
line1 = line1.strip()
print ("line1 = ",line1)
N_l_CC, N_b_CC, juk1, junk2, junk3 = line1.split(" ")
ccfile.close()

N_l_CC = int(N_l_CC)
N_b_CC = int(N_b_CC)

iafile = open(Data_allIa)
line1 = iafile.readline()
line1 = line1.strip()
print ("line1 = ",line1)
N_l_Ia, N_b_Ia, j1, j2, j3 = line1.split(" ")
iafile.close()

N_l_Ia = int(N_l_Ia)
N_b_Ia = int(N_b_Ia)

lat_CC, long_CC, dP_dOmega_CC = np.loadtxt(Data_allCC, dtype='float', delimiter=' ', skiprows=1, usecols=(0,1,2), unpack=True)

lat_Ia, long_Ia, dP_dOmega_Ia = np.loadtxt(Data_allIa, dtype='float', delimiter=' ', skiprows=1, usecols=(0,1,2), unpack=True)

lat_CC = np.reshape(lat_CC, (N_l_CC, N_b_CC))
long_CC = np.reshape(long_CC, (N_l_CC, N_b_CC))
dP_dOmega_CC = np.reshape(dP_dOmega_CC, (N_l_CC, N_b_CC))

lat_Ia = np.reshape(lat_Ia, (N_l_Ia, N_b_Ia))
long_Ia = np.reshape(long_Ia, (N_l_Ia, N_b_Ia))
dP_dOmega_Ia = np.reshape(dP_dOmega_Ia, (N_l_Ia, N_b_Ia))


############################################################
#####  Green 2019 Supernova Remnant Catalog

Data_Green = 'Green2019.tsv'

snrfile = open(Data_Green)

long_SNR_deg, lat_SNR_deg = np.loadtxt(Data_Green, dtype='float', delimiter=';', skiprows=65, usecols=(0,1), unpack=True)

snrfile.close()

lat_SNR = lat_SNR_deg * np.pi/180.
long_SNR = long_SNR_deg * np.pi/180.
for j in range(len(long_SNR)):
    if (long_SNR[j] > np.pi):
        long_SNR[j] = long_SNR[j] - 2.*np.pi

############################################################
## full Aitoff projection
############################################################


fig2 = plt.figure(figsize=(15.,8))

#ax2 = plt.subplot(111,projection="aitoff")
ax2 = plt.subplot(111,projection="mollweide")
#ax1.set_ylim(-np.pi/4.,np.pi/4.)

ax2.grid(True)

#levs = [0.1,0.5,0.95]
#levs = [0.5,2.,8.]

#colors_cc = ['cornflowerblue','royalblue','mediumblue']
#colors_cc = ['blue','mediumblue','darkblue']
#colors_ia = ['lightcoral','indianred','firebrick']
#colors_ia = ['red','firebrick','maroon']

########## Type Ia ##########

alpha_transp = 0.6

l_radian = long_Ia * np.pi/180.
b_radian = lat_Ia * np.pi/180.

P_max_Ia = np.max(dP_dOmega_Ia)
print ("P_max_Ia = %.2e" % (P_max_Ia))

cs = ax2.contourf(-l_radian,b_radian,dP_dOmega_Ia,levs,colors=colors_ia,alpha=alpha_transp)


cs = ax2.imshow

########## core collapse ##########

l_radian = long_CC * np.pi/180.
b_radian = lat_CC * np.pi/180.

P_max_CC = np.max(dP_dOmega_CC)

print ("P_max_CC = %.2e" % (P_max_CC))

cs = ax2.contourf(-l_radian,b_radian,dP_dOmega_CC,levs,colors=colors_cc,alpha=alpha_transp)

cs = ax2.imshow



ax2.plot(-l_ceq_rad, b_ceq_rad, 'y-', linewidth = 3)

ax2.set_xticklabels([150,120,90,60,30,0,-30,-60,-90,-120,-150])
ax2.set_yticklabels([-75,-60,-45,-30,-15,0.,15,30,45,60])


#####

P_SNR = kde.gaussian_kde([-long_SNR, lat_SNR])

GC = np.array([0.,0.])
GACp = np.array([180.,0.])*np.pi/180.
GACm = np.array([180.,0.])*np.pi/180.
print ("P_SNR(GC,GACp,GACm) = (%.2e, %.2e, %.2e)" % (P_SNR(GC),P_SNR(GACp), P_SNR(GACm)))

PlotSNRKDE = True
PlotSNRKDE = False

if (PlotSNRKDE):
    xi, yi = np.mgrid[-180:180:360*4*1j, -90:90:180*4*1j]*np.pi/180. #Resolution of image
    zi = P_SNR(np.vstack([xi.flatten(), yi.flatten()]))

    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greys)
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.binary)
    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.jet)    


############################################################

plt.scatter(-long_SNR, lat_SNR, facecolors='k', edgecolors='black', zorder=5,alpha=0.5, marker=".")



plt.scatter(-hist_Ia_l, hist_Ia_b, facecolors=Ia_color, edgecolors='y', zorder=10, marker=(5,1), s=200)
for j in range(3):
    plt.annotate(hist_Ia_names[j], (-hist_Ia_l_lab[j],hist_Ia_b_lab[j]), \
                xytext=(-hist_Ia_tx[j],hist_Ia_ty[j]), color='r', zorder=10, fontsize=15)
    
plt.scatter(-hist_CC_l, hist_CC_b, facecolors='b', marker=(5,1), zorder=10, s=200, edgecolors='y')
for j in range(2):
    plt.annotate(hist_CC_names[j], (-hist_CC_l_lab[j],hist_CC_b_lab[j]), \
                xytext=(-hist_CC_tx[j],hist_CC_ty[j]), color='darkblue', zorder=10, fontsize=15)



plt.scatter(-maybe_CC_l, maybe_CC_b, facecolors=maybe_CC_color, marker=(5,1), zorder=10, s=200, edgecolors='royalblue')
for j in range(3):
    plt.annotate(maybe_CC_names[j], (-maybe_CC_l_lab[j],maybe_CC_b_lab[j]), \
                xytext=(-maybe_CC_tx[j],maybe_CC_ty[j]), color='cornflowerblue', zorder=10, fontsize=15)

    
plt.tick_params(axis='both', which='major', labelsize=15)


plt.title("Milky Way Supernova Remnants", weight='bold', fontsize=20)

plt.text(20.*np.pi/180., 67.*np.pi/180., r"celestial equator", color='y', fontsize=15)

l_NCP = 122.932
b_NCP = 27.1284

l_SCP = 302.932-360.
b_SCP = -27.1284

plt.text(-l_NCP*np.pi/180., b_NCP*np.pi/180., r"NCP", color='y', fontsize=15, multialignment='center')
plt.text(-l_SCP*np.pi/180., b_SCP*np.pi/180., r"SCP", color='y', fontsize=15, multialignment='center')

ylimits = np.array([-45.,+45.]) * (np.pi/180.0)
#plt.ylim(ylimits)



figbasename = "SNR_sky"

figname_png = figbasename+".png"
fig2.savefig(figname_png)
print ("plot written to: %s" % figname_png)

figname_eps = figbasename+".eps"
fig2.savefig(figname_eps)
print ("plot written to: %s" % figname_eps)

