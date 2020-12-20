##########################################################
##
## SNSummary



import numpy as np

import matplotlib.pyplot as plt

import matplotlib.colorbar as colorbar
import matplotlib.colors
import matplotlib.patches as mpatches

from matplotlib.lines import Line2D
import scipy.integrate as integrate

import matplotlib.ticker as ticker

import astropy as ap
from astropy.coordinates import SkyCoord

from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter

from math import gamma

from scipy.stats import kde

import time




################################################################################
##### Galactic to Mollweide coordinate conversions


def GalCoords_to_Moll(ell_rad,bee_rad):

    if (np.abs(np.abs(bee_rad)-np.pi/2.) < 1.e-3):
        theta = bee_rad
        N_rep = -1

    else:

        theta_0 = bee_rad

        eps = 1.e-4
        
        theta_n = 0.
        theta_n1 = theta_0
        
        N_rep = 0
        while (np.abs(theta_n - theta_n1) > eps*np.abs(theta_n1)):
            N_rep = N_rep + 1
            theta_n = theta_n1
            theta_n1 = theta_n - (2.*theta_n+np.sin(2.*theta_n)-np.pi*np.sin(bee_rad))/(2.+2.*np.cos(2.*theta_n))

        theta = theta_n1

    #print ("phi, theta, N = %.3f %.3f %i" % (bee_rad, theta, N_rep))

    
    x_Moll = 2.*np.sqrt(2.)*ell_rad * np.cos(theta) / np.pi
    y_Moll = np.sqrt(2.) * np.sin(theta)

    MollCoords = ([x_Moll,y_Moll])

    return MollCoords


def Moll_to_GalCoords(x_M,y_M):

    theta_M = np.arcsin(y_M/np.sqrt(2.))

    ell_rad = np.pi*x_M/(2.*np.sqrt(2.)*np.cos(theta_M))

    bee_rad = np.arcsin((2.*theta_M+np.sin(2.*theta_M))/np.pi)

    GalCoords = ([ell_rad,bee_rad])

    return GalCoords



################################################################################


Data_visCC = "./DATA/SNVisProb_CC_CC_V2.0_M-15.0_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"
Data_visIa = "./DATA/SNVisProb_Ia_Ia_V2.0_M-14.7_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"


Data_visCC = "./DATA/SNVisProb_CC_CC_V2.0_M-15.0_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res181x361_prob.dat"
Data_visIa = "./DATA/SNVisProb_Ia_Ia_V2.0_M-14.7_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res181x361_prob.dat"

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
## Mollweide projection
############################################################

fig = plt.figure(figsize=(15.,8))


#ax1 = plt.subplot(111,projection="aitoff")
ax1 = plt.subplot(111)
x_max, y_dum = GalCoords_to_Moll(np.pi, 0.)
x_min, y_dum = GalCoords_to_Moll(-np.pi, 0.)
ax1.set_xlim(x_max*1.01, x_min*1.01)
b_max = 30. * np.pi/180.
b_min = - b_max
x_dum, y_max = GalCoords_to_Moll(0., b_max)
x_dum, y_min = GalCoords_to_Moll(0., b_min)
ax1.set_ylim(y_min*1.01, y_max*1.01)

#ax1.grid(True)

levs = [0.1,0.5,0.95]
levs = [0.5,1.,2.]
levs = [0.5,4.,32.,100.]
levs = [0.5,2.,8.,32.]
levs = [0.1,0.3,1.,3.]
levs = [0.2,0.6,2.,6.]
levs = [0.3,1.,3.,10.]
levs = [0.25,1.,3.,10.]
levs = ([0.25, 1.])

colors_cc = ['cornflowerblue','royalblue','mediumblue','darkblue']
colors_ia = ['lightcoral','indianred','firebrick','maroon']

colors_cc = ['blue','mediumblue','darkblue','black']
colors_ia = ['red','firebrick','maroon','black']

colors_cc = ['cornflowerblue', 'blue','mediumblue','darkblue']
colors_ia = ['lightcoral', 'red','firebrick','maroon']
colors_ia = ['red','firebrick','maroon']

#colors_cc = ['cornflowerblue', 'blue']
#colors_ia = ['lightcoral', 'red']


#levs = [np.exp(-2.),np.exp(-0.5)]
#colors_cc = ['royalblue','mediumblue']
#colors_ia = ['indianred','firebrick']




########## Type Ia ##########

# find the 68% CL

P_Ia_68 = 0.13
P_Ia_95 = 0.0070
P_Ia_0 = P_Ia_68

N_Ia_vis_in = 0
N_Ia_vis_tot = dP_dOmega_Ia.size
I_in = 0.
I_tot = 0.
for i in range(N_l_Ia):
    for j in range(N_b_Ia):
        I_tot = I_tot + dP_dOmega_Ia[i,j]
        if (dP_dOmega_Ia[i,j] > P_Ia_0):
            I_in = I_in + dP_dOmega_Ia[i,j]
            N_Ia_vis_in = N_Ia_vis_in + 1
            

#f_Ia_vis_in = N_Ia_vis_in/(N_Ia_vis_tot+0.)
f_Ia_vis_in = I_in/I_tot

print ("visible Ia: P=%.3f contains f_Ia_vis = %.4f" % \
       (P_Ia_0, f_Ia_vis_in))

levs_Ia = ([P_Ia_95, P_Ia_68, 100.])

############################################################

alpha_transp = 0.6

l_radian = long_Ia * np.pi/180.
b_radian = lat_Ia * np.pi/180.

P_max_Ia = np.max(dP_dOmega_Ia)
print ("P_max_Ia = %.2e" % (P_max_Ia))

x_Ia = 0.*l_radian
y_Ia = 0.*b_radian

print ("l_radian size:", l_radian.shape[0], l_radian.shape[1])
print ("loop begins")
for i in range(l_radian.shape[0]):
    for j in range(l_radian.shape[1]):
        #print (l_radian[i,j], b_radian[i,j])
        x_Ia[i,j], y_Ia[i,j] = GalCoords_to_Moll(l_radian[i,j], b_radian[i,j])

        
cs_Ia = ax1.contourf(x_Ia,y_Ia,dP_dOmega_Ia,levels=levs_Ia,colors=colors_ia,alpha=alpha_transp,)

cs_Ia = ax1.imshow

########## core collapse ##########

# find the 68% CL

P_CC_68 = 0.17
P_CC_95 = 0.0143
P_CC_0 = P_CC_68

N_CC_vis_in = 0
N_CC_vis_tot = dP_dOmega_CC.size
I_in = 0.
I_tot = 0.
for i in range(N_l_Ia):
    for j in range(N_b_Ia):
        I_tot = I_tot + dP_dOmega_CC[i,j]
        if (dP_dOmega_CC[i,j] > P_CC_0):
            I_in = I_in + dP_dOmega_CC[i,j]
            N_CC_vis_in = N_CC_vis_in + 1
            

#f_Ia_vis_in = N_Ia_vis_in/(N_Ia_vis_tot+0.)
f_CC_vis_in = I_in/I_tot

print ("visible CC: P=%.3f contains f_CC_vis = %.4f" % \
       (P_CC_0, f_CC_vis_in))

levs_CC = ([P_CC_95, P_CC_68, 100.])

#####

l_radian = long_CC * np.pi/180.
b_radian = lat_CC * np.pi/180.

x_CC = 0.*l_radian
y_CC = 0.*b_radian

for i in range(l_radian.shape[0]):
    for j in range(l_radian.shape[1]):
        x_CC[i,j], y_CC[i,j] = GalCoords_to_Moll(l_radian[i,j], b_radian[i,j])


P_max_CC = np.max(dP_dOmega_CC)
print ("P_max_CC = %.2e" % (P_max_CC))

cs_CC = ax1.contourf(x_CC,y_CC,dP_dOmega_CC,levels=levs_CC,colors=colors_cc,alpha=alpha_transp)

cs_CC = ax1.imshow

#facecols = [colors_cc[1],colors_cc[0],colors_ia[1],colors_ia[0]]
#cmap = matplotlib.colors.ListedColormap(facecols)
#proxy = [ax1.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) for pc in cs_CC.collections]
#proxy = [plt.Rectangle((0,0),1,1,fc = cmap)]
#plt.legend(proxy, ['CCSN 68%','CCSN 95%','SNIa 68%','SNIa 95%'])

#hc,lc = cs_CC.legend_elements("dP_dOmega_CC")
#hi,li = cs_Ia.legend_elements("dP_dOmega_Ia")
#ax1.legend(hc+hi, lc+li)

cc68patch = mpatches.Patch(color=colors_cc[1], label='CC 68%')
cc95patch = mpatches.Patch(color=colors_cc[0], label='CC 95%')
ia68patch = mpatches.Patch(color=colors_ia[1], label='Ia 68%')
ia95patch = mpatches.Patch(color=colors_ia[0], label='Ia 95%')
ceqline = Line2D([0], [0], color='y', linewidth=3, linestyle='-',label=r'Celestial Equator')
plt.rc('legend', fontsize=16)
legloc = ([0.05,0.03])
ax1.legend(handles=[cc68patch,cc95patch,ia68patch,ia95patch,ceqline], loc=legloc, framealpha=1.)



ScaleLabel = r'Probability $P/P_{\rm max}$'
    
print ("contour plotted")


ax1.set_xlabel(r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold")

ax1.set_ylabel(r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold")

plt.text(0., 1.25*y_min, r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold",horizontalalignment="center")
plt.text(1.15*x_max, 0., r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold",verticalalignment="center",rotation=90)


hist_Ia_l = np.array([327.6-360., 4.5, 120.1]) * np.pi/180.
hist_Ia_b = np.array([+14.6, +6.8, +1.4]) * np.pi/180.


hist_Ia_names = np.array(["SN1006", "SN1604 (Kepler)", "SN1572 (Tycho)"])
hist_Ia_l_lab = np.array([-20., 42., 122.]) * np.pi/180.
hist_Ia_b_lab = np.array([+30., 20., 7.]) * np.pi/180.
hist_Ia_tx = np.array([-17., 24., 122.]) * np.pi/180.
hist_Ia_ty = np.array([17., 10., 5.]) * np.pi/180.

x_hIa = 0.*hist_Ia_l
y_hIa = 0.*hist_Ia_b
x_hIa_lab = 0.*hist_Ia_l
y_hIa_lab = 0.*hist_Ia_b
x_hIa_tx = 0.*hist_Ia_l
y_hIa_ty = 0.*hist_Ia_b

for j in range(hist_Ia_l.size):
    x_hIa[j], y_hIa[j] = GalCoords_to_Moll(hist_Ia_l[j], hist_Ia_b[j])
    x_hIa_lab[j], y_hIa_lab[j] = GalCoords_to_Moll(hist_Ia_l_lab[j], hist_Ia_b_lab[j])
    x_hIa_tx[j], y_hIa_ty[j] = GalCoords_to_Moll(hist_Ia_tx[j], hist_Ia_ty[j])



##### plot historical SNIa
Ia_color = 'r'
plt.scatter(x_hIa, y_hIa, facecolors=Ia_color, edgecolors='y', zorder=10, marker=(5,1), s=200)
for j in range(3):
    plt.annotate(hist_Ia_names[j], (x_hIa_lab[j],y_hIa_lab[j]), \
                xytext=(x_hIa_tx[j],y_hIa_ty[j]), color='r', zorder=10, fontsize=15)

##### plot historical CC    
hist_CC_l = np.array([-175.4, 130.7]) * np.pi/180.
hist_CC_b = np.array([-5.8, +3.1]) * np.pi/180.


hist_CC_names = np.array(["SN1054 (Crab)", "SN1181"])
hist_CC_l_lab = np.array([-124., 153.]) * np.pi/180.
hist_CC_b_lab = np.array([-9., 7.]) * np.pi/180.
hist_CC_tx = np.array([-124.,159.]) * np.pi/180. 
hist_CC_ty = np.array([-6., 7.]) * np.pi/180.

x_hCC = 0.*hist_CC_l
y_hCC = 0.*hist_CC_b
x_hCC_lab = 0.*hist_CC_l
y_hCC_lab = 0.*hist_CC_b
x_hCC_tx = 0.*hist_CC_l
y_hCC_ty = 0.*hist_CC_b

for j in range(hist_CC_l.size):
    x_hCC[j], y_hCC[j] = GalCoords_to_Moll(hist_CC_l[j], hist_CC_b[j])
    x_hCC_lab[j], y_hCC_lab[j] = GalCoords_to_Moll(hist_CC_l_lab[j], hist_CC_b_lab[j])
    x_hCC_tx[j], y_hCC_ty[j] = GalCoords_to_Moll(hist_CC_tx[j], hist_CC_ty[j])

plt.scatter(x_hCC, y_hCC, facecolors='b', marker=(5,1), zorder=10, s=200, edgecolors='y')
for j in range(2):
    plt.annotate(hist_CC_names[j], (x_hCC_lab[j],y_hCC_lab[j]), \
                xytext=(x_hCC_tx[j],y_hCC_ty[j]), color='darkblue', zorder=10, fontsize=15)

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
maybe_CC_ty = np.array([-8.,-4.,-7.]) * np.pi/180.

x_mCC = 0.*maybe_CC_l
y_mCC = 0.*maybe_CC_b
x_mCC_lab = 0.*maybe_CC_l
y_mCC_lab = 0.*maybe_CC_b
x_mCC_tx = 0.*maybe_CC_l
y_mCC_ty = 0.*maybe_CC_b

for j in range(maybe_CC_l.size):
    x_mCC[j], y_mCC[j] = GalCoords_to_Moll(maybe_CC_l[j], maybe_CC_b[j])
    x_mCC_lab[j], y_mCC_lab[j] = GalCoords_to_Moll(maybe_CC_l_lab[j], maybe_CC_b_lab[j])
    x_mCC_tx[j], y_mCC_ty[j] = GalCoords_to_Moll(maybe_CC_tx[j], maybe_CC_ty[j])


maybe_CC_color = 'lightskyblue'
maybe_CC_color = 'deepskyblue'

plt.scatter(x_mCC, y_mCC, facecolors=maybe_CC_color, marker=(5,1), zorder=10, s=200, edgecolors='y')
for j in range(3):
    plt.annotate(maybe_CC_names[j], (x_mCC_lab[j],y_mCC_lab[j]), \
                xytext=(x_mCC_tx[j],y_mCC_ty[j]), color=maybe_CC_color, zorder=10, fontsize=15)


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

x_ceq = np.zeros(200)
y_ceq = np.zeros(200)

for jj in range(200):
    l_ceq_rad[jj] = l_ceq[sort_inds[jj]] * np.pi/180.
    b_ceq_rad[jj] = b_ceq[sort_inds[jj]] * np.pi/180.
    x_ceq[jj], y_ceq[jj] = GalCoords_to_Moll(l_ceq_rad[jj], b_ceq_rad[jj])
    
    
#l_ceq_rad = l_ceq * np.pi/180.
#b_ceq_rad = b_ceq * np.pi/180.
ax1.plot(x_ceq, y_ceq, 'y-', linewidth = 3)

#l_grid = ([150.,120.,90.,60.,30.,0.,-30.,-60.,-90.,-120.,-150.]) * np.pi/180.
l_grid_deg = np.array([180.,150.,120.,90.,60.,30.,0.,-30.,-60.,-90.,-120.,-150.,-180.])
l_grid = l_grid_deg * np.pi/180.
#b_grid = ([-45.,-30.,-15.,0.,15.,30.,45.]) * np.pi/180.
b_grid_deg = np.array([-30.,-15.,0.,15.,30.]) 
b_grid = b_grid_deg * np.pi/180.

b_range = np.linspace(y_min, y_max, 100)
l_range = np.linspace(-np.pi, np.pi, 100)

x_vals = 0.*b_range
y_vals = 0.*b_range

for i in range(l_grid.size):
    for j in range(b_range.size):
        x_vals[j], y_vals[j] = GalCoords_to_Moll(l_grid[i], b_range[j])
    lwid = 0.5
    if ((i == 0) or (i == l_grid.size-1)):
        lwid = 3
    ax1.plot(x_vals,y_vals,'k-',linewidth=lwid)
    ax1.text(x_vals[i], 1.04*y_vals[0], r"%i" % l_grid_deg[i], horizontalalignment="center",verticalalignment="bottom",fontsize=20)
    
for i in range(b_grid.size):
    for j in range(l_range.size):
        x_vals[j], y_vals[j] = GalCoords_to_Moll(l_range[j], b_grid[i])
    lwid = 0.5
    if ((i == 0) or (i == b_grid.size-1)):
        lwid = 3
    ax1.plot(x_vals,y_vals,'k-',linewidth=lwid)
    ax1.text(1.02*x_vals[-1], y_vals[i], r"%i" % b_grid_deg[i], horizontalalignment="right",verticalalignment="center",fontsize=20)
    

ax1.set_xticklabels([150,120,90,60,30,0,-30,-60,-90,-120,-150])
ax1.set_yticklabels([-75,-60,-45,-30,-15,0.,15,30,45,60])

ax1.axis('off')


############################################################


plt.tick_params(axis='both', which='major', labelsize=15)


plt.title("Naked-Eye Supernovae", weight='bold', fontsize=20)

band = "V"
m_lim = 2.0

x_m,y_m = GalCoords_to_Moll(170.*np.pi/180.,25.*np.pi/180.)
#plt.text(x_m, y_m ,r"$%s_{\rm max} > %.1f$" % (band,m_lim), color='yellow', fontsize=20)

x_m,y_m = GalCoords_to_Moll(15.*np.pi/180., 28.*np.pi/180.)
#plt.text(x_m, y_m, r"celestial equator", color='y', fontsize=15)

l_NCP = 122.932
b_NCP = 27.1284

x_NCP, y_NCP = GalCoords_to_Moll(l_NCP*np.pi/180., b_NCP*np.pi/180.)

l_SCP = 302.932-360.
b_SCP = -27.1284

x_SCP, y_SCP = GalCoords_to_Moll(l_SCP*np.pi/180., b_SCP*np.pi/180.)


plt.text(x_NCP, y_NCP, r"NCP", color='y', fontsize=15, multialignment='center')
plt.text(x_SCP, y_SCP, r"SCP", color='y', fontsize=15, multialignment='center')


figbasename = "SNvis_sky_Moll"

figname_png = figbasename+".png"
fig.savefig(figname_png)
print ("plot written to: %s" % figname_png)

figname_eps = figbasename+".eps"
fig.savefig(figname_eps)
print ("plot written to: %s" % figname_eps)



############################################################
##### Intrinsic Probability and SNR Locations
############################################################


Data_allCC = "./DATA/SNVisProb_CC_V65.0_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"
Data_allIa = "./DATA/SNVisProb_Ia_V65.0_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res91x61_prob.dat"

Data_allCC = "./DATA/SNProb_CC_CC_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res181x361_prob.dat"
Data_allIa = "./DATA/SNProb_Ia_Ia_zsun20_Rthin2.9_hthin95_Rthick2.4_hthick800_res181x361_prob.dat"

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
## full Mollweide projection
############################################################


fig2 = plt.figure(figsize=(15.,8))

#ax2 = plt.subplot(111,projection="aitoff")
#ax2 = plt.subplot(111,projection="mollweide")
ax2 = plt.subplot(111)
ax2.set_xlim(x_max*1.01, x_min*1.01)
ax2.set_ylim(y_min*1.01, y_max*1.01)

ax2.axis('off')

for i in range(l_grid.size):
    for j in range(b_range.size):
        x_vals[j], y_vals[j] = GalCoords_to_Moll(l_grid[i], b_range[j])
    lwid = 0.5
    if ((i == 0) or (i == l_grid.size-1)):
        lwid = 3
    ax2.plot(x_vals,y_vals,'k-',linewidth=lwid)
    ax2.text(x_vals[i], 1.04*y_vals[0], r"%i" % l_grid_deg[i], horizontalalignment="center",verticalalignment="bottom",fontsize=20)

    
for i in range(b_grid.size):
    for j in range(l_range.size):
        x_vals[j], y_vals[j] = GalCoords_to_Moll(l_range[j], b_grid[i])
    lwid = 0.5
    if ((i == 0) or (i == b_grid.size-1)):
        lwid = 3
    ax2.plot(x_vals,y_vals,'k-',linewidth=lwid)
    ax2.text(1.02*x_vals[-1], y_vals[i], r"%i" % b_grid_deg[i], horizontalalignment="right",verticalalignment="center",fontsize=20)


#ax2.grid(True)

#levs = [0.1,0.5,0.95]
#levs = [0.5,2.,8.]
levs = ([1.,100.])

#colors_cc = ['cornflowerblue','royalblue','mediumblue']
#colors_cc = ['blue','mediumblue','darkblue']
#colors_ia = ['lightcoral','indianred','firebrick']
#colors_ia = ['red','firebrick','maroon']

########## Type Ia ##########
# find the 68% CL

P_Ia_68 = 0.94
P_Ia_95 = 0.050
P_Ia_0 = P_Ia_95

N_Ia_all_in = 0
N_Ia_all_tot = dP_dOmega_Ia.size
I_in = 0.
I_tot = 0.
for i in range(N_l_Ia):
    for j in range(N_b_Ia):
        I_tot = I_tot + dP_dOmega_Ia[i,j]
        if (dP_dOmega_Ia[i,j] > P_Ia_0):
            I_in = I_in + dP_dOmega_Ia[i,j]
            N_Ia_all_in = N_Ia_all_in + 1
            

f_Ia_all_in = I_in/I_tot

print ("intrinsic Ia: P=%.3f contains f_Ia_all = %.4f" % \
       (P_Ia_0, f_Ia_all_in))

levs_Ia = ([P_Ia_95, P_Ia_68, 100.])

############################################################

alpha_transp = 0.6

l_radian = long_Ia * np.pi/180.
b_radian = lat_Ia * np.pi/180.

x_Ia = 0.*l_radian
y_Ia = 0.*b_radian

print ("l_radian size:", l_radian.shape[0], l_radian.shape[1])
print ("loop begins")
for i in range(l_radian.shape[0]):
    for j in range(l_radian.shape[1]):
        #print (l_radian[i,j], b_radian[i,j])
        x_Ia[i,j], y_Ia[i,j] = GalCoords_to_Moll(l_radian[i,j], b_radian[i,j])


P_max_Ia = np.max(dP_dOmega_Ia)
print ("P_max_Ia = %.2e" % (P_max_Ia))

cs = ax2.contourf(x_Ia,y_Ia,dP_dOmega_Ia,levels=levs_Ia,colors=colors_ia,alpha=alpha_transp)


cs = ax2.imshow

########## core collapse ##########

# find the 68% CL

P_CC_68 = 4.65
P_CC_95 = 0.205
P_CC_0 = P_CC_68

N_CC_all_in = 0
N_CC_all_tot = dP_dOmega_CC.size
I_in = 0.
I_tot = 0.
for i in range(N_l_Ia):
    for j in range(N_b_Ia):
        I_tot = I_tot + dP_dOmega_CC[i,j]
        if (dP_dOmega_CC[i,j] > P_CC_0):
            I_in = I_in + dP_dOmega_CC[i,j]
            N_CC_all_in = N_CC_all_in + 1
            

#f_Ia_vis_in = N_Ia_vis_in/(N_Ia_vis_tot+0.)
f_CC_all_in = I_in/I_tot

print ("intrinsic CC: P=%.3f contains f_CC_all = %.4f" % \
       (P_CC_0, f_CC_all_in))

levs_CC = ([P_CC_95, P_CC_68, 100.])

##################################################

l_radian = long_CC * np.pi/180.
b_radian = lat_CC * np.pi/180.

x_CC = 0.*l_radian
y_CC = 0.*b_radian

for i in range(l_radian.shape[0]):
    for j in range(l_radian.shape[1]):
        x_CC[i,j], y_CC[i,j] = GalCoords_to_Moll(l_radian[i,j], b_radian[i,j])


P_max_CC = np.max(dP_dOmega_CC)

print ("P_max_CC = %.2e" % (P_max_CC))

cs = ax2.contourf(x_CC,y_CC,dP_dOmega_CC,levels=levs_CC,colors=colors_cc,alpha=alpha_transp)

cs = ax2.imshow



ax2.plot(x_ceq, y_ceq, 'y-', linewidth = 3)

ax2.set_xticklabels([150,120,90,60,30,0,-30,-60,-90,-120,-150])
ax2.set_yticklabels([-75,-60,-45,-30,-15,0.,15,30,45,60])


#####

P_SNR = kde.gaussian_kde([long_SNR, lat_SNR])

GC = np.array([0.,0.])
GACp = np.array([180.,0.])*np.pi/180.
GACm = np.array([180.,0.])*np.pi/180.
print ("P_SNR(GC,GACp,GACm) = (%.2e, %.2e, %.2e)" % (P_SNR(GC),P_SNR(GACp), P_SNR(GACm)))

PlotSNRKDE = False
PlotSNRKDE = True

if (PlotSNRKDE):
    xi, yi = np.mgrid[180:-180:360*4*1j, -90:90:180*4*1j]*np.pi/180. #Resolution of image
    #print ("KDE")
    #print (xi[0,0], yi[0,0], xi[-1,-1], yi[-1,-1])
    zi = P_SNR(np.vstack([xi.flatten(), yi.flatten()]))

    xMi = 0. * xi
    yMi = 0. * yi

    for i in range(xi.shape[0]):
        for j in range(xi.shape[1]):
            xMi[i,j], yMi[i,j] = GalCoords_to_Moll(xi[i,j], yi[i,j])




    ### uncomment this to make plot
    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.binary)

    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greys)
    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.jet)

    z_max = np.max(zi)

    print ("z_max = ",z_max)
    
    alpha_95 = 4.
    alpha_68 = 1.25
    P_lim_95 = z_max * np.exp(-alpha_95)
    P_lim_68 = z_max * np.exp(-alpha_68)

    levs_SNR = ([P_lim_95,P_lim_68])
    colors_SNR = (['gray','black'])
    linesty_SNR = (['dashed','solid'])

    plt.contour(xMi, yMi, zi.reshape(xMi.shape), levels=levs_SNR, colors=colors_SNR, linestyles=linesty_SNR)

    P_lim = P_lim_68


    N_SNR_tot = long_SNR.size

    N_SNR_in = 0
    for j in range(N_SNR_tot):
        P_obs = P_SNR([long_SNR[j], lat_SNR[j]])
        if (P_obs > P_lim):
            N_SNR_in = N_SNR_in + 1

    f_in = N_SNR_in/(N_SNR_tot+0.)

    print ("inside P_0/P_max = %.3f contour: %i of %i SNRs for f_in = %.3f" % (P_lim/z_max,N_SNR_in,N_SNR_tot,f_in))

    print ("Historical Ia")
    N_Ia_in = 0
    for j in range(hist_Ia_l.size):
        P_obs = P_SNR([hist_Ia_l[j], hist_Ia_b[j]])
        if (P_obs > P_lim):
            N_Ia_in = N_Ia_in + 1
            print ("%s: in" % (hist_Ia_names[j]))
        else:
            print ("%s: out" % (hist_Ia_names[j]))
            
    print ("Historical CC")
    N_CC_in = 0
    for j in range(hist_CC_l.size):
        P_obs = P_SNR([hist_CC_l[j], hist_CC_b[j]])
        if (P_obs > P_lim):
            N_CC_in = N_CC_in + 1
            print ("%s: in" % (hist_CC_names[j]))
        else:
            print ("%s: out" % (hist_CC_names[j]))
            



############################################################

x_SNR = 0. * long_SNR
y_SNR = 0. * long_SNR

for j in range(long_SNR.size):
    x_SNR[j], y_SNR[j] = GalCoords_to_Moll(long_SNR[j], lat_SNR[j])

plt.scatter(x_SNR, y_SNR, facecolors='k', edgecolors='black', zorder=5,alpha=0.5, marker=".")



plt.scatter(x_hIa, y_hIa, facecolors=Ia_color, edgecolors='y', zorder=10, marker=(5,1), s=200)
for j in range(3):
    plt.annotate(hist_Ia_names[j], (x_hIa_lab[j],y_hIa_lab[j]), \
                xytext=(x_hIa_tx[j],y_hIa_ty[j]), color='r', zorder=10, fontsize=15)
    
plt.scatter(x_hCC, y_hCC, facecolors='b', marker=(5,1), zorder=10, s=200, edgecolors='y')
for j in range(2):
    plt.annotate(hist_CC_names[j], (x_hCC_lab[j],y_hCC_lab[j]), \
                xytext=(x_hCC_tx[j],y_hCC_ty[j]), color='darkblue', zorder=10, fontsize=15)



plt.scatter(x_mCC, y_mCC, facecolors=maybe_CC_color, marker=(5,1), zorder=10, s=200, edgecolors='y')
for j in range(3):
    plt.annotate(maybe_CC_names[j], (x_mCC_lab[j],y_mCC_lab[j]), \
                xytext=(x_mCC_tx[j],y_mCC_ty[j]), color='cornflowerblue', zorder=10, fontsize=15)

    
plt.tick_params(axis='both', which='major', labelsize=15)


plt.title("Milky Way Supernova Remnants", weight='bold', fontsize=20)

#plt.text(20.*np.pi/180., 67.*np.pi/180., r"celestial equator", color='y', fontsize=15)

plt.xlabel(r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold")

plt.ylabel(r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold")

plt.text(0., 1.25*y_min, r"Galactic longitude ${\ell}$ [deg]",fontsize=20,weight="bold",horizontalalignment="center")
plt.text(1.15*x_max, 0., r"Galactic latitude ${b}$ [deg]",fontsize=20,weight="bold",verticalalignment="center",rotation=90)

legloc = ([0.05,0.03])
ax1.legend(handles=[cc68patch,cc95patch,ia68patch,ia95patch], loc=legloc, framealpha=1.)



l_NCP = 122.932
b_NCP = 27.1284

l_SCP = 302.932-360.
b_SCP = -27.1284


plt.text(x_NCP, y_NCP, r"NCP", color='y', fontsize=15, multialignment='center')
plt.text(x_SCP, y_SCP, r"SCP", color='y', fontsize=15, multialignment='center')

ylimits = np.array([-45.,+45.]) * (np.pi/180.0)
#plt.ylim(ylimits)



cc68patch = mpatches.Patch(color=colors_cc[1], label='CC 68%')
cc95patch = mpatches.Patch(color=colors_cc[0], label='CC 95%')
ia68patch = mpatches.Patch(color=colors_ia[1], label='Ia 68%')
ia95patch = mpatches.Patch(color=colors_ia[0], label='Ia 95%')
ceqline = Line2D([0], [0], color='y', linewidth=3, linestyle='-',label=r'Celestial Equator')
SNR68line = Line2D([0], [0], color='black', linewidth=3, linestyle='solid',label=r'SNRs 68%')
SNR95line = Line2D([0], [0], color='gray', linewidth=3, linestyle='dashed',label=r'SNRs 95%')
plt.rc('legend', fontsize=16)
legloc = ([0.05,0.01])
ax2.legend(handles=[cc68patch,cc95patch,ia68patch,ia95patch,SNR68line,SNR95line,ceqline], loc=legloc, framealpha=1.)


figbasename = "SNR_sky_Moll"

figname_png = figbasename+".png"
fig2.savefig(figname_png)
print ("plot written to: %s" % figname_png)

figname_eps = figbasename+".eps"
fig2.savefig(figname_eps)
print ("plot written to: %s" % figname_eps)

