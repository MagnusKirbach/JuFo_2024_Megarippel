import numpy as np
import matplotlib.pyplot as plt
import math as math
from mpl_toolkits.mplot3d import Axes3D
import random as rd
from matplotlib.collections import PolyCollection
from scipy.signal import find_peaks, savgol_filter, peak_prominences, peak_widths

def heaviside(x,flip):
	if x < flip:
		return 1
	else:
		return 0

Anzahl_durchläufe = 50
frequenz = 100
d_max = 8 # mm
t_max = 200/frequenz * Anzahl_durchläufe # s
delta_t = 0.4 # Schrittweite in s
h = 0.01 # Schrittweite in mm

standardabweichung = 0.5 # mm
erwartungswert = 1  # mm

d_min = 0.01 # mm
d_tau_0 = 3 # mm
d_tau_streuung = 0.75 # mm

n_d = int(round(d_max/h))
n_t = int(round(t_max/delta_t))
d_hilfe = np.linspace(d_min,d_max,n_d) # Korngrößen
t_hilfe = np.linspace(0,t_max,n_t) # Zeit

d_tau = np.zeros(n_t) # mm
Schnelligkeit = 1

R = np.zeros([n_t,n_d]) # ist die Erosions-Funktion, korngrößenabhängig, kleinere Körner werden schneller transportiert

mein_sinus = np.linspace(0,2*np.pi,n_t)
A = 1 # Ist die Konstante, mit der die Heaviside-Funktion multipliziert wird

for zeit in range(n_t):
    d_tau[zeit] = d_tau_0 + d_tau_streuung * np.sin(Anzahl_durchläufe*mein_sinus[zeit])
    for grain in range(n_d):
        R[zeit,grain] = np.sqrt(A)*(1-d_hilfe[grain]/d_tau[zeit])*heaviside(d_hilfe[grain],d_tau[zeit])
        
d,t = np.meshgrid(d_hilfe,t_hilfe)
kgv_oberfläche = np.zeros([n_t,n_d])

sigma = np.sqrt(math.log(1+(standardabweichung/erwartungswert)**2))
mü = math.log((erwartungswert**4)/(erwartungswert**2+standardabweichung**2))/2

log_wahrscheinlichkeit = np.zeros(n_d)

for i in range(n_d):
    log_wahrscheinlichkeit[i] = (1/(d_hilfe[i]*sigma*np.sqrt(2*np.pi)))*np.exp(-((np.log(d_hilfe[i])-mü)**2)/(2*sigma**2))

kgv_boden = log_wahrscheinlichkeit
R_mittel = 0
summe = 0
summe_danach = 0
puffer = 0.01

for k in range(n_d):
    kgv_oberfläche[0,k] = log_wahrscheinlichkeit[k]

for time in range(n_t-1):

    for korn in range(n_d):
        R_mittel += R[time,korn]*kgv_oberfläche[time,korn]*h

    for korn in range(n_d):

        kgv_oberfläche[time+1,korn] = kgv_oberfläche[time,korn] + delta_t * ( -R[time,korn]*kgv_oberfläche[time,korn] + R_mittel*kgv_boden[korn] )
        
        summe += kgv_oberfläche[time+1,korn]*h
        
    for j in range(n_d):
        kgv_oberfläche[time+1,j] += (1-summe)/(n_d*h)
        summe_danach += kgv_oberfläche[time+1,j]*h
    
    if ((summe_danach >= 1 - puffer) and (summe_danach <= 1 + puffer)): 
    	pass

    else:
        print("Summe danach:",summe_danach)
        	
    summe = 0
    summe_danach = 0
    R_mittel = 0

for x in range(n_t):
    for y in range(n_d):
        if x%99 == 0:
            pass
        else:
            kgv_oberfläche[x,y] = 0

def polygon_under_graph(x, y):
   return [(x[0], 0.), *zip(x, y), (x[-1], 0.)]

x = d[1,:] 
n_ts=round(n_t*delta_t)
times = range(0, n_ts)

verts = [polygon_under_graph(x, kgv_oberfläche[l*round(n_t/n_ts), :]) for l in times]

poly = PolyCollection(verts, alpha=.8, edgecolors='k')

ax = plt.figure(dpi = 200).add_subplot(projection='3d')
ax.add_collection3d(poly, zs=times, zdir='y')
ax.set(xlim=(0, d_max), ylim=(n_ts, 0), zlim=(0, 1.1), 
      xlabel='d in mm', ylabel='t in s', zlabel='p in 1/m')

plt.show()
