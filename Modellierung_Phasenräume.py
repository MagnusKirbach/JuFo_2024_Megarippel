import numpy as np
import math as math
import random as rd
from scipy.signal import find_peaks, savgol_filter

def heaviside(x,flip):
	if x < flip:
		return 1
	else:
		return 0

anzahl_durchläufe = 50

def daten_kgv_oberflaeche(d_tau_0,d_tau_streuung,standardabweichung,erwartungswert,r_tau,frequenz):
    
    d_max = 6 # mm
    t_max = 200/frequenz * anzahl_durchläufe # s
    h = 0.01 # Schrittweite in mm
    delta_t = 0.4 # Schrittweite in s    

    d_min = 0.01 # mm
    d_tau_0 = d_tau_0 # mm
    d_tau_streuung = d_tau_streuung # mm

    n_d = int(round(d_max/h))
    n_t = int(round(t_max/delta_t))
    d_hilfe = np.linspace(d_min,d_max,n_d) # Korngrößen
    t_hilfe = np.linspace(0,t_max,n_t) # Zeit
    
    d_tau = np.zeros(n_t)  # mm
    R = np.zeros([n_t,n_d]) # ist die Erosions-Funktion, korngrößenabhängig, kleinere Körner werden schneller transportiert
    mein_sinus = np.linspace(0,2*np.pi,n_t)
    A = r_tau

    for zeit in range(n_t):
        d_tau[zeit] = d_tau_0 + d_tau_streuung * np.sin(anzahl_durchläufe*mein_sinus[zeit])
        for grain in range(n_d):
            R[zeit,grain] = np.sqrt(A)*(1-d_hilfe[grain]/d_tau[zeit])*heaviside(d_hilfe[grain],d_tau[zeit])

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
        
        if summe_danach >= 1 - puffer and summe_danach <= 1 + puffer: 
            pass

        else:
            print("Summe danach:",summe_danach)
        summe = 0
        summe_danach = 0
        R_mittel = 0

    return kgv_oberfläche, n_t

frequenz = np.linspace(1,15,15)
EX = np.linspace(0.1,3,30)
DX = np.linspace(0.1,1,10)
d_tau_streuung = np.linspace(0,1.5,31)
d_tau_null = np.linspace(0.1,5.9,29)
r_tau = np.array([0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
linspace_weil = np.linspace(1,5,15)
r_tau = np.append(r_tau, linspace_weil)

speicher = np.zeros([len(frequenz), len(DX)])

for i in range(len(frequenz)):
    for j in range(len(DX)):
        funktion = daten_kgv_oberflaeche(3,0.75,DX[j],0.25,1,frequenz[i])
        rohe_daten = funktion[0]
        n_t = funktion[1]

        raw_ausgabe_peaks = find_peaks(rohe_daten[len(rohe_daten)-1,:],prominence=1e-03, width=4)
        raw_ausgabe_peaks_davor = find_peaks(rohe_daten[len(rohe_daten)-1-int(n_t/anzahl_durchläufe),:],  prominence=1e-03, width=4)
        
        letzte_peaks_hilfe = raw_ausgabe_peaks[0]
        davor_peaks_hilfe = raw_ausgabe_peaks_davor[0]

        if ((0.75 < np.abs(3-0.25)*0.7) and (0.25+DX[j] < 3-0.75)):
            speicher[i,j] = len(raw_ausgabe_peaks[0])
        else:
            speicher[i,j] = -len(raw_ausgabe_peaks[0]) 
        
        if len(letzte_peaks_hilfe) != len(davor_peaks_hilfe):
            print("nicht im Gleichgewicht, keine gleiche Anzahl von Peaks")
            speicher[i,j] = 8
        else:
            for k in range(len(letzte_peaks_hilfe)):
                if np.abs(letzte_peaks_hilfe[k]-davor_peaks_hilfe[k]) > 10:
                    print("nicht zwangsläufig im Gleichgewicht, Überprüfung notwendig")
                    speicher[i,j] = 8
        print(raw_ausgabe_peaks[0], DX[j], frequenz[i])

print(speicher)
print("Phasenraum mit Frequenz (x) und DX (y), 8 ungleichgewicht, - davor böses Regime, EX auf 0.25, d_tau_streeung auf 0.75, d_tau_0 auf 3, d_max auf 6, r auf 1 ")
speicher = np.savetxt("DX_EX_schlecht_schlecht.csv", speicher, fmt="%d", delimiter=",")
