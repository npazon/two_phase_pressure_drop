# -*- coding: utf-8 -*-
"""
Created on Thu Apr  15 18:01:25 2021

@author: Nejka Pazon
"""

import numpy as np

def dokumentacija():
    dok = open('dok_homogeni.txt', 'r', encoding='utf8')
    return print(dok.read())

def lastnosti_toka(alpha, m_tot, rho_L, mi_L, rho_G, mi_G, D):
    '''Vrne slovar z lastnostmi toka v osnovnih merskih enotah,
    srednjo gostoto v kg/m^3, srednjo dinamično viskoznost v kg/(m*s)
    in Reynoldsovo število.
    
    - alpha je delež faze   
    - m_tot je celoten masni pretok homogenega toka v kg/(m^2*s)    
    - rho_L in rho_G sta gostoti kapljevite in plinaste faze v kg/m^3   
    - mi_L in mi_G sta dinamični viskoznosti kapljevite in plinaste faze v kg/(m*s)  
    - D je premer cevi v m'''
    A = np.pi*D**2/4 # konstanten pretočni presek cevi [m2]
    x = alpha*(rho_G/rho_L)/(1-alpha+alpha*(rho_G/rho_L)) # kvaliteta toka
    rho = (rho_L*rho_G)/(x*rho_L+(1-x)*rho_G) # srednja gostota [kg/m3]
    mi = (mi_L*mi_G)/(x*mi_L+(1-x)*mi_G) # srednja dinamična viskoznost (McAdamsova korelacija) [kg/ms]
    w = m_tot/rho # hitrost toka v cevi [m/s]
    w_G = w_L = w # predpostavka homogenega toka
    (A_G, A_L) = (alpha*A, (1-alpha)*A) # pretočni presek posamezne faze [m2]
    (m_G, m_L) = (rho_G*w_G*A_G, rho_L*w_L*A_L) # masni pretok plinaste in kapljevite faze [kg/s]
    V_tot = (m_G+m_L)/rho # volumski pretok [m3/s]
    beta = alpha # delež volumskega toka
    (V_G, V_L) = (beta*V_tot, (1-beta)*V_tot) # volumski pretok plinaste in kapljevite faze [m3/s]
    (j_G, j_L) = (V_G/A, V_L/A) # povprečena hitrost plinaste in kaplevite faze [m/s]
    Re = m_tot*D/mi #Reynoldsovo število

    lastnosti = {'rho':rho, 'mi':mi, 'A':A, 'A_G':A_G, 'A_L':A_L,
    'x':x,'w':w, 'm_G':m_G, 'm_L':m_L, 'j_G':j_G, 'j_L':j_L, 'Re':Re}

    return lastnosti, rho, mi, Re
    
def f_Darcy(D, e, Re, tokovni_vzorec=False):
    '''Vrne Dracyjev faktor trenja.
    V primeru trubulentnega toka faktor trenja izračuna po Colebrookovi enačbi.
    V primeru mehurčkastega ali rahlo obročastega toka izračuna faktor trenja po korekcijski enačbi.
         
    - D je premer cevi v m
    - e je hrapavost stene cevi v m   
    - Re je Reynoldsovo število    
    - tokovni_vzorec = True pomeni mehurčkasti ali rahlo obročasti tok'''
    if Re<2300:
        f_D = 64/Re
    else:
        f_0 = 64/2300
        f_D = 0
        err = 10**(-7)
        err_i = np.abs(f_0 - f_D)
        i = 1
        while np.abs(err_i >= err):
            f_D = (-2.0*np.log10((e/D)/3.7 + 2.51/(Re*np.sqrt(f_0))))**(-2)
            err_i = np.abs(f_0 - f_D)
            i = i + 1
            f_0 = f_D
    if tokovni_vzorec == True: # korekcija faktorja trenja za mehurčkasti ali rahlo obročasti tok
        f_0 = f_D/4
        err = 10**(-7)
        f_F = 0
        err_i = np.abs(f_0 - f_F)
        i = 1
        while np.abs(err_i >= err):
            f_F = (3.48-4*np.log10(2*(e/D) + 9.35/(Re*np.sqrt(f_0))))**(-2)
            err_i = np.abs(f_0 - f_F)
            i = i + 1
            f_0 = f_F
        f_D = f_F*4

    return f_D

def tlačne_izgube(m_tot, rho, theta, D, f_D):
    '''Vrne tlačne izgube zaradi trenja in gravitacijskega pospeška v cevi v Pa/m.
    
    - m_tot je celoten masni pretok homogenega toka v kg/(m^2*s)  
    - rho je srednja gostota homogenega toka v kg/m^3   
    - theta je kot, ki ga cev zaklepa z vertikalo v rad     
    - D je premer cevi v m     
    - f_D je Darcyjev faktor trenja'''
    g = 9.81 # gravitacijski pospešek [m/s2]
    return f_D * m_tot**2 / (2 * D *rho) + rho * g * np.cos(theta)

def tlačni_padec(dp, L):
    '''Vrne tlačni padec v cevi na dolžini L v Pa.

    - dp_F so tlačne izgube zaradi trenja v cevi v Pa/m    
    - dp_G so tlačne izgube v cevi zaradi gravitacijskega pospeška v Pa/m'''
    return dp * L