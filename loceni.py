# -*- coding: utf-8 -*-
"""
Created on Thu Apr  15 18:01:25 2021

@author: Nejka Pazon
"""

import numpy as np

def dokumentacija():
    dok = open('dok_ločeni.txt', 'r', encoding='utf8')
    return print(dok.read())

def lastnosti_toka(beta, m_tot, rho_L, mi_L, rho_G, mi_G, D):
    '''Vrne srednjo gostoto v kg/m^3, povprečeni hitrosti posameznih faz v m/s,
    in Reynoldsovi števili za samo eno od faz v cevi.
    
    - beta je delež volumskega pretoka   
    - m_tot je masni pretok dvofaznega toka v kg/s   
    - rho_L in rho_G sta gostoti kapljevite in plinaste faze v kg/m^3   
    - mi_L in mi_G sta dinamični viskoznosti kapljevite in plinaste faze v kg/(m*s)  
    - D je premer cevi v m'''
    A = np.pi*D**2/4 # konstanten pretočni presek cevi [m2]
    rho = beta*rho_G+(1-beta)*rho_L # srednja gostota [kg/m3]
    V_tot = m_tot/rho # volumski pretok [m3/s]
    (V_G, V_L) = (beta*V_tot, (1-beta)*V_tot) # volumski pretok plinaste in kapljevite faze [m3/s]
    (m_G, m_L) = (V_G*rho_G, V_L*rho_L) # masni pretok plinaste in kapljevite faze [kg/s]
    (j_G, j_L) = (V_G/A, V_L/A) # povprečena hitrost plinaste in kapljevite faze [m/s]
    (Re_G, Re_L) = (j_G*rho_G*D/mi_G, j_L*rho_L*D/mi_L) #Reynoldsovo število samo kapljevite in samo plinaste faze v cevi 

    lastnosti = {'rho':rho, 'A':A, 'V_tot':V_tot, 'm_G':m_G, 'm_L':m_L,
    'V_G':V_G, 'V_L':V_L, 'j_G':j_G, 'j_L':j_L, 'Re_G':Re_G, 'Re_L':Re_L}
    
    return lastnosti, rho, j_G, j_L, Re_G, Re_L

def f_Darcy_faza(D, e, Re):
    '''Vrne Dracyjev faktor trenja.
    V primeru trubulentnega toka faktor trenja izračuna po Colebrookovi enačbi.
         
    - D je premer cevi v m
    - e je hrapavost stene cevi v m   
    - Re je Reynoldsovo število    '''
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
    
    return f_D

def tlačne_izgube_faza(j, rho, theta, D, f_D):
    '''Vrne tlačne izgube zaradi trenja in gravitacijskega pospeška v primeru le ene faze v cevi v Pa/m.
    
    - j je povprečena hitrost faze v m/s 
    - rho je gostota faze v kg/m^3   
    - theta je kot, ki ga cev zaklepa z vertikalo v rad     
    - D je premer cevi v m     
    - f_D je Darcyjev faktor trenja za fazo'''
    g = 9.81 # gravitacijski pospešek [m/s2]
    return f_D * rho * j**2 / (2 * D) + rho * g * np.cos(theta)

def tlačne_izgube_LM(dpf_G, dpf_L, Re_G, Re_L):
    '''Vrne tlačne izgube zaradi trenja v cevi po Lockhart-Martinelli korelaciji v Pa/m.
    
    - dpf_G in dpf_L so tlačne izgube zaradi trenja v primeru le ene faze v cevi v Pa/m
    - Re_G in Re_L sta Reynoldsovi števili v primeru le ene faze v cevi'''
    C = 0
    psi = (dpf_L/dpf_G)**(1/2) #Lockhart-Martinelli parameter
    if Re_L < 2300:
        if Re_G < 2300:
            C = 5
        elif Re_G > 2300:
            C = 12
    elif Re_L > 2300:
        if Re_G < 2300:
            C = 10
        elif Re_G > 2300:
            C = 20
    phi_G = 1+C*psi+psi**2 #multiplikator
    
    return dpf_G*phi_G

def tlačni_padec(dp_LM, L):
    '''Vrne tlačni padec v cevi na dolžini L v Pa.

    - dp_F so tlačne izgube zaradi trenja v cevi v Pa/m    
    - L je dolžina, na kateri nas zanima tlačni padec v cevi'''
    return dp_LM * L