import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spint

from GeoDef import *
from AutoGenGeo import *
from RysGeo import *
from Alok import *
from FunBaz import *
from Aij import *
from RysRozw import *

if __name__ == '__main__':
    
    # Preprocessing
    WEZLY, ELEMENTY, WB = GeoDef()
    n=np.shape(WEZLY)[0]

    c = 0 
    f = lambda x: 0*x # wymuszenie

    x_a =  0 
    x_b =  1
    n = 4
    
    
    WEZLY, ELEMENTY = AutoGenGeo(x_a,x_b,n)
    # warunki brzegowe
    WB    = [{"ind": 1, "typ":'D', "wartosc":1}, 
              {"ind": n, "typ":'D', "wartosc":2}]
    
    #rysujGeometrię
    RysGeo(WEZLY, ELEMENTY, WB)
    #alokacja
    A,b = Alok(n)

    stopien_fun_bazowych = 1
    phi, dphi = FunBaz(stopien_fun_bazowych)
    
    #PROCESSING
    
    liczbaElementow = np.shape(ELEMENTY)[0]
    
    for ee in np.arange(0, liczbaElementow ):
        
        elemRowInd = ee
        elemGlobalInd = ELEMENTY[ee,0]        
        elemWezel1 = ELEMENTY[ee,1]     
        elemWezel2 = ELEMENTY[ee,2]    
        indGlobalneWezlow = np.array([elemWezel1, elemWezel2 ])
    
        x_a = WEZLY[ elemWezel1-1 ,1]
        x_b = WEZLY[ elemWezel2-1 ,1]
    
    
        Ml = np.zeros( [stopien_fun_bazowych+1, stopien_fun_bazowych+1] )
        
        J = (x_b-x_a)/2
        
        m = 0; n = 0 ;
        Ml[m,n] = J * spint.quad( Aij(dphi[m], dphi[n], c, phi[m],phi[n]), -1, 1)[0]
        
        m = 0; n = 1 ;
        Ml[m,n] = J * spint.quad( Aij(dphi[m], dphi[n], c, phi[m],phi[n]), -1, 1)[0]
        
        m = 1; n = 0 ;
        Ml[m,n] = J * spint.quad( Aij(dphi[m], dphi[n], c, phi[m],phi[n]), -1, 1)[0]
        
        m = 1; n = 1 ;
        Ml[m,n] = J * spint.quad( Aij(dphi[m], dphi[n], c, phi[m],phi[n]), -1, 1)[0]
                
        
        
        
        
        A[np.ix_(indGlobalneWezlow-1, indGlobalneWezlow-1  ) ] =  \
        A[np.ix_(indGlobalneWezlow-1, indGlobalneWezlow-1  ) ] + Ml
        
    print(WB)
        
    
    # Uwzględnienie warunków brzegowych
    
    if WB[0]['typ'] == 'D':
        ind_wezla = WB[0]['ind']
        wart_war_brzeg = WB[0]['wartosc']
        
        iwp = ind_wezla - 1
        
        WZMACNIACZ = 10**14
        
        b[iwp] = A[iwp,iwp]*WZMACNIACZ*wart_war_brzeg
        A[iwp, iwp] = A[iwp,iwp]*WZMACNIACZ
        
        
    if WB[1]['typ'] == 'D':
        ind_wezla = WB[1]['ind']
        wart_war_brzeg = WB[1]['wartosc']
        
        iwp = ind_wezla - 1
        
        WZMACNIACZ = 10**14
        
        b[iwp] = A[iwp,iwp]*WZMACNIACZ*wart_war_brzeg
        A[iwp, iwp] = A[iwp,iwp]*WZMACNIACZ        
    
    
    if WB[0]['typ'] == 'N':
        print('Nie zaimplementowano jeszcze. Zad.dom')
    
    
    if WB[1]['typ'] == 'N':
        print('Nie zaimplementowano jeszcze. Zad.dom')    
    
    
    # Rozwiazanie ukladu
    
    u = np.linalg.solve(A,b)
    
    RysRozw(WEZLY, ELEMENTY, WB, u)
    
    
    
    
    