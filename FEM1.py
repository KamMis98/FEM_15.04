# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 15:08:54 2021
@author: Kamil
"""
'''
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % PRE - PROCESSING % % % % %
%% 1(a) inicjalizacja parametrow sterujacych programem CZESC 1
[ parametry_sterujace ] = inicjalizacja_parametrow_sterujacych ( ... ) ;
%% 1(b) definicja parametrow fizycznych i geometrycznych obszaru , warunkow brzegowych CZESC 1
% - definicja przedzialu ,
% - liczba wezlow / elementow modelu w przypadku dyskretyzacji jednorodnej
% - zadanie funkcji wymuszajacej , parametrow rownania rozniczkowego
% - ...
[ parametry_geom_i_fiz ] = definicja_parametrow_geom_i_fiz ( ... ) ;
%% 1(b/c) CZESC 1
[ WEZLY , ELEMENTY , WAR_BRZEGOWE ] = definicja_zmiennych_przechowujacych_informacje_o_geometrii (
parametry_geom_i_fiz , ... )
% - / ew. odczytanie geometrii z pliku
%% 1(d) prezentacja geometrii zagadnienia CZESC 2
rysuj_geometrie ( WEZLY , ELEMENTY , WAR_BRZEGOWE ) ;
% 1(e) utworzenie macierzy wypelnionych zerami CZESC 3
[A , b] = alokacja_pamieci_na_zmienne_globalne ( liczba_wezlow ) ;
%% 2(a) CZESC 3
[ lokalne_fun_ksztaltu , pochodne_lokalnych_fun_ksztaltu ] = definicja_funkcji_ksztaltu ( ... ) ;
% % % % % PROCESSING % % % % %
for k = 1: liczba_elementow_skonczonych
%% 2(b) CZESC 4
[ nr_glob_elem , nr_glob_wezlow_elem , wspolrzedne_wezlow , jakobian ] =
zgromadzenie_informacji_dotyczacych_elementu (k , WEZLY , ELEMENTY , ...) ;
M = obliczenie_lokalnej_macierzy_opisujacej_element ( ) ;
%% 2(c) CZESC 5
A = agregacja_macierzy_lokalnej_w_macierzy_globalnej (M , nr_glob_wezlow_elem , ...) ;
b = obliczanie_elementow_wektora_prawej_strony ( ) ;
end % for k = 1: liczba_elementow_skonczonych
%% 2(d) CZESC 6
[A , b] = uwzglednienie_warunkow_brzegowych ( WAR_BRZEGOWE , ... ) ;
%% 2(e) CZESC 7
a = rozwiazanie_url (A ,b )
% % % % % POST - PROCESSING % % % % %
%% 3(a) obrobka_wynikow ???? CZESC 8
%% 3(b) prezentacja graficzna rozwiazania
rysuj_rozwiazanie ( WEZLY , ELEMENTY , a)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
'''
1. Nowe repozytorium

### PREPROCESSING
2. Skrypt gĹowny - inicjalizacja parametrĂłw:
    * Parametry sterujace ( nabla2 u + c u = f ) -> c=0, f=0 (stale -> docelowo -> funkcje lambda)   
    * Definicja geometrii:
        * parametry geometrii (poczatek/koniec przedzialu)    x_0 = 0; x_p = 1   
        * liczba wezlow   
        * tablice przechowujace informacje o wezĹach, elementach, warunkach brzegowych  (Tablice typu `array` z pakietu `numpy`)
        
            WEZLY
            
            |l.p. | x_i |
            |:---:|:---:|
            | 1   | 0   |
            | 2   | 1   |
            | 3   |0.5  |
            | 4   |0.75 |
         
            ELEMENTY           
            
            |l.p. | iwp | iwk |
            |:---:|:---:|:---:|
            |1    |  1  |  3  |
            |2    |  4  |  2  |
            |3    |  3  |  4  |

            # warunki brzegowe:
            
            twb_L = 'D' 
            twb_P = 'D'
            
            wwb_L = 0
            wwb_P = 1            
         
3. Funkcja generujaca geometriÄ w sposĂłb automatyczny
    
    n = 100
    
    WEZLY, ELEMENTY = generujTabliceGeometrii(x_0, x_p, n)

4. Funkcja przedstawiajÄ
ca geometriÄ na podstawie zadanych wezlow, elementĂłw...

5. Utworzenie zmiennych przechowujacych macierz wspĂłĹczynnikĂłw i wektor prawej strony ukĹadu.

6. ..... Lokalne funkcje ksztaĹtu... ? - funkcje lambda:
   - wer. 1 - funkcje liniowe
   - wer. 2 (trudniejsza) generowanie funkcji z bazy Lagrange'a dla zadanego k

### PROCESSING ...

ZADANIE 1 - ZROBIONE
'''
import numpy as np
import matplotlib.pyplot as plt

'''
ZADANIE 2
'''
#parametry sterujące
c=0
f=0

#geometria
x_a=0 #początek przedziału
x_b=1 #koniec przedziału
#x_k - węzeł końcowy
#x_p - węzeł początkowy
n=4   #liczba węzłów, numer węzła; x_i - współrzędna węzła

#tablice
wezly=np.array([[1, 0],
               [2, 1],
               [3, 0.5],
               [4, 0.75]])

elementy=np.array([[1, 1, 3],
                  [2, 4, 2],
                  [3, 3, 4]])

#warunki brzegowe
twb_L = 'D'
twb_P = 'D'
wwb_L = 0
wwb_P = 1

warunki_brzeg=np.array([x_a, x_b])

#def tablica_wezlow(x_i, n):
 #   temp = (x_a-x_b)/(n-1)
  #  wezly = np.array([wezly, elementy])

'''
ZADANIE 3
'''
def generujTabliceGeometrii(x_a, x_b, n):
     temp = (x_a-x_b)/(n-1)
     wezly = np.array([x_a, x_b])
     for i in range(1, n-1, 1):
       wezly = np.block([wezly, i*temp+x_a])
     return wezly

print(generujTabliceGeometrii(1,2,5))

def rysujGeometrie(wezly, elementy):
    ilosc_wezlow=wezly.size()  #sprawdź czy to ta funkcja da iloć wezłów
    wektor_zer=np.zeros()
    x=wezly[:,1] #x to druga kolumna tablicy węzły
   # y=elementy(:,)
    plt.plot( wezly[:,1], wektor_zer,)
    plt.plot(x,y,"r:",linewidth=6)

'''
ZADANIE 4
'''
for i in range(1, n):
        print("-")
        print("|")
        
#print("|")



