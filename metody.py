import gpxpy
import pandas as pd
from math import sin, cos, sqrt, tan, pi, atan
from datetime import timedelta
import matplotlib.pyplot as plt


def Vincenty(fi_a,la_a,fi_b,la_b,a = 6378137.000, e2 = 0.00669438002290):
    """
    #INPUT:
         fi_a [rad]
         la_a [rad]
         fi_b [rad]
         la_b [rad]
    #OUTPUT:
        s_AB [m]
    """
    
    if fi_a == fi_b and la_a == la_b:
        s_AB = 0
        
    elif fi_a == fi_b and la_a != la_b: 
        N = a/sqrt(1-e2*(sin(fi_a))**2)
        dL = abs(la_a - la_b)  
        s_AB = dL*N*cos(fi_a)
        
    else:
        b = a*sqrt(1-e2)
        f = 1 - (b/a)
        delta_la = la_b - la_a
        U_a = atan((1-f)*tan(fi_a))
        U_b = atan((1-f)*tan(fi_b))
        L = delta_la
        
        while True:
            sin_sigma = sqrt((cos(U_b)*sin(L))**2 + (cos(U_a)*sin(U_b) - sin(U_a)*cos(U_b)*cos(L))**2)
            cos_sigma = sin(U_a)*sin(U_b) + cos(U_a)*cos(U_b)*cos(L)
            sigma = atan(sin_sigma/cos_sigma)
            sin_alfa = (cos(U_a)*cos(U_b)*sin(L))/(sin_sigma)
            cos2_alfa = 1 - (sin_alfa)**2
            cos2_sigma_m = cos_sigma - (2*sin(U_a)*sin(U_b))/(cos2_alfa)
            C = (f/16)*cos2_alfa*(4+f*(4-3*cos2_alfa))
            Ls = L
            L = delta_la + (1-C)*f*sin_alfa*(sigma + C*sin_sigma*(cos2_sigma_m + C*cos_sigma*(-1+2*(cos2_sigma_m)**2)))
            
            if (L-Ls)<(0.000001/206265):
                break
            
        u2 = ((a**2 - b**2)/b**2)*cos2_alfa
        A = 1 + (u2/16384)*(4096 + u2*(-768 + u2*(320 - 175*u2)))
        B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
        delta_sigma = B*sin_sigma*(cos2_sigma_m + (1/4)*B*(cos_sigma*(-1 + 2*(cos2_sigma_m)**   2) - (1/6)*B*cos2_sigma_m*(-3 + 4*(sin_sigma)**2)*(-3 + 4*(cos2_sigma_m)**2)))
        s_AB = b*A*(sigma - delta_sigma)

    return s_AB

def wczytanie_pliku(plik):
    """
    #INPUT:
        plik z rozszerzeniem gpx
    #OUTPUT:
        df[szerokosc[rad], dlugosc[rad], h[m], czas[data, h:min:sec]]
        suma_odl3D [m]
        czas [h:min:sec]
        suma_delta_h_dod [m]
        suma_delta_h_ujm [m]
        suma_delta_h [m]
        sr_predkosc [m/s]
        sr_nachylenie [%]     
        zapis do pliku 'INF2_Proj3.txt'
    """
    lon = []
    lat = []
    el = []
    dates = []
    
    with open(plik, "r") as gpx_file: #otworzenie pliku
        gpx = gpxpy.parse(gpx_file)
        
    for track in gpx.tracks:
        for seg in track.segments:
            for point in seg.points:
                lon.append(point.longitude)
                lat.append(point.latitude)
                el.append(point.elevation)
                
                if point.time !=None:
                    point.time = point.time.replace(tzinfo=None)
                    dates.append(point.time)
                    
                else:
                    dates = []
                    
    if dates == []:
        df = pd.DataFrame({"lat":lat, "lon":lon, "el":el})
        
    else:    
        df = pd.DataFrame({"lat":lat, "lon":lon, "el":el, "time":dates})
        
    df["lat"] = df["lat"]*pi/180
    df["lon"] = df["lon"]*pi/180
    df.dropna(subset=["el"], axis=0, inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    if 'time' in df.columns:
        df.dropna(subset=["time"], axis=0, inplace=True)
        
    else:
        pass

    df.reset_index(drop=True, inplace=True)    
        
    suma_odl2D = 0
    suma_odl3D = 0
    suma_seconds = 0
    suma_delta_h = 0
    suma_delta_h_dod = 0
    suma_delta_h_ujm = 0
    t = []
    odl = []
    v = []
    
    with open('INF2_Proj3.txt', 'w+') as plik:
        plik.write(88*"-")
        plik.write("\n|{:^11}|{:^10}|{:^10}|{:^10}|{:^10}|{:^10}|\n".format("Odcinek", "Odleglosc [m]", "Przewyzszenie [m]", "Czas [sec]", "Predkosc [km/h]", "Nachylenie [%]"))
            
        for i in df.index.values:
            if i + 1 > max(df.index.values):
                break
            else:
                odl_2D = Vincenty(df.iloc[i, 0], df.iloc[i, 1], df.iloc[i + 1, 0], df.iloc[i + 1, 1]) #odleglosc 2D
            #       ee.append(odl)
            suma_odl2D += odl_2D # suma odleglosci 2D
            delta_h = df.iloc[i + 1, 2] - df.iloc[i, 2] #przewyzszenie
            
            if delta_h > 0:
                suma_delta_h_dod += delta_h
            elif delta_h < 0:
                suma_delta_h_ujm += delta_h
            else:
                pass
            
            suma_delta_h += delta_h
            odl_3D = sqrt(odl_2D**2 + delta_h**2) # odleglosc 3D
            
            
            if odl_3D == 0:
                nachylenie = 0
            else:
                nachylenie = round((delta_h/odl_3D)*100, 3)
                
            suma_odl3D += odl_3D # suma odleglosci 3D
            odl.append(suma_odl3D)
            
            if 'time' in df.columns:
                delta_t = df.iloc[i + 1, 3] - df.iloc[i, 3]
                seconds = delta_t.total_seconds() # czas w sekundach
                suma_seconds += seconds # suma sekund
                czas = timedelta(seconds=suma_seconds) # czas godz:min:sec
                t.append(suma_seconds/60)
                plik.write(88*"-")
                
                if seconds != 0:
                    predkosc = (odl_3D/seconds)*3.6
                    v.append(predkosc)
                    plik.write("\n|{:^5}-{:^5}|{:^13}|{:^17}|{:^10}|{:^14}|{:^14}|\n".format(i, i + 1, round(odl_3D, 3), round(delta_h, 3), seconds, round(predkosc,3), nachylenie))
                    
                else:
                    predkosc = 0 #439 dla krk1 ???
                    v.append(predkosc)
                    plik.write("\n|{:^5}-{:^5}|{:^13}|{:^17}|{:^10}|{:^14}|{:^14}|\n".format(i, i + 1, round(odl_3D, 3), round(delta_h, 3), seconds, predkosc, nachylenie))
            
                plik.write(88*"-")
                
            else:
                plik.write("\n|{:^5}-{:^5}|{:^13}|{:^17}|{:^10}|{:^14}|{:^14}|\n".format(i, i + 1, round(odl_3D, 3), round(delta_h, 3), "Brak inf", "Brak inf", nachylenie))
    
    if suma_seconds == 0:
        sr_predkosc = "Brak informacji"
        czas = "Brak informacji"
        
    else:
        sr_predkosc = round((suma_odl3D/suma_seconds)*3.6, 2)
        
    sr_nachylenie = round((suma_delta_h/suma_odl3D)*100, 2)
    pods = "Całkowita odległosć: {} m\nCzas: {}\nCałkowite przewyższenie:{} m\nPrzewyższenie w górę: {} m\nPrzewyższenie w dół: {} m\nŚrednia prędkosć: {} km/h\nŚrednie nachylenie: {} %".format(round(suma_odl3D, 3), czas, round(suma_delta_h, 3), round(suma_delta_h_dod, 3), round(suma_delta_h_ujm, 3), sr_predkosc, sr_nachylenie)
    
    
    if 'time' in df.columns:
        dane = (df, pods, odl, t, v)
    else:
        dane = (df, pods, odl)
        
    return dane


def wykres(plik):
    dane = wczytanie_pliku(plik)
    x = list(dane[0].index.values)
    dane[4].insert(0, 0)
    graph = plt.figure(figsize=(30,15))
#    
    ax1 = graph.add_subplot(311)    
    ax1.plot(x,dane[0]["el"]) 
    ax1.set_ylabel('H [m]')
    
    ax2 = graph.add_subplot(312)      
    ax2.scatter(dane[3],dane[2], s=3) 
    ax2.set_ylabel('S [m]')
    
    ax3 = graph.add_subplot(313)
    ax3.plot(x, dane[4]) 
    ax3.set_ylabel('V [km/h]')
    plt.savefig('wykresy.png')

    return graph
    



















