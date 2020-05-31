from math import pi
from kivy.app import App
from kivy.lang import Builder
from kivy.uix.floatlayout import FloatLayout
from kivy.garden.mapview import MapMarker, MarkerMapLayer
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.properties import ObjectProperty
from kivy.uix.popup import Popup
import matplotlib
matplotlib.use("module://kivy.garden.matplotlib.backend_kivy")
from kivy.garden.matplotlib.backend_kivyagg import FigureCanvasKivyAgg 
import matplotlib.pyplot as plt
import metody
import os

class LoadDialog(FloatLayout):
    
    load   = ObjectProperty(None)
    cancel = ObjectProperty(None)   
    
class Okno(FloatLayout):
    pass       

class WindowManager(ScreenManager):
    pass

class WiecejWykresow(Screen):
    wykresy = ObjectProperty()
    
    def wykres(self):
        filename = "nowyplik.gpx" 
        self.fig = metody.wykres(filename)
        self.cnv = FigureCanvasKivyAgg(self.fig)
        self.wykresy.add_widget(self.cnv)
        self.cnv.draw()


class WiecejStatystyk(Screen):
    tab = ObjectProperty()
    
    def statystyka(self):
        
        with open("INF2_Proj3.txt", "r") as plik:
            self.tab.text = "{}".format(plik.read())
        

class MainWindow(Screen):
    
    my_map     = ObjectProperty()
    data_layer = ObjectProperty()
    tekst      = ObjectProperty()
    wykres     = ObjectProperty()
    file       = ObjectProperty()

    def __init__(self, **kwargs):
        super(MainWindow, self).__init__(**kwargs)
        self.my_map.map_source = "thunderforest-landscape"
        
    def generuj_trase(self):
        
        if self.file is None:
            self.show_popup()
            
        else:
            dane = metody.wczytanie_pliku(self.file)
            self.draw_route(dane[0]) 
            self.rysuj_wykres(dane) 
            self.tekst.text = dane[1] 
            fi_sr = dane[0]["lat"].mean()*180/pi
            la_sr = dane[0]["lon"].mean()*180/pi
            self.my_map.center_on(fi_sr, la_sr)
            self.my_map.zoom = 12


        
    def rysuj_wykres(self, dane):
        self.fig = plt.figure() 
        self.ax1 = self.fig.add_subplot(111)
        self.ax1.set_title("Trasa") 
        self.ax1.scatter([dane[0]["lon"]*180/pi], [dane[0]["lat"]*180/pi], s=6) #rysowanie trasy
        self.cnv = FigureCanvasKivyAgg(self.fig) 
        self.wykres.add_widget(self.cnv)
        self.cnv.draw() #rysowanie wykresu
      
    #rysowanie trasy   
    def draw_route(self, df):
        self.data_layer = MarkerMapLayer() #utworzenie warstwy, creating instance of class MarkerMapLayer 
        self.my_map.add_layer(self.data_layer) #dodanie warstwy do mapy
        
        for i in range(df.shape[0]):
            self.mark_point(lat = float((df.iloc[i, 0])*180/pi), lon = float(df.iloc[i, 1]*180/pi), layer=self.data_layer) #narysowanie trasy
            
    #rysowanie znacznika - w tym przypadku kropki       
    def mark_point(self, lat, lon, layer=None, markerSource="dot.png"):
        
        if lat != None and lon != None:
            marker = MapMarker(lat=lat, lon=lon, source=markerSource ) #utworzenie markera
            self.my_map.add_marker(marker, layer=layer) #dodanie markera do mapy
    


    def usun_trase(self):
        
        if self.file is None:
            self.show_popup()
            
        else:
            self.my_map.remove_layer(self.data_layer) #usuniecie trasy
        
    def dismiss_popup(self):
        self._popup.dismiss()
        
    def show_popup(self):
        show = Okno() 
        
        # Create the popup window
        popupWindow = Popup(title="Blad!", content=show, size_hint=(None,None),size=(400,200))
        
                
        popupWindow.open() # show the popup 
                
    def load(self, path, filename):
        self.file = filename[0] 
        
        with open(os.path.join(path, filename[0])) as stream:
            with open('nowyplik.gpx', 'w+') as plik:
                plik.write(stream.read())
                
        self.dismiss_popup()
        
    

    def show_load(self):
        content = LoadDialog(load=self.load, cancel=self.dismiss_popup)
        self._popup = Popup(title="Wczytanie pliku", content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()


kv = Builder.load_file("main.kv")
sm = WindowManager()

screens = [MainWindow(name="main"), WiecejWykresow(name="ww"), WiecejStatystyk(name="ws")]
for screen in screens:
    sm.add_widget(screen)
    
sm.current = "main"

class MyApp(App):
    def build(self):
        return sm

        
if __name__ == '__main__':
    MyApp().run()