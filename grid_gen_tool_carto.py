# BaseMap example by geophysique.be
# tutorial 08
import cartopy as Basemap
import matplotlib.pyplot as plt
import numpy as np
from math import sin, cos, sqrt, atan2, radians 

class Click():
    def __init__(self, ax, func, button=1):
        self.ax=ax
        self.func=func
        self.button=button
        self.press=False
        self.move = False
        self.c1=self.ax.figure.canvas.mpl_connect('button_press_event', self.onpress)
        self.c2=self.ax.figure.canvas.mpl_connect('button_release_event', self.onrelease)
        self.c3=self.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmove)

    def onclick(self,event):
        if event.inaxes == self.ax:
            if event.button == self.button:
                self.func(event)
    def onpress(self,event):
        self.press=True
    def onmove(self,event):
        if self.press:
            self.move=True
    def onrelease(self,event):
        if self.press and not self.move:
            self.onclick(event)
        self.press=False; self.move=False

def calc_dist(glon1,glat1,glon2,glat2):
    # Distance Calculation
    lat1 = radians(glat1)
    lon1 = radians(glon1)
    lat2 = radians(glat2)
    lon2 = radians(glon2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    azimuth = atan2(dlon,dlat)*(180./np.pi)

    distance = R * c
    nlats = int(np.ceil(distance/grid_res))
    return distance, azimuth

 
def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    # Radian conversion
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.

    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
def great(m, startlon, startlat, azimuth,*args, **kwargs):
    glon1 = startlon
    glat1 = startlat
    glon2 = glon1
    glat2 = glat1
 
    step = 50
 
    glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    if azimuth-180 >= 0:
        while glon2 <= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    elif azimuth-180 < 0:
        while glon2 >= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)

def grid_gen(event):
    global MEEP, lat1, lon1, lats, lons, azimuth

    # PLOT FIRST POINT
    if MEEP<1:
       lon1 = event.xdata; lat1 = event.ydata
       print(lon1,lat1)
       m.plot(lon1,lat1,'bo')
       plt.draw()
       MEEP+=1

    # CALCULATE AND PLOT SECOND ALONG  GC LINE
    elif MEEP<2:
       distance, azimuth = calc_dist(lon1,lat1,event.xdata,event.ydata)
       nlats = int(np.ceil(distance/grid_res))
       lats = np.zeros((nlats,1))
       lons = np.zeros((nlats,1))

       for nt in range(nlats):
           lons[nt], lats[nt], baz = shoot(lon1, lat1, azimuth, nt*grid_res)
           if m.is_land(lons[nt],lats[nt]):
              m.plot(lons[nt],lats[nt],'go',zorder=3)
           else:
              m.plot(lons[nt],lats[nt],'bo',zorder=3)

       m.drawgreatcircle(lons[0], lats[0], lons[-1], lats[-1],del_s=10,color='k', lw=2.)
       plt.draw()
       MEEP+=1
    
    # CALCULATE SHORTEST PERP. DISTANCE
    elif MEEP<3:
         min_dist = 9999
         for nl in range(len(lats)):
             distance,_ = calc_dist(lons[nl],lats[nl],event.xdata,event.ydata)
             if distance<min_dist:
                min_dist=distance
         nlons = int(np.ceil(min_dist/grid_res))

         lats = np.tile(lats,nlons)
         lons = np.tile(lons,nlons)

         for nt in range(len(lats)):
             for nn in range(nlons):
                 # Already have at initial line
                 if nn>0:
                    lons[nt,nn], lats[nt,nn], baz = shoot(lons[nt,0], lats[nt,0], azimuth+90, nn*grid_res)

         m.plot(lons,lats,'bo')
         plt.draw() 
         # SAVE AS 2D ARRAY... pre land-masking  
         np.savez('test_grid', lats, lons)
         MEEP+=1
MEEP=0
# approximate radius of earth in km
R = 6378.13
 
fig = plt.figure(figsize=(11.7,8.3))
plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
ax = plt.subplot(111)
 
m = Basemap(resolution='l')
m.drawmapboundary(fill_color='azure')
m.fillcontinents(color='palegoldenrod',lake_color='azure')
m.drawcoastlines()

grid_res = input("Grid resolution in km: ") 

click = Click(ax, grid_gen)
#fig.canvas.mpl_connect("button_press_event", onclick)
plt.show()
