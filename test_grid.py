import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

fid = nc.Dataset('ocean_hgrid.nc')
x = fid.variables['x'][:]
y = fid.variables['y'][:]

plt.figure()
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l')
xnew, ynew =m(x,y)
m.pcolor(xnew,ynew,y)
m.drawmapboundary(fill_color='azure')
m.fillcontinents(color='palegoldenrod',lake_color='azure')
m.drawcoastlines()
m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1])

plt.colorbar()
plt.show()
