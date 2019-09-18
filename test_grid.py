import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


fid = nc.Dataset('ocean_hgrid.nc')
x = fid.variables['x'][:]
y = fid.variables['y'][:]

plt.figure()
m = Basemap()

m.pcolor(x,y,y)
m.drawcoastlines()
plt.colorbar()
plt.show()
