import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import ESMF
import os

def write_nc(lon_array,lat_array,var):
                #fid = nc.Dataset('/Users/elizabethdrenkard/external_data/ERAinterim/land_sea_mask/lsm_LENS_roms.nc','w', format='NETCDF3_CLASSIC')
                fid = nc.Dataset('bathym2.nc','w', format='NETCDF3_CLASSIC') 
                # dimensions
                fid.createDimension('nyp', lat_array.shape[0])
                fid.createDimension('nxp', lon_array.shape[1])
                # variables
                latitudes  = fid.createVariable('y', 'f8', ('nyp','nxp',))
                longitudes = fid.createVariable('x', 'f8', ('nyp','nxp',))
                variable   = fid.createVariable('h', 'f8', ('nyp','nxp',))
                # data
                #print var.shape, lat_array.shape, lon_array.shape
                latitudes[:]  = lat_array
                longitudes[:] = lon_array
                variable[:] = var 

                # attributes
                longitudes.units = "degrees_east"
                longitudes.valid_min = lon_array.min()
                longitudes.valid_max = lon_array.max()
                longitudes.long_name = "longitude"

                latitudes.units = "degrees_north"
                latitudes.valid_min = lat_array.min()
                latitudes.valid_max = lat_array.max()
                latitudes.long_name = "latitude"


                variable.long_name = 'bathymetry'
                variable.units = 'meters'
                variable.coordinates = "x y"
                variable.valid_range = var.min() , var.max()

                # close
                fid.close()
                return None

#Define source field information
topo_dir = '/work/Liz.Drenkard/bathy_interp/etopo/'

topo = np.loadtxt(os.path.join(topo_dir,'etopo20data'))
lons = np.loadtxt(os.path.join(topo_dir,'etopo20lons'))
lats = np.loadtxt(os.path.join(topo_dir,'etopo20lats'))

# Bilinear interpolation of source lat/lon
# NOTE: I think this is fine for longitudes but not sure latitude in degrees.
s_corn_lat = (lats[:-1] + lats[1:])/2
s_corn_lon = (lons[:-1] + lons[1:])/2
corn_topo = topo[1:-1,1:-1]

sourcegrid = ESMF.Grid(np.array(corn_topo.shape), staggerloc = ESMF.StaggerLoc.CORNER, \
                       coord_sys = ESMF.CoordSys.SPH_DEG)

source_lon = sourcegrid.get_coords(0, staggerloc=ESMF.StaggerLoc.CORNER)
source_lat = sourcegrid.get_coords(1, staggerloc=ESMF.StaggerLoc.CORNER)

X1, Y1 = np.meshgrid(s_corn_lon,s_corn_lat)

source_lon[...] = X1
source_lat[...] = Y1

sourcefield = ESMF.Field(sourcegrid, name = 'ETOPO_bathy')
srcfracfield = ESMF.Field(sourcegrid, 'srcfracfield')

sourcefield.data[...] = corn_topo


# DEFINE DESTINATION GRID

# define corner corrdinates via bilinear interpolation of hgrid points 
fid = nc.Dataset('/work/eay/grid_gen/ocean_hgrid_final2.nc')
d_lat = fid.variables['y'][:]
d_lon = fid.variables['x'][:]


d_corn_lat = (d_lat[:-1,:-1] + d_lat[:-1,1:] + d_lat[1:,1:] + d_lat[1:,:-1])/4
d_corn_lon = (d_lon[:-1,:-1] + d_lon[:-1,1:] + d_lon[1:,1:] + d_lon[1:,:-1])/4

destgrid = ESMF.Grid(np.array((d_corn_lat.shape[0]-1,d_corn_lat.shape[1]-1)), \
           staggerloc = ESMF.StaggerLoc.CORNER, \
           coord_sys = ESMF.CoordSys.SPH_DEG)

dest_lon = destgrid.get_coords(0,staggerloc=ESMF.StaggerLoc.CORNER)
dest_lat = destgrid.get_coords(1,staggerloc=ESMF.StaggerLoc.CORNER)

dest_lon[...] = d_corn_lon
dest_lat[...] = d_corn_lat

destfield = ESMF.Field(destgrid, name = 'MOM6_regional')

# DEFINE INTERPOLATION FUNCTION
regrid = ESMF.Regrid(sourcefield, destfield, regrid_method = ESMF.RegridMethod.CONSERVE, 
         src_mask_values=np.array([0], dtype=np.int32), src_frac_field=srcfracfield,
         norm_type=ESMF.NormType.FRACAREA, unmapped_action = ESMF.UnmappedAction.IGNORE)

destfield = regrid(sourcefield, destfield)

write_nc(d_lon[1:-1,1:-1],d_lat[1:-1,1:-1],destfield.data)

plt.pcolormesh(d_corn_lon, d_corn_lat,destfield.data)

plt.show()

# SAVE BATHYMETRY FIELD AS NETCDF FILE

