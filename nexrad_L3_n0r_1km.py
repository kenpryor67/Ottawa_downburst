# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 11:33:01 2017

@author: Ken.Pryor
"""

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm

plt.rcParams["figure.figsize"] = [8,8]

def read_Z_plot(ncf):
    nc_fid = Dataset(ncf, 'r')
    Z = nc_fid.variables["Reflectivity"][:]  # shape lat, lon as shown above
    X = nc_fid.variables["x"][:]
    Y = nc_fid.variables["y"][:]
    lats = nc_fid.variables['lat'][:]  # extract/copy the data
    lons = nc_fid.variables['lon'][:]
    nc_fid.close()
    return Z, X, Y, lats, lons

Z_file = 'Level3_Composite_n0r_1km_20170927_1805.gini.nc4'
Z, X, Y, lats, lons = read_Z_plot(Z_file)
Z= Z[0,:,:]
print('Z shape', Z.shape, Z)
print('X shape, Y shape', X.shape, Y.shape)
print('lats shape, lons shape', lats.shape, lons.shape)
makepos = np.full((801,775),-1)
Z_pos = np.multiply(Z,makepos)
Z_pos = Z_pos/10000
print('Z_pos shape', Z_pos.shape, Z_pos)
np.delete(Z,[0,1,2])
np.delete(X,[0,1,2])
np.delete(Y,[0,1,2])

# CREATE A MAP

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(width=600000,height=600000,
            resolution='l',projection='stere',\
            lat_ts=45,lat_0=45,lon_0=-76.)
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,Z_pos,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(43.,48.,2.),labels=[True,False,False,False],color="white")
m.drawmeridians(np.arange(-78.,-73.,2.),labels=[False,False,False,True],color="white")
m.drawcountries(color="white")
lon = -75.80307
lat = 45.37165
x, y = m(lon,lat)
print(x,y)
np.delete(Z_pos,[0,1,2])
np.delete(lats,[0,1])
np.delete(lons,[0,1])
plt.annotate('B', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color="white")
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('Reflectivity (Z) 1805 UTC 27 September 2017')
plt.savefig("Z_0927_1805.png",bbox_inches='tight')
plt.show()

Z_file = 'Level3_Composite_n0r_1km_20170927_1905.gini.nc4'
Z, X, Y, lats, lons = read_Z_plot(Z_file)
Z= Z[0,:,:]
print('Z shape', Z.shape, Z)
print('X shape, Y shape', X.shape, Y.shape)
print('lats shape, lons shape', lats.shape, lons.shape)
makepos = np.full((801,775),-1)
Z_pos = np.multiply(Z,makepos)
Z_pos = Z_pos/10000
print('Z_pos shape', Z_pos.shape, Z_pos)
np.delete(Z,[0,1,2])
np.delete(X,[0,1,2])
np.delete(Y,[0,1,2])

# CREATE A MAP

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(width=600000,height=600000,
            resolution='l',projection='stere',\
            lat_ts=45,lat_0=45,lon_0=-76.)
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,Z_pos,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(43.,48.,2.),labels=[True,False,False,False],color="white")
m.drawmeridians(np.arange(-78.,-73.,2.),labels=[False,False,False,True],color="white")
m.drawcountries(color="white")
lon = -75.80307
lat = 45.37165
x, y = m(lon,lat)
print(x,y)
np.delete(Z_pos,[0,1,2])
np.delete(lats,[0,1])
np.delete(lons,[0,1])
plt.annotate('B', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color="white")
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('Reflectivity (Z) 1905 UTC 27 September 2017')
plt.savefig("Z_0927_1905.png",bbox_inches='tight')
plt.show()

Z_file = 'Level3_Composite_n0r_1km_20170927_1910.gini.nc4'
Z, X, Y, lats, lons = read_Z_plot(Z_file)
Z= Z[0,:,:]
print('Z shape', Z.shape, Z)
print('X shape, Y shape', X.shape, Y.shape)
print('lats shape, lons shape', lats.shape, lons.shape)
makepos = np.full((801,775),-1)
Z_pos = np.multiply(Z,makepos)
Z_pos = Z_pos/10000
print('Z_pos shape', Z_pos.shape, Z_pos)
np.delete(Z,[0,1,2])
np.delete(X,[0,1,2])
np.delete(Y,[0,1,2])

# CREATE A MAP

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(width=600000,height=600000,
            resolution='l',projection='stere',\
            lat_ts=45,lat_0=45,lon_0=-76.)
m.drawmapboundary(fill_color='0.0')
im1 = m.pcolormesh(lons,lats,Z_pos,shading='gouraud',cmap=plt.cm.nipy_spectral,latlon=True) # shading can also be "flat"
m.drawparallels(np.arange(43.,48.,2.),labels=[True,False,False,False],color="white")
m.drawmeridians(np.arange(-78.,-73.,2.),labels=[False,False,False,True],color="white")
m.drawcountries(color="white")
lon = -75.80307
lat = 45.37165
x, y = m(lon,lat)
print(x,y)
np.delete(Z_pos,[0,1,2])
np.delete(lats,[0,1])
np.delete(lons,[0,1])
plt.annotate('B', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color="white")
cb = m.colorbar(im1,"right", size="5%", pad="10%")
ax.set_title('Reflectivity (Z) 1910 UTC 27 September 2017')
plt.savefig("Z_0927_1910.png",bbox_inches='tight')
plt.show()