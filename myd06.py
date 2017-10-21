"""
Copyright (C) 2014 The HDF Group
Copyright (C) 2014 John Evans

This example code illustrates how to access and visualize a LAADS MOD06 swath
file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python MOD06_L2_Cloud_Optical_Thickness.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

The netcdf library must be compiled with HDF4 support in order for this example
code to work.  Please see the README for details.
"""
from __future__ import print_function, division
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

plt.rcParams["figure.figsize"] = [8,8]

USE_NETCDF4 = False

def run(FILE_NAME):
    
    GEO_FILE_NAME = 'MYD03.A2017270.1805.006.2017271151004.hdf'
    #GEO_FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], GEO_FILE_NAME)
    DATAFIELD_NAME = 'Cloud_Effective_Radius'
   
    from pyhdf.SD import SD, SDC
    hdf = SD(FILE_NAME, SDC.READ)

    # Read dataset.
    data2D = hdf.select(DATAFIELD_NAME)
    data = data2D[:,:].astype(np.double)
    print(data.shape,data)
    hdf_geo = SD(GEO_FILE_NAME, SDC.READ)
    print(hdf_geo)
    # Read geolocation dataset from MOD03 product.
    lat = hdf_geo.select('Latitude')
    latitude = lat[:,:]
    lon = hdf_geo.select('Longitude')
    longitude = lon[:,:]
    print(latitude.shape,latitude)
    print(longitude.shape,longitude)
    
    # Retrieve attributes.
    attrs = data2D.attributes(full=1)
    lna=attrs["long_name"]
    long_name = lna[0]
    aoa=attrs["add_offset"]
    add_offset = aoa[0]
    fva=attrs["_FillValue"]
    _FillValue = fva[0]
    sfa=attrs["scale_factor"]
    scale_factor = sfa[0]        
    vra=attrs["valid_range"]
    valid_min = vra[0][0]        
    valid_max = vra[0][1]        
    ua=attrs["units"]
    units = ua[0]

    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor 
    data = np.ma.masked_array(data, np.isnan(data))
    
    fig = plt.figure()      
    # Render the plot in a south plar stereographic projection.
    m = Basemap(width=800000,height=800000,
            resolution='l',projection='stere',\
            lat_ts=45,lat_0=45,lon_0=-76.)
   
    m.drawcoastlines(linewidth=1)
    m.drawparallels(np.arange(40., 50., 2.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-80, -70., 2.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, rasterized=True, cmap=plt.cm.nipy_spectral)
    lon = -75.80307
    lat = 45.37165
    x, y = m(lon,lat)
    print(x,y)
    plt.annotate('B', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color='white')
    cb = m.colorbar()
    cb.set_label(units)
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n'.format(basename), fontsize=12)
    #fig = plt.gcf()
    plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)   

    GEO_FILE_NAME = 'MYD06_L2.A2017270.1805.006.2017271155804.hdf'
    DATAFIELD_NAME = 'Cloud_Top_Temperature'
   
    from pyhdf.SD import SD, SDC
    hdf = SD(FILE_NAME, SDC.READ)

    # Read dataset.
    data2D = hdf.select(DATAFIELD_NAME)
    data = data2D[:,:].astype(np.double)
    print(data.shape,data)
    hdf_geo = SD(GEO_FILE_NAME, SDC.READ)
    print(hdf_geo)
    # Read geolocation dataset from MOD03 product.
    lat = hdf_geo.select('Latitude')
    latitude = lat[:,:]
    lon = hdf_geo.select('Longitude')
    longitude = lon[:,:]
    print(latitude.shape,latitude)
    print(longitude.shape,longitude)
    
    # Retrieve attributes.
    attrs = data2D.attributes(full=1)
    lna=attrs["long_name"]
    long_name = lna[0]
    aoa=attrs["add_offset"]
    add_offset = aoa[0]
    fva=attrs["_FillValue"]
    _FillValue = fva[0]
    sfa=attrs["scale_factor"]
    scale_factor = sfa[0]        
    vra=attrs["valid_range"]
    valid_min = vra[0][0]        
    valid_max = vra[0][1]        
    ua=attrs["units"]
    units = ua[0]

    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor 
    data = np.ma.masked_array(data, np.isnan(data))
    
    fig = plt.figure()      
    # Render the plot in a south plar stereographic projection.
    m = Basemap(width=800000,height=800000,
            resolution='l',projection='stere',\
            lat_ts=45,lat_0=45,lon_0=-76.)
   
    m.drawcoastlines(linewidth=1)
    m.drawparallels(np.arange(40., 50., 2.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-80, -70., 2.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, rasterized=True, cmap=plt.cm.nipy_spectral)
    lon = -75.80307
    lat = 45.37165
    x, y = m(lon,lat)
    print(x,y)
    plt.annotate('B', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color='white')
    cb = m.colorbar()
    cb.set_label(units)
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n'.format(basename), fontsize=12)
    #fig = plt.gcf()
    plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)   
    
    GEO_FILE_NAME = 'MYD06_L2.A2017270.1805.006.2017271155804.hdf'
    DATAFIELD_NAME = 'Brightness_Temperature'
   
    from pyhdf.SD import SD, SDC
    hdf = SD(FILE_NAME, SDC.READ)

    # Read dataset.
    data2D = hdf.select(DATAFIELD_NAME)
    data = data2D[:,:].astype(np.double)
    print(data.shape,data)
    hdf_geo = SD(GEO_FILE_NAME, SDC.READ)
    print(hdf_geo)
    # Read geolocation dataset from MOD03 product.
    lat = hdf_geo.select('Latitude')
    latitude = lat[:,:]
    lon = hdf_geo.select('Longitude')
    longitude = lon[:,:]
    print(latitude.shape,latitude)
    print(longitude.shape,longitude)
    
    # Retrieve attributes.
    attrs = data2D.attributes(full=1)
    lna=attrs["long_name"]
    long_name = lna[0]
    aoa=attrs["add_offset"]
    add_offset = aoa[0]
    fva=attrs["_FillValue"]
    _FillValue = fva[0]
    sfa=attrs["scale_factor"]
    scale_factor = sfa[0]        
    vra=attrs["valid_range"]
    valid_min = vra[0][0]        
    valid_max = vra[0][1]        
    ua=attrs["units"]
    units = ua[0]

    invalid = np.logical_or(data > valid_max,
                            data < valid_min)
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = (data - add_offset) * scale_factor 
    data = np.ma.masked_array(data, np.isnan(data))
    print(data.shape,data)
    data = data[0,:,:]
    
    fig = plt.figure()      
    # Render the plot in a south plar stereographic projection.
    m = Basemap(width=800000,height=800000,
            resolution='l',projection='stere',\
            lat_ts=45,lat_0=45,lon_0=-76.)
   
    m.drawcoastlines(linewidth=1)
    m.drawparallels(np.arange(40., 50., 2.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-80, -70., 2.), labels=[0, 0, 0, 1])
    m.pcolormesh(longitude, latitude, data, latlon=True, rasterized=True, cmap=plt.cm.nipy_spectral)
    lon = -75.80307
    lat = 45.37165
    x, y = m(lon,lat)
    print(x,y)
    plt.annotate('B', xy=(x, y), xycoords='data', xytext=(x, y), textcoords='data', color='white')
    cb = m.colorbar()
    cb.set_label(units)
    basename = os.path.basename(FILE_NAME)
    plt.title('{0}\n'.format(basename), fontsize=12)
    #fig = plt.gcf()
    plt.show()
    pngfile = "{0}.py.png".format(basename)
    fig.savefig(pngfile)   

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'MYD06_L2.A2017270.1805.006.2017271155804.hdf'
    """
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass
    """
    run(hdffile)
    
