
from paths import *

from osgeo import gdal
from osgeo import osr
import os
import rasterio
from rasterio.mask import mask

import numpy as np
import pandas as pd
import rasterio
import geopandas
import warnings
import shapely


earth_radius = 6377.6
deg = earth_radius * np.pi / 180  # 111.3 km
arcmin = deg / 60  # 1.8 km
arcsec = deg / 3600.  # 31 m
cell_size = 9.  # km

def _names_match(name, mtn_name):
    ref = name.lower()
    query = mtn_name.lower()
    if ' ' in ref:
        if (ref.split(' ')[-1][:4] not in 
            ['moun', 'mts.', 'plat', 'isle', 'isla', 'peni', 'shie'] or
            ref.split(' ')[0] == 'new'):
            ref, query = ref.split(' ')[-1], query.split(' ')[-1]
    return ref[:4] == query[:4]

def _bin_map(matrix, factor=5, value='relief'):
    mat = np.where(matrix > 0., matrix, np.nan)
    height, width = mat.shape
    res_h, res_w = int(height / factor), int(width / factor)
    tail_h, tail_w = height % factor, width % factor
    in_mat = mat[:height - tail_h, :width - tail_w].reshape(res_h, factor, res_w, factor)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if value == 'mean':
            return np.nanmean(np.nanmean(in_mat, axis=1), axis=2)
        elif value == 'std':
            return np.sqrt(np.nanmean(np.nanmean(in_mat ** 2, axis=1), axis=2) -
                           np.nanmean(np.nanmean(in_mat, axis=1), axis=2) ** 2)
        elif value == 'min':
            return np.nanmin(np.nanmin(in_mat, axis=1), axis=2)
        elif value == 'max':
            return np.nanmax(np.nanmax(in_mat, axis=1), axis=2)
        else:  # value == 'relief'
            return (np.nanmean(np.nanmean(in_mat, axis=1), axis=2) -
                     np.nanmin(np.nanmin(in_mat, axis=1), axis=2))


def _get_map(mt_name, topo, moho):
    if type(mt_name) == int:
        name = shp.name[mt_name]
        geom = shp.geometry[mt_name]
    else:
        for name, geom in zip(shp.name, shp.geometry):
            if _names_match(name, mt_name):
                break
    
    geom_map = [shapely.geometry.mapping(geom)]
    maps_topo = mask(topo, geom_map, crop=True)[0][0]

    maps_moho, _ = mask(moho, geom_map, crop=True)

    return name, np.where((maps_topo > 0), maps_topo, np.nan), geom, maps_moho

def _solve(p, A):
    Delta = -16*A + p**2
    if Delta < 0:
        return np.nan, np.nan
    else:
        return p/4 + np.sqrt(Delta)/4, p/4 - np.sqrt(Delta)/4

def _cal_eros(maps):  # From Vance et al. (2003)
    factor = int(cell_size / arcmin) + 1
    rlf = _bin_map(maps, factor=factor, value='relief')
    eros_map = 1e-2 * 10 ** (rlf / 500.)
    return np.nanmean(eros_map), eros_map


def _loadxyz2grd(fname,order,dlon,dlat):
    """
    load .xyz file to grid
    input:
    fname: load file name
    order: the squence of lons and lats saving in .xyz 
    dlon:
    dlat: 
    
    data: elevation data, an array in size of (n_lat,n_lon) 
    latRange: range of latitude, an array as [minlat,maxlat]
    lonRange: range of longitude, an array as [minlon,maxlon]
    dtype: dtype in gdal, as gdal.GDT_Byte or gdal.GDT_Float32
    """ 
    data = np.loadtxt(fname)
    if order == "lonlat":
        lons,lats,z = data[:,0],data[:,1],data[:,2]
    if order == "latlon":
        lats,lons,z = data[:,0],data[:,1],data[:,2]
    minlon,maxlon,minlat,maxlat = lons.min(),lons.max(),lats.min(),lats.max()
    extent = (minlon,maxlon,minlat,maxlat)
    nlon = int((maxlon-minlon)/dlon+1)
    nlat = int((maxlat-minlat)/dlat+1)
    size = nlon*nlat
    data = np.zeros((nlat,nlon))
    for i in range(size):
        lon_idx = np.around((lons[i]-minlon)/dlon).astype("int")
        lat_idx = np.around((lats[i]-maxlat)/dlat).astype("int")
        data[lat_idx,lon_idx]=z[i]
    return data,extent

def _array2geotiff(fname, data, extent, dtype):   
    """
    save GeoTiff file from the array of dem data
    input:
    fname: save file name
    data: elevation data, an array in size of (n_lat,n_lon) 
    latRange: range of latitude, an array as [minlat,maxlat]
    lonRange: range of longitude, an array as [minlon,maxlon]
    dtype: dtype in gdal, as gdal.GDT_Byte or gdal.GDT_Float32
    """   
    nx = data.shape[1]
    ny = data.shape[0]
    xmin,xmax,ymin,ymax = extent
    dx = (xmax - xmin) / float(nx)
    dy = (ymax - ymin) / float(ny)
    geotransform = (xmin, dx, 0, ymax, 0, -dy)
    dst = gdal.GetDriverByName('GTiff').Create(fname, nx, ny, 1, dtype)
    dst.SetGeoTransform(geotransform) 
    dst.GetRasterBand(1).WriteArray(data)
    srs = osr.SpatialReference()
    #srs.ImportFromEPSG(4326)
    srs.ImportFromEPSG(3857) # WGS84 lat/long
    #srs.ImportFromEPSG(9122) # WGS84 lat/long
    dst.SetProjection(srs.ExportToWkt())  
    dst.FlushCache() 


def glob_info(mt_name):
    """
    Parameters:
        mt_name: str
            The name of the mountain range. Naive fuzzy search supported
            ("Zagros", "zagros", "Zagr").

    Returns:
        A dict. Its keywords contain
            'name' (the name in the database),
            'length' (length in km),
            'width': (width in km),
            'height': (average height relative to sea level in km),
            'erosion': (average erosion rate in km / Myr).
            'moho': (average moho relative to sea level in km)
    """
    name, maps_topo, geom, maps_moho = _get_map(mt_name, topo, moho)
    L, w = _solve(p=geom.length*deg, A=geom.area*deg**2)

    h_ave = np.nanmean(maps_topo) / 1e3

    eros, eros_map = _cal_eros(maps_topo)

    moho_real = np.extract(maps_moho[0,:,:] != 0., maps_moho[0,:,:])
    moho_ave = moho_real.mean()

    return {'name': name, 'length': L, 'width': w, 'height': h_ave, 'erosion': eros, 'moho': moho_ave}


shp = geopandas.read_file('./ne_50m_geography_regions_polys/ne_50m_geography_regions_polys.shp')
topo = rasterio.open('./ETOPO1_Bed_g_geotiff.tif')

fname = './crust1.0/depthtomoho.xyz'
dlon_crust1pt0_moho, dlat_crust1pt0_moho = 1, 1
data_crust1pt0_moho, extent_crust1pt0_moho = _loadxyz2grd(fname,"lonlat", dlon_crust1pt0_moho, dlat_crust1pt0_moho)

outputfile = './outpu_figs/'
outputPath = os.path.join(os.path.abspath("."),outputfile) 
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

fname_tif = outputfile + 'crust1pt0_moho.tif'
data = np.flipud(data_crust1pt0_moho)
extent = extent_crust1pt0_moho
_array2geotiff(fname_tif, data, extent, gdal.GDT_Float32)
moho = rasterio.open(fname_tif)

print(glob_info("Tibet"))
# print(glob_info("Himalayas")) # The length of the Himalaya-Tibetan plateau is adopted using the length of ‘HIMALAYAS’. 
print(glob_info("Alps"))
print(glob_info("Zagros"))

