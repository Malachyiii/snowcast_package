# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 15:51:35 2022

@author: mmoran
"""
import ee
import signal
import shapely
import platform
import rasterio
import rasterio.features
import rioxarray
import planetary_computer
import xrspatial
import numpy

import urllib
from PIL import Image
from datashader.transfer_functions import shade, stack
from pystac_client import Client
from shapely import wkt
from datetime import datetime, timedelta
from pandas import to_datetime
from numpy import divide


if platform.system() == 'Windows':
    pass

else:
    class TimeoutException(Exception):   # Custom exception class
        pass
    
    def timeout_handler(signum, frame):   # Custom signal handler
        raise TimeoutException
    
    signal.signal(signal.SIGALRM, timeout_handler)
    
##############Pull Modis images###########################################

def pull_MODIS_image(geometry, date, modis, buffer_percent = 0.0):
    
    
    try:
        if not isinstance(date, datetime):
            date = to_datetime(date)
    except Exception as e:
        print(e)
        print("You did not enter a date for the date variable")
        return
    
    try:
        if not isinstance(geometry, shapely.geometry.base.BaseGeometry):
            geometry = wkt.loads(geometry)
    except Exception as e:
        print(e)
        print("You did not enter a geometry for the geometry variable")
        return
    



    #We need a start date and an end date. Just like a regular python slice, 
    #the end date is not included, so by using a 1 day frame, I am actually limiting
    #the range to up to the day in question
    start_date = date - timedelta(days = 7)
    end_date = date + timedelta(days = 1)
    
    aoi = ee.Geometry.Polygon(list(shapely.affinity.scale(geometry,
                                                        xfact=(1.+ buffer_percent), 
                                                        yfact=(1.+ buffer_percent)).exterior.coords))


    #First I get the image collection from the MODIS data, filter it only to the days in question
    #and select my bands, then sort so the most recent day in the group is at the top
    try:
     point = ee.ImageCollection(f'MODIS/006/{modis}') \
                .filterDate(start_date, end_date) \
                .filterBounds(aoi) \
                .filter(ee.Filter.notNull(['system:index'])) \
                .select(['NDSI_Snow_Cover', 'Snow_Albedo_Daily_Tile', 'NDSI']) \
                .sort('system:index', False) \
                .first() \
                .unmask(0)
      
    except Exception as e:
      print(e)
      print("Some Error with Image Collection")
      print("Have you logged in to google earth engine?")
      return 
        
    
    # Get individual band arrays and build them into an RGB image
    SnowCover = point.select('NDSI_Snow_Cover')
    Albedo = point.select('Snow_Albedo_Daily_Tile')
    Raw = point.select('NDSI').divide(100)

    point = point.addBands(SnowCover.rename('SnowCover'))
    point = point.addBands(Albedo.rename('Albedo'))
    point = point.addBands(Raw.rename('Raw'))
  
    rgb = ['SnowCover', 'Albedo', 'Raw']
    
    
    #get URL for image
    url = point.select(rgb).clip(aoi).getThumbURL({'min': 0, 'max':100, 'region': aoi, 'format': 'jpg'})
    
    try:
       r = urllib.request.urlopen(url)
       im = Image.open(r)
       return numpy.array(im)
      
    except Exception as e:
       print(e)



##############PULL MODIS LIST#####################################

def pull_MODIS_list(geometry, date, modis, signal_timer = 5):
    
    '''
    This function pulls the data from the specified MODIS snowcover satellite
    as three columns of data. One for The snow cover, one for the Albedo, and one
    for the unmodified NDSI measurment
    '''
    try:
        if not isinstance(date, datetime):
            date = to_datetime(date)
    except Exception as e:
        print(e)
        print("You did not enter a date for the date variable")
        return
    
    try:
        if not isinstance(geometry, shapely.geometry.base.BaseGeometry):
            geometry = wkt.loads(geometry)
    except Exception as e:
        print(e)
        print("You did not enter a geometry for the geometry variable")
        return
    
    
    try:
      Collection = ee.ImageCollection(f'MODIS/006/{modis}') \
                  .select(['NDSI_Snow_Cover', 'Snow_Albedo_Daily_Tile', 'NDSI'])
      
    except Exception as e:
      print(e)
      print("Some Error with Image Collection")
      print("Have you logged in to google earth engine?")
      return 
  


    still_working = True
    while still_working:
      if platform.system() == 'Windows':
          print("Running on a Windows system, Timeout recovery is not available")
      else:
        signal.alarm(signal_timer)
      
      try:
        row = [geometry, date]
    
        #We need a start date and an end date. Just like a regular python slice, 
        #the end date is not included, so by using a 1 day frame, I am actually limiting
        #the range to only the day in question
        start_date = date - timedelta(days = 7)
        end_date = date + timedelta(days = 1)
    
        #First I get the image collection from the MODIS data, filter it only to the days in question
        #and select my bands, then sort so the most recent day in the group is at the top
        DatedCollection = Collection.filter(ee.Filter.date(start_date, end_date)) \
                                    .filter(ee.Filter.notNull(['system:index'])) \
                                    .sort('system:index', False)
    
        #Because the image collection is limited to a single day, there is only one image
        #So I just take it
        point = DatedCollection.first().unmask(0)
    
        aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
    
        bands = point.reduceRegion(reducer = ee.Reducer.mean(),
        geometry= aoi)
    
        bands = bands.toArray(['NDSI_Snow_Cover', 'Snow_Albedo_Daily_Tile', 'NDSI']).getInfo()
        bands = divide(bands, [100,100,10000] )
    
        row.extend(bands)
      
      except TimeoutException:
        print("Request Timeout, Retrying")
      
      except Exception as e:
        print(e)
        print("Some other Error")
      else: 
        if platform.system() == 'Windows':
          pass
        else:
          signal.alarm(0)
        still_working = False

    return row

######################### COPERNICUS########################
def get_copernicus(geometry, signal_timer =5):
    
    try:
        if not isinstance(geometry, shapely.geometry.base.BaseGeometry):
            geometry = wkt.loads(geometry)
    except Exception as e:
        print(e)
        print("You did not enter a geometry for the geometry variable")
        return
    

    still_working = True
    while still_working:
      try:
        client = Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        ignore_conformance=True,
        )
    
      except Exception as e:
        print(e)
        print("Some Error with Image Collection")
      else: 
        still_working = False  

    error_num = 0
    still_working = True
    while still_working:
      
      if platform.system() == 'Windows':
          print("Running on a Windows system, Timeout recovery is not available")
      else:
        signal.alarm(signal_timer)
      
        
      try:    
    
        # Adapted from https://planetarycomputer.microsoft.com/dataset/cop-dem-glo-90#Example-Notebook :
          if error_num == 0:
              area_of_interest = {
              "type": "Polygon",
              "coordinates": [geometry.exterior.coords],
              }
                
              bbox = rasterio.features.bounds(area_of_interest)
              search = client.search(
                collections=["cop-dem-glo-30"],
                bbox = bbox,
                )
                
              min_lon = geometry.bounds[0]
              min_lat = geometry.bounds[1]
              max_lon = geometry.bounds[2]
              max_lat = geometry.bounds[3]
            
          elif error_num > 1:
              area_of_interest = {
              "type": "Polygon",
              "coordinates": [list(shapely.affinity.scale(geometry,
                                                          xfact=(1.+ error_num), 
                                                          yfact=(1.+ error_num)).exterior.coords)],
              }
              print(f"Boundaries scaled by {(error_num*100)} %")
              bbox = rasterio.features.bounds(area_of_interest)
              search = client.search(
                collections=["cop-dem-glo-30"],
                intersects= area_of_interest,
                )
            
              min_lon = shapely.affinity.scale(geometry,
                                                          xfact=(1.+ error_num), 
                                                          yfact=(1.+ error_num)).bounds[0]
              min_lat = shapely.affinity.scale(geometry,
                                                          xfact=(1.+ error_num), 
                                                          yfact=(1.+ error_num)).bounds[1]
              max_lon = shapely.affinity.scale(geometry,
                                                          xfact=(1.+ error_num), 
                                                          yfact=(1.+ error_num)).bounds[2]
              max_lat = shapely.affinity.scale(geometry,
                                                          xfact=(1.+ error_num), 
                                                          yfact=(1.+ error_num)).bounds[3]
    
          
    
          
          items = list(search.get_items())
    
          signed_asset = planetary_computer.sign(items[0].assets["data"])
          
          data = (
              rioxarray.open_rasterio(signed_asset.href)
              .squeeze()
              .drop("band")
          )
    
    
          mask_lon = (data.x >= min_lon) & (data.x <= max_lon)
          mask_lat = (data.y >= min_lat) & (data.y <= max_lat)
    
    
          cropped_data = data.where(mask_lon & mask_lat, drop=True)
    
          #hillshade = xrspatial.hillshade(cropped_data)
          img = stack(shade(cropped_data, cmap="red"), 
                      shade(xrspatial.slope(cropped_data), cmap="blue"),
                      shade(xrspatial.aspect(cropped_data), cmap="green"))
          
      
        
      except Exception as e:
        print("Some other Error")
        print(e)
        error_num += 1
      
      except TimeoutException:
          print("Request Timeout")
            
      else: 
        if platform.system() == 'Windows':
            pass
        else: 
            signal.alarm(0)
        still_working = False

    return numpy.array(img.to_pil())



##########Sentinel 1 #########################################################

def pull_Sentinel1(geometry, date):
    
    try:
        if not isinstance(date, datetime):
            date = to_datetime(date)
    except Exception as e:
        print(e)
        print("You did not enter a date for the date variable")
        return
    
    try:
        if not isinstance(geometry, shapely.geometry.base.BaseGeometry):
            geometry = wkt.loads(geometry)
    except Exception as e:
        print(e)
        print("You did not enter a geometry for the geometry variable")
        return
    
   #define area of interest by coordinates
    aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
   #print(aoi)
    start_date = date - timedelta(days = 80)
    end_date = date + timedelta(days = 1)

    try:
      #print("calculating")
      # Sentinel-1 image filtered on date range and on aoi
        se2 = ee.ImageCollection('COPERNICUS/S1_GRD')\
          .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))\
          .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))\
          .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'HV').Not())\
          .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'HH').Not())\
          .filterDate(start_date, end_date)\
          .filterBounds(aoi)\
          .sort('system:time_start', False)\
          .first()

        #Create a ratio band of VV/VH
        VVVH = (se2.select('VV').divide(se2.select('VH')))

        se2 = se2.addBands(VVVH.rename('VVVH'))

        rgb = ['VV', 'VH', 'VVVH']


        url = se2.select(rgb).clip(aoi).getThumbURL({'min': -50, 'max': 1, 'region': aoi, 'format': 'jpg'})

        #now I open the url and download the image to the specified file location
        #response = requests.get(url, stream=True).content        
        
        try:
            r = urllib.request.urlopen(url)
            im = Image.open(r)
            return numpy.array(im)
        except Exception as e:
            print(e)

    except Exception as e:
        print(e)
        print("Have you logged in to google earth engine?")
        return


#############Sentinel 2a ########################################################
def pull_Sentinel2a(geometry, date):

    try:
        if not isinstance(date, datetime):
            date = to_datetime(date)
    except Exception as e:
        print(e)
        print("You did not enter a date for the date variable")
        return
    
    try:
        if not isinstance(geometry, shapely.geometry.base.BaseGeometry):
            geometry = wkt.loads(geometry)
    except Exception as e:
        print(e)
        print("You did not enter a geometry for the geometry variable")
        return
    

    aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
    #print(aoi)
    start_date = date - timedelta(days = 80)
    end_date = date + timedelta(days = 1)
    
    try:
        #print("calculating")
        # Sentinel-2 image filtered on date range and on aoi
        se2 = ee.ImageCollection('COPERNICUS/S2')\
            .filterDate(start_date, end_date)\
            .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 10))\
            .filterBounds(aoi) \
            .select(['B8A', 'B11', 'B12'])\
            .sort('system:time_start', False)\
            .first()
    
          #Create a ratio band of (B8a-B11)/(B8a+B11)
    
        BRatio = se2.expression(
              '((B8 - B11)/(B8 + B11))*100', {
              'B8': se2.select('B8A'),
              'B11': se2.select('B11')
              })
    
        se2 = se2.addBands(BRatio.rename('BRatio'))
    
    
        rgb = ['BRatio', 'B11', 'B12']
    
    
        url = se2.select(rgb).clip(aoi).getThumbURL({'min': -500, 'max':500, 'region': aoi, 'format': 'jpg'})
    
    
    
        try:
            r = urllib.request.urlopen(url)
            im = Image.open(r)
            return numpy.array(im)
        except Exception as e:
            print(e)
    
    except Exception as e:
        print(e)
        print("Have you logged in to google earth engine?")
        return
    
    
##################Sentinel 2b################################

def pull_Sentinel2b(geometry, date):
    
    try:
        if not isinstance(date, datetime):
            date = to_datetime(date)
    except Exception as e:
        print(e)
        print("You did not enter a date for the date variable")
        return
    
    try:
        if not isinstance(geometry, shapely.geometry.base.BaseGeometry):
            geometry = wkt.loads(geometry)
    except Exception as e:
        print(e)
        print("You did not enter a geometry for the geometry variable")
        return

    #define area of interest by coordinates
    aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
    #print(aoi)
    start_date = date - timedelta(days = 80)
    end_date = date + timedelta(days = 1)
    
    try:
        #print("calculating")
        # Sentinel-2 image filtered on date range and on aoi
        se2 = ee.ImageCollection('COPERNICUS/S2')\
            .filterDate(start_date, end_date)\
            .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 10))\
            .filterBounds(aoi) \
            .select(['B8', 'B4', 'B2'])\
            .sort('system:time_start', False)\
            .first()
    
          #Create a ratio band of (B8a-B11)/(B8a+B11)
    
        BRatio = se2.expression(
              '((B8 - B4)/(B8 + B4))*10000', {
              'B8': se2.select('B8'),
              'B4': se2.select('B4')
              })
    
        se2 = se2.addBands(BRatio.rename('BRatio'))
    
    
        rgb = ['BRatio', 'B2', 'B4']
    
    
        url = se2.select(rgb).clip(aoi).getThumbURL({'min': -100, 'max':5000, 'region': aoi, 'format': 'jpg'})
    
    
    
        try:
            r = urllib.request.urlopen(url)
            im = Image.open(r)
            return numpy.array(im)
        
        except Exception as e:
            print(e)
    
    except Exception as e:
        print(e)
        print("Have you logged in to google earth engine?")
        return
    

############# testing###########
#from shapely.geometry import Polygon

#polygon = Polygon([[-119.49, 38.22], [-119.59, 38.22], [-119.59, 38.32], [-119.49, 38.32], [-119.49, 38.22]])

#print(pull_MODIS_list(polygon, '2018-10-10', 'MOD10A1'))

#print(get_copernicus(polygon))

#print(pull_Sentinel1(polygon, '2018-10-10'))

#print(pull_Sentinel2a(polygon, '2018-10-10'))

#print(pull_Sentinel2b(polygon, '2018-10-10'))

#print(pull_MODIS_image(polygon, '2018-12-12', 'MOD10A1', buffer_percent=0.0))