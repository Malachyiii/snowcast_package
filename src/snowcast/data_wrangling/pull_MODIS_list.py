# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 15:51:35 2022

@author: mmoran
"""
import ee
import signal
import shapely
import platform
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
      ee.Initialize()
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
    
        #We need a start date and an end date. Just to ensure we pull something
        start_date = date - timedelta(days = 7)
        end_date = date + timedelta(days = 1)
    
        #First I get the image collection from the MODIS data, filter it only to the days in question
        #and select my bands, then sort so the most recent day in the group is at the top
        DatedCollection = Collection.filter(ee.Filter.date(start_date, end_date)) \
                                    .filter(ee.Filter.notNull(['system:index'])) \
                                    .sort('system:index', False)
    
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