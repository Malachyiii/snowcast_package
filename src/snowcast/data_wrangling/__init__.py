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
import geopy
import numpy as np
from PIL import Image
from datashader.transfer_functions import shade, stack
from pystac_client import Client
from shapely import wkt
from datetime import datetime, timedelta
from pandas import to_datetime, date_range, DataFrame
from geopandas import GeoDataFrame
from geotiff import GeoTiff
import geopy.distance as distance
from numpy import divide
from shapely.geometry import Polygon
import seaborn as sns
import matplotlib.pylab as plt


if platform.system() == 'Windows':
    pass

else:
    class TimeoutException(Exception):   # Custom exception class
        pass
    
    def timeout_handler(signum, frame):   # Custom signal handler
        raise TimeoutException
    
    signal.signal(signal.SIGALRM, timeout_handler)
    
############## chop_aso to divide ASO up into gridsquares#######################################

def chop_aso(tiff_image_name, groundtruth =  False):
    '''
    

    Parameters
    ----------
    tiff_image_name : str
        Path to the ASO Geotiff image.
    groundtruth : Boolean, optional
        True/False Do you want actual observed values. The default is False.

    Returns
    -------
    df : pd.DataFrame
        Chop_aso is a function that takes in a tiff image from the Airborne Snow
        Observatory in the 50m resolution and cuts it up into a dataframe of 
        1kmx1km chunks. The output is a df with 3 columns. An arbitrary cell_id
        number, the geometry of the square as a shapely polygon, and Snow Water
        Equivalent (SWE). If groundtruth is True the SWE is the mean observed SWE
        for that cell. Otherwise it is a 0 (to be filled in later)
        
    
    '''
    df = []
    
    #Reading the image and transforming it into an array
    geo_tiff = GeoTiff(tiff_image_name)

    zarr_array = geo_tiff.read()

    im_array = np.array(zarr_array)
    
    #We pad the array so it can be divided into equal chunks
    pad_after0 = (20 - (im_array.shape[0] % 20))
    pad_after1 = (20 - (im_array.shape[1] % 20))
    im_array = np.pad(im_array, 
           pad_width = ((0,pad_after0), (0,pad_after1)), 
           constant_values = -9999)
    
    y = 0
    
    #Because we know our resolution we can iterate across the array by pixel
    for i in range(0, im_array.shape[1], 20):
        for j in range(0, im_array.shape[0], 20):
            
            d = distance.distance(kilometers=1)

            #We grab the coordinates of the pixel and the top left corner and 
            #draw a square
            #Going clockwise, from upper-left to lower-left, lower-right...
            coords = geo_tiff.get_wgs_84_coords(i,j)
            p1 = geopy.Point((coords[1], coords[0]))
            p2 = d.destination(point=p1, bearing=180)
            p3 = d.destination(point=p2, bearing=90)
            p4 = d.destination(point=p3, bearing=0)

            #Now we create our polygon and clip it out of an image, we fill
            #in the mask with np.nan
            points = [(p.longitude, p.latitude) for p in [p1,p2,p3,p4]]
            polygon = Polygon(points)

            area_box = [(p1.longitude, p1.latitude), (p3.longitude, p3.latitude)]
            array = geo_tiff.read_box(area_box)
            array[array == -9999] = np.nan
            
            #Now we fill in the SWE and Cell_id based on whether we are looking
            #the ground truth or not
            try:
                if groundtruth:
                    mean = np.nanmean(array)*39.3701 #metric to inches conversion
                else:
                    mean = np.nanmean(array)*0 #returns np.nan if emptyslice and 0 otherwise
                    
            except Exception as e:
                print(e)
                print("Filling in with nan")
                mean = np.nan
                
            if groundtruth:
                cell_id = tiff_image_name[0:-4]+str(y)
            else:
                cell_id = y
            
            row= [cell_id, mean, polygon]
            df.append(row)
            y+=1
        
    df = GeoDataFrame(DataFrame(data = df, columns = ['cell_id', 'SWE', 'geometry']))
    
    return df
    
    
    
##############Pull Modis images###########################################

def pull_MODIS_image(geometry, date, modis, buffer_percent = 0.0):
    '''
    
    
    Parameters
    ----------
    geometry : shapely.geometry, string
        Geometry of the area for which we would like to pull an image.
    date : date, string
        Date for which we would like to pull imagery. Can be a date or string
        defining the date.
    modis : string
        Name of the Modis satellite for which we would like to pull.
    buffer_percent : float, optional
        Float definining how much we would like to buffer around the given
        geometry by 1 = 100% buffer. The default is 0.0.

    Returns
    -------
    np.array
        Returns a 3D Numpy Array of all three bands scaled to be between 0 and 100.

    '''
    
    #Checking for correct input types and converting if otherwise
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
    #the end date is not included, so by using a 1 day delta, I am actually limiting
    #the range to up to the day in question and up to 7 days before
    start_date = date - timedelta(days = 7)
    end_date = date + timedelta(days = 1)
    
    #Creating our polygon with buffer if requested
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
    
    
    Parameters
    ----------
    geometry : shapely.geometry, string
        Geometry of the area for which we would like to pull an image.
    date : date, string
        Date for which we would like to pull imagery. Can be a date or string
        defining the date.
    modis : string
        Name of the Modis satellite for which we would like to pull.
    signal_timer : integer, optional
        Number of seconds to wait before timeout. The default is 5.

    Returns
    -------
    row : list
        This function pulls the data from the specified MODIS snowcover satellite
        as three columns of data. One for The snow cover, one for the Albedo, and one
        for the unmodified NDSI measurment, all normalized to between 0 and 1
        
    '''
    
    #checking input variable types
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
    
    #Connecting to the google earth engine image collection
    try:
      Collection = ee.ImageCollection(f'MODIS/006/{modis}') \
                  .select(['NDSI_Snow_Cover', 'Snow_Albedo_Daily_Tile', 'NDSI'])
      
    except Exception as e:
      print(e)
      print("Some Error with Image Collection")
      print("Have you logged in to google earth engine?")
      return 
  

    #Now we try this until we have a successful connection or timeout. 
    #Not available on windows
    still_working = True
    while still_working:
      if platform.system() == 'Windows':
          print("Running on a Windows system, Timeout recovery is not available")
      else:
        signal.alarm(signal_timer)
      
      try:
        row = [geometry, date]
    
        #We need a start date and an end date. Just like a regular python slice, 
        #the end date is not included, so by using a 1 day timedelta, I am actually limiting
        #the range to only the day in question and up to 7 days before
        start_date = date - timedelta(days = 7)
        end_date = date + timedelta(days = 1)
    
        #First I get the image collection from the MODIS data, filter it only to the days in question
        #and select my bands, then sort so the most recent day in the group is at the top
        DatedCollection = Collection.filter(ee.Filter.date(start_date, end_date)) \
                                    .filter(ee.Filter.notNull(['system:index'])) \
                                    .sort('system:index', False)
    
        #Because it is sorted, we can take the first
        point = DatedCollection.first().unmask(0)
        
        #We clip it to the geometry and pull the correct bands
        aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
    
        bands = point.reduceRegion(reducer = ee.Reducer.mean(),
        geometry= aoi)
    
        #Normalizes bands to between 0 and 1
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
    '''
    

    Parameters
    ----------
    geometry : shapely.geometry, string
        Geometry of the area for which we would like to pull an image.
    signal_timer : integer, optional
        Number of seconds to wait before timeout. The default is 5.

    Returns
    -------
    np.array
        Given a geometry, pulls the Copernicus DEM 30m resolution image for that
        geometry. Returns a 3D Numpy Array of Elevation, slope and aspect.

    '''
    
    #Checking for correct image type
    try:
        if not isinstance(geometry, shapely.geometry.base.BaseGeometry):
            geometry = wkt.loads(geometry)
    except Exception as e:
        print(e)
        print("You did not enter a geometry for the geometry variable")
        return
    
    
    #Here we establish a connection to planetary computer. We don't move on
    #until we have connected
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
      #implementing a timeout timer if system allows
      if platform.system() == 'Windows':
          print("Running on a Windows system, Timeout recovery is not available")
      else:
        signal.alarm(signal_timer)
      
        
      try:    
    
        # Adapted from https://planetarycomputer.microsoft.com/dataset/cop-dem-glo-90#Example-Notebook :
          if error_num == 0:
              #searches for an image that contains our bounding box
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
              #if we have multiple errors we scale up the request area and
              #look for an intersection.
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
              #The amount of scaling goes up with the number of timeouts until
              #an image is found
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
    
          
    
          #Here we bringin the asset
          items = list(search.get_items())
    
          signed_asset = planetary_computer.sign(items[0].assets["data"])
          
          data = (
              rioxarray.open_rasterio(signed_asset.href)
              .squeeze()
              .drop("band")
          )
    
          #We apply a mask to cover up everything but the exact area we want
          mask_lon = (data.x >= min_lon) & (data.x <= max_lon)
          mask_lat = (data.y >= min_lat) & (data.y <= max_lat)
    
    
          cropped_data = data.where(mask_lon & mask_lat, drop=True)
    
          #No we stack the image using datashader to give the elevation on the
          #red band, slope on the blue and aspect on the green
          img = stack(shade(cropped_data, cmap="red"), 
                      shade(xrspatial.slope(cropped_data), cmap="blue"),
                      shade(xrspatial.aspect(cropped_data), cmap="green"))
          
      
        
      except Exception as e:
        print("Copernicus Error, Trying Again")
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

    #bring everything together into an image and convert to a 3 channel np array
    img = img.to_pil()
    img=img.convert('RGB')
    return numpy.array(img).reshape(img.size[0], img.size[1], 3)



##########Sentinel 1 #########################################################

def pull_Sentinel1(geometry, date):
    '''
    

    Parameters
    ----------
    geometry : shapely.geometry, string
        Geometry of the area for which we would like to pull an image.
    date : date, string
        Date for which we would like to pull imagery. Can be a date or string
        defining the date.

    Returns
    -------
    np.array
        Returns a 3D numpy array in the shape of the image with one band of
        sentinel 1 on each channel

    '''
    
    #Testing for correct inputs
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
    
    #define area of interest by coordinates and construct a date range to look
    #over
    aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
    start_date = date - timedelta(days = 80)
    end_date = date + timedelta(days = 1)

    try:
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
        
        try:
            r = urllib.request.urlopen(url)
            img = Image.open(r)
            return numpy.array(img).reshape(img.size[0], img.size[1], 3)
        except Exception as e:
            print(e)

    except Exception as e:
        print(e)
        print("Have you logged in to google earth engine?")
        return


#############Sentinel 2a ########################################################
def pull_Sentinel2a(geometry, date):
    '''
    

    Parameters
    ----------
    geometry : shapely.geometry, string
        Geometry of the area for which we would like to pull an image.
    date : date, string
        Date for which we would like to pull imagery. Can be a date or string
        defining the date.

    Returns
    -------
    np.array
        Returns a 3D numpy array in the shape of the image with bands
        representing the geologic bands of sentinel 2 on each channel

    '''
    
    #Testing for correct inputs
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
    
    #define our area of interest by coordinates and construct a date range to
    #look over
    aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
    start_date = date - timedelta(days = 80)
    end_date = date + timedelta(days = 1)
    
    try:
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
            img = Image.open(r)
            return numpy.array(img).reshape(img.size[0], img.size[1], 3)
        except Exception as e:
            print(e)
    
    except Exception as e:
        print(e)
        print("Have you logged in to google earth engine?")
        return
    
    
##################Sentinel 2b################################

def pull_Sentinel2b(geometry, date):
    '''
    

    Parameters
    ----------
    geometry : shapely.geometry, string
        Geometry of the area for which we would like to pull an image.
    date : date, string
        Date for which we would like to pull imagery. Can be a date or string
        defining the date.

    Returns
    -------
    np.array
        Returns a 3D numpy array in the shape of the image with bands
        representing the vegetation from sentinel 2 on each channel

    '''
    
    #Testing for correct inputs
    
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

    #define area of interest by coordinates and define a date range to look in
    aoi = ee.Geometry.Polygon(list(geometry.exterior.coords))
    start_date = date - timedelta(days = 80)
    end_date = date + timedelta(days = 1)
    
    try:
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
    
    
        #pull in the image by url and save
    
        try:
            r = urllib.request.urlopen(url)
            img = Image.open(r)
            return numpy.array(img).reshape(img.size[0], img.size[1], 3)
        
        except Exception as e:
            print(e)
    
    except Exception as e:
        print(e)
        print("Have you logged in to google earth engine?")
        return
    
############# GRIDMET Weather Data##########################################
def pull_GRIDMET(geometry, date, num_days_back = 7):
    '''
    

    Parameters
    ----------
    geometry : shapely.geometry, string
        Geometry of the area for which we would like to pull an image.
    date : date, string
        Date for which we would like to pull data. Can be a date or string
        defining the date.
    num_days_back : integer, optional
        How many days worth of history to pull. The default is 7.

    Returns
    -------
    df : pd.DataFrame
        Returns a data frame with 7 columns, geometry, date, precipitation,
        wind direction in degrees, minimum temperature, maximum temperature,
        and wind velocity.

    '''
    
    #testing for correct input types
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
    
    #creating a data frame into which to copy our data.
    daterange = date_range(end = date, periods = num_days_back).tolist()
    geometries = [geometry for x in range(num_days_back)]
    
    df = DataFrame({"geometry":geometries, 
                    "date":daterange,
                    "precip":["NULL"]*num_days_back, 
                    "wind_dir":["NULL"]*num_days_back, 
                    "temp_min":["NULL"]*num_days_back, 
                    "temp_max":["NULL"]*num_days_back, 
                    "wind_vel":["NULL"]*num_days_back})
    ##Main loop that iterates over areas and data in the rows of the data frame


    for i in range(len(df)):
    #define area of interest by coordinates
        aoi = ee.Geometry.Polygon(list(df.geometry[i].exterior.coords))
        start_date = df.date[i] - timedelta(days=1)
        end_date = df.date[i]


        try:
            #Brining in the image collection from Google earth engine
            lst = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')\
                .filterDate(start_date, end_date)\
                .filterBounds(aoi)\
                .select('pr', 'th', 'tmmn', 'tmmx', 'vs')

            #Calculating 5 variables of interest
            precip = round(lst.mean().sample(aoi, 1000).first().get('pr').getInfo(),2)
            wind_dir = round(lst.mean().sample(aoi, 1000).first().get('th').getInfo(),2)
            temp_min = round(lst.mean().sample(aoi, 1000).first().get('tmmn').getInfo(),2)
            temp_max = round(lst.mean().sample(aoi, 1000).first().get('tmmx').getInfo(),2)
            wind_vel = round(lst.mean().sample(aoi, 1000).first().get('vs').getInfo(),2)
            
            
            df.iloc[i,2::] = [precip, wind_dir, temp_min, temp_max, wind_vel]
        
        except Exception as e:
            print(e)
            df.iloc[i,2::] = ['NULL', 'NULL', 'NULL', 'NULL', 'NULL']

    return df
        
#################### Stitch Our dataframe back into an image ############
def stitch_aso(reference, df, date):
    '''
    

    Parameters
    ----------
    reference : string
        String path to a reference image of the target tiff
    df : TYPE
        A dataframe in the same structure and order as the one put out by
        chop_aso().
    date : string, date
        A date to add to the image names.

    Returns
    -------
    im_array : np.array
        An array of the same size and shape as the original tiff but with the
        values of the SWE column of the dataframe. Also saves a csv and a .jpeg
        in the working directory

    '''
    geo_tiff = GeoTiff(reference)

    zarr_array = geo_tiff.read()

    im_array = np.array(zarr_array)
    
    pad_after0 = (20 - (im_array.shape[0] % 20))
    pad_after1 = (20 - (im_array.shape[1] % 20))
    im_array = np.pad(im_array, 
           pad_width = ((0,pad_after0), (0,pad_after1)), 
           constant_values = -9999)
    
    im_array[im_array == -9999] = np.nan
    im_array = im_array*0 +1
    
    
    
    y = 0

    for i in range(0, im_array.shape[1], 20):
        for j in range(0, im_array.shape[0], 20):
            im_array[j:j+20, i:i+20] = im_array[j:j+20, i:i+20] * df.SWE[y]
            y+=1


    ax = sns.heatmap(im_array, vmin = 0, vmax = 100, cmap = "mako_r", yticklabels=False, xticklabels=False)
    plt.savefig(f"{reference[0:-7]}_{date}_prediction.png")
    np.savetxt(f"{reference[0:-7]}_{date}_prediction.csv", im_array, delimiter=",")
           
    
    return im_array



############# testing###########

def testing():
    '''
    

    Returns
    -------
    None.
        Runs a test of all the data pulling functions based on a geometry and date

    '''
    
    polygon = Polygon([[-119.49, 38.22], [-119.59, 38.22], [-119.59, 38.32], [-119.49, 38.32], [-119.49, 38.22]])
    
    print("Testing Modis Image Pull")
    print(pull_MODIS_image(polygon, '2018-12-12', 'MOD10A1', buffer_percent=0.0))
    print("Testing Modis List Pull")
    print(pull_MODIS_list(polygon, '2020-10-10', 'MOD10A1'))
    print("Testing Copernicus Image Pull")
    print(get_copernicus(polygon))
    print("Testing Sentinel1 Image Pull")
    print(pull_Sentinel1(polygon, '2018-10-10'))
    print("Testing Sentinel2a Image Pull")
    print(pull_Sentinel2a(polygon, '2018-10-10'))
    print("Testing Sentinel2b Image Pull")
    print(pull_Sentinel2b(polygon, '2018-10-10'))
    print("Testing GRIDMET Weather Pull")
    print(pull_GRIDMET(polygon, '2018-10-10'))