# Overview

This is a package to help automate the pulling of Satellite Imagery and Tabular data that is useful for the training of neural networks to predict Snow Water Equivalent. It can also be easily adapted to pull imagery for any of the areas these satellite assets cover.

This package was created to solve a specific task, but could easily prove useful for other projects. If you find you would like to modify, add, change, or fix any function, feel free to fork this repo and make your own changes!

# Installation

There are many dependencies for this package, the most problematic of which is GDAL. Attempt to install the package first via: 

```
pip install snowcast-Malachyiii
from snowcast import data_wrangling
```

If this does not work, the most likely solution is to first install the dependencies yourself.

Functions can be used by prepending the `data_wrangling.<function_name>` form

## Windows Users

If you would like to install this package to a Windows machine, you **must** first install GDAL from a wheel. GDAL and Windows don't play nice unfortunately.

# Functions

This package contains:

## chop_aso

The intent of chop_aso is to take in a .tiff file, and cut it up into a dataframe of 1kmx1km geometries that can be fed into the later functions. It could easily be modified for any .tif file or target variable as long as you know two things:

1. The pixel resolution of the image
2. The size of the target geometry

It is meant to work with tiffs that are lat/long, but could also be modified to take other projections.

`chop_aso(tiff_image_name, groundtruth =  False):`
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

## pull_MODIS_image

Not used in our actual prediction module, but useful so we decided to keep it. This pulls the imagery for either of the two MODIS satellites (MOD10A1 and MYD10A1) and outputs the image from that satellite, cropped to the given geometry. If it cannot find an image from the specified day, it looks up to 7 days back

`pull_MODIS_image(geometry, date, modis, buffer_percent = 0.0):`
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
## pull_MODIS_list

This pulls the imagery for either of the two MODIS satellites (MOD10A1 and MYD10A1) and converts the data from that image, cropped to the given geometry, outputting the result as a table. If it cannot find an image from the specified day, it looks up to 7 days back

`pull_MODIS_list(geometry, date, modis, signal_timer = 5):`
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

## get_copernicus

This function pulls an image from the Copernicus 30m resolution DEM for the specified geometry. A date is not necessary because the DEM does not have date dependent image collections. The initial pull brings in the elevation, and the function datashader to create a slope and an aspect layer, outputing all three of these as a 3 channel numpy array

`get_copernicus(geometry, signal_timer =5):`
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
## Sentinel Series

The Sentinel functions are all roughly the same. The look on google earth engine for an image of the given geometry, on the specific satellite from the particular day, looking backwards in time if an image isn't available that day. The output is either all the bands (sentinel 1) or a subset of bands that are specifically useful for snow forecasting. In the case of Sentine2a thats geologic features, and in the case of Sentinel2b that's vegetation features. Any of these functions could be easily modified to pull in other bands as necessary.
 
`pull_Sentinel1(geometry, date):`
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
 
`pull_Sentinel2a(geometry, date):`
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

`pull_Sentinel2b(geometry, date):`
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

## pull_GRIDMET

This function is intended to be used to create an LSTM style neural network. It pulls a data frame of the weather data from a specific geometry for given number of days before the supplied date. This function could easily be modified to pull in other weather data, as there are many variables available.

`def pull_GRIDMET(geometry, date, num_days_back = 7):`
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
## stitch_aso
 
This function takes the dataframe we have created, and stitches it back together into an image and a 1 channel numpy array of predictions. It is very important that the order of the dataframe not be changed, as it depends on math rather than geometries to work. All points not inside the shape in question should be coded as numpy nan's. Any changes made in chop_aso to account for resolution differences **must also be made here**


`stitch_aso(reference, df, date):`
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

## testing

Testing provides an all in one function to test your image pulls to ensure your google earth engine and planetary computer connections are working

testing():
    '''
    

    Returns
    -------
    None.
        Runs a test of all the data pulling functions based on a geometry and date

    '''
    
# Wrap-Up

If you have any questions about this module, or wish to use it for your project, please email me at malachy.j.moran@gmail.com and cite this project in your paper's works cited.
