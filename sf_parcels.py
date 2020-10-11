import pandas as pd
import geopandas as gpd
from shapely.ops import cascaded_union
import shapely.wkt
from utils.parcel import get_parcel_obj_list


try:
    biggened_parcels = gpd.read_file('data/sf/sf_built_parcels')
    print('read in clean data...')
except:
    full_data = pd.read_csv('data/sf/LandUse2016.csv')
    full_data.the_geom = full_data.the_geom.apply(shapely.wkt.loads)
    full_data = gpd.GeoDataFrame(full_data, geometry='the_geom')
    full_data.crs = {'init': 'epsg:4326'}
    print('read in data...')

    no_parks = full_data[~full_data.LANDUSE.isin(['OpenSpace', 'Right of Way', 'VACANT'])]
    # get rid of treasure island
    no_parks = no_parks[~((no_parks.the_geom.centroid.y > 37.8) & (no_parks.the_geom.centroid.x > -122.39))]
    #get rid of a lake merced development
    no_parks = no_parks[~((no_parks.the_geom.centroid.y < 37.725) & (no_parks.the_geom.centroid.x < -122.49))]
    no_parks = no_parks[(no_parks.the_geom.centroid.x < -122.35)]
    print('excluded land uses and reduced area...')

    biggened_parcels = no_parks.the_geom.apply(lambda shape: shape.buffer(0.000001))
    biggened_parcels = biggened_parcels.to_crs({'init': 'epsg:7131'})
    print('converted crs...')
    
    biggened_parcels.to_file('data/sf/sf_built_parcels')

no_park_blocks = cascaded_union(biggened_parcels)
print('identified block polygons...')

parcel_objs = get_parcel_obj_list(biggened_parcels, no_park_blocks, biggened_parcels.crs, min_dist=4,
                                  struct_yr_built=no_parks.YRBUILT)
print('created parcel objects...')

parcel_objs.to_pickle('data/sf/sf_parcels.pkl.zip', 'zip')
print('done!')