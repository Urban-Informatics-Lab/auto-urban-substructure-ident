import pandas as pd
import geopandas as gpd
from shapely.geometry import *
from shapely.ops import cascaded_union
from utils.parcel import get_parcel_obj_list


def split_polygon(parcel_geometry):
    if type(parcel_geometry) == MultiPolygon:
        return list(parcel_geometry)
    elif type(parcel_geometry) != Polygon:
        raise Exception(type(parcel_geometry))
    else:
        return [parcel_geometry]


houston_utm_crs = {'init': 'EPSG:32615'}

try:
    no_mp_data = gpd.read_file('data/houston/innerloop_built_parcels')
except:
    innerloop = gpd.read_file('data/houston/innerloop.geojson')
    innerloop.to_crs(houston_utm_crs, inplace=True)
    innerloop = innerloop.iloc[0, 0]

    full_data = gpd.read_file('data/houston/houston_lu2018.geojson')
    print('read in data...')
    harris_data = full_data.loc[full_data.COUNTY=='HARRIS']
    noshape_mask = harris_data.geometry.apply(lambda x: type(x) is not Polygon and
                                               type(x) is not MultiPolygon)
    clean_harris_data = harris_data.loc[~noshape_mask]
    clean_harris_data.to_crs(houston_utm_crs, inplace=True)
    print('converted crs...')
    in_beltway = clean_harris_data.geometry.apply(lambda parcel: parcel.centroid.within(innerloop))

    beltway_data = clean_harris_data.loc[in_beltway]
    
    vacant_county = ['1000', '1002', '2000', '2002', '4200', '4300', '4400', '4600', '4700', '4730', '4740', 
                     '4750', '4760', '4770', '7000', '7300', '9999', 'SF5', 'SF6', '4510', '4530']
    parks_county = ['4930', '4387', '4389', '4520', '2004', '9200', '9920']
    parking_county = ['4338', '4339', '4372', '4329', '4395', '4368']
    infra_county = ['4762', '1105', '1106', '4371', '7001', '4720', '4500', '4752', '4900', '4765']
    farms_county = ['9011', '9047', '9910', '9063', '9003', '9007', '9015', '9019', '9023', '9027', '9031', '9035', '9039',
                    '9043', '9051', '9055', '9059', '9067', '9071', '9075', '9079', '9083', '9100']
    questionable_county = ['4546', '4544']
    
    cleaned_beltway_data = beltway_data.loc[~beltway_data.LANDUSE_CD.isin(vacant_county + parks_county + farms_county
                                                                      + parking_county + infra_county + 
                                                                      questionable_county)]
    print('excluded land uses and reduced area...')
    
    no_mp_data = pd.DataFrame(cleaned_beltway_data.geometry.apply(split_polygon).tolist(), 
                          index=cleaned_beltway_data.STATE_CLASS).stack()
    print('split multipolygons...')
    no_mp_data = no_mp_data.reset_index()[[0, 'STATE_CLASS']]
    no_mp_data.columns = ['geometry', 'STATE_CLASS']
    no_mp_data = gpd.GeoDataFrame(no_mp_data, geometry='geometry', crs=cleaned_beltway_data.crs)
    
    
    no_mp_data.to_file('data/houston/innerloop_built_parcels')
print('reduced data to inner loop...')
    
biggened_parcels = no_mp_data.geometry.apply(lambda shape: shape.buffer(0))
blocks = cascaded_union(biggened_parcels)
print('identified block polygons...')
parcel_objs = get_parcel_obj_list(biggened_parcels, blocks, biggened_parcels.crs, min_dist=4,
                                  land_use=no_mp_data.STATE_CLASS)
print('created parcel objects...')
parcel_objs.to_pickle('data/houston/inner_houston_parcels.pkl.zip', 'zip')
print('done!')