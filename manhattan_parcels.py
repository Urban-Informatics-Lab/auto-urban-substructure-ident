import numpy as np
import geopandas as gpd
from shapely.ops import cascaded_union
from utils.parcel import get_parcel_obj_list

try:
    cleaned_mn = gpd.read_file('data/nyc/mn')
    nycutmcrs = {'init': 'epsg:6347'}
except:
    full_data = gpd.read_file('data/nyc/pluto18v2.1')
    print('read in data...')

    # limit to manhattan
    mn = full_data.loc[full_data.Borough.isin(['MN'])]

    # exclude land uses:
    #  9: 'Open Space & Outdoor Recreation' # Parks 
    # 10: 'Parking Facilities', # Parking
    # 11: 'Vacant Land' # Vacant
    mn = mn.loc[mn.LandUse.astype(np.float).isin([1,2,3,4,5,6,8])]

    # covert to UTM coordinates
    nycutmcrs = {'init': 'epsg:6347'}
    mn.to_crs(nycutmcrs, inplace=True)

    # remove islands
    cleaned_mn = mn.loc[(mn.geometry.centroid.y > 4505050) &
                        ~((mn.geometry.centroid.x > 590500) & (mn.geometry.centroid.y < 4517500)) & 
                        ~mn.Address.apply(lambda x: 'loop' in x.lower() or 'main' in x.lower() 
                                          or 'river road' in x.lower())]

    cleaned_mn.to_file('data/nyc/mn_built_lots')
print('reduced data to Manhattan...')

# create parcel objects
mn_blocks = cascaded_union(cleaned_mn.geometry)
print('identified block polygons...')
parcel_objs = get_parcel_obj_list(cleaned_mn.geometry, mn_blocks, nycutmcrs, min_dist=4,
                                  struct_yr_built=cleaned_mn.YearBuilt)
print('created parcel objects...')
parcel_objs.to_pickle('data/nyc/manhattan_parcels.pkl.zip', 'zip')
print('done!')