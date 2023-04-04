# %% Import libraries
import geopandas as gpd
import pandas as pd
import requests
from bs4 import BeautifulSoup
pd.set_option('display.max_columns', None)
# %% Load Study Area
study_area = gpd.read_file("cnddb_data/office.geojson")

# transform to UTM
study_area = study_area.to_crs("EPSG:3310")  # Same CRS as CNDDB

# %% load quads and find all touching quads
quads = gpd.read_file("cnddb_data/quad75.shp")

# find intersecting quad(s) with study area
center_quad = quads[quads.intersects(study_area.unary_union)]

# find all quads touching the center quad
surronding_quads = quads[quads.touches(center_quad.unary_union)]

# combine center_quad and surronding_quads into one geodataframe
search_quads = gpd.GeoDataFrame(pd.concat([center_quad, surronding_quads]))

search_area = search_quads.to_crs("EPSG:4326")  # convert to WGS84 for GBIF
search_area = search_area.unary_union.wkt  # convert to WKT for GBIF
# TODO add option to include additional quads by name or code
# %% Load CNDDB data and clean
cnddb = gpd.read_file("cnddb_data/cnddb.shp", mask=search_quads)


def assign_taxon_category(row):
    plants = ['Dicots', 'Monocots', 'Bryophytes',
              'Gymnosperms', 'Herbaceous', 'Lichens', 'Ferns']

    animals = ['Reptiles', 'Mammals', 'Birds', 'Fish', 'Mollusks',
               'Crustaceans', 'Insects', 'Amphibians', 'Arachnids']

    if row['TAXONGROUP'] in plants:
        return 'Plant'
    elif row['TAXONGROUP'] in animals:
        return 'Animal'
    else:
        return 'Other'


cnddb['TABLECATEGORY'] = cnddb.apply(assign_taxon_category, axis=1)

cnddb = cnddb[cnddb['TABLECATEGORY'] != 'Other']
# %% Load CNPS data and filter by search quads


def fetch_CNPS_table(search_quads):
    quad_string = '&quad='
    quad_list = search_quads['QUADCODE'].tolist()
    for quad in quad_list:
        quad_string = quad_string + quad + ':'
    search_url = 'https://rareplants.cnps.org/Search/result?&crpr=1B:2B:4' + quad_string

    for i in range(3):
        try:
            response = requests.get(search_url)
            if response.status_code == 200:
                break
        except requests.exceptions.RequestException as e:
            print(e)
            continue

    soup = BeautifulSoup(response.text, 'html.parser')
    table = soup.find('table', attrs={'id': 'resultList'})
    rows = table.find_all('tr')
    table = soup.find('table', attrs={'id': 'resultList'})
    rows = table.find_all('tr')
    headers = [header.text for header in rows[0].find_all('th')]
    headers = [header.strip() for header in headers]
    data = []
    for row in rows[1:]:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        data.append([ele for ele in cols])
    return pd.DataFrame(data, columns=headers)


cnps = fetch_CNPS_table(search_quads)
# %% Build search list for GBIF by combining CNDDB and CNPS data
cnddb_list = cnddb['SNAME'].unique().tolist()
# remove "pop. X" where X can be any number in CNDDB list but keep rest of name
cnddb_list = [name.split('pop. ')[0].strip() for name in cnddb_list]

search_list = cnps['Scientific Name'].unique().tolist()

for name in cnddb_list:
    if name not in search_list:
        search_list.append(name)

# %% Search GBIF for species in search_list


data = []
for name in search_list:
    for i in range(3):
        print(name)
        try:
            search_url = 'http://api.gbif.org/v1/occurrence/search?scientificName={name}&geometry={search_area}&limit=1000'.format(
                name=name, search_area=search_area)
            response = requests.get(search_url)
            if response.status_code == 200:
                break
        except requests.exceptions.RequestException as e:
            print(e)
            continue
    try:
        gbif_data = response.json()
        # include support for multiple pages of results
        while not gbif_data['endOfRecords']:
            gbif_df = pd.json_normalize(gbif_data['results'])
            gbif_df['search_name'] = name
            data.append(gbif_df)
            offset = gbif_data['offset'] + gbif_data['limit']
            search_url = 'http://api.gbif.org/v1/occurrence/search?scientificName={name}&geometry={search_area}&limit=1000&offset={offset}'.format(
                name=name, search_area=search_area, offset=offset)
            response = requests.get(search_url)
            gbif_data = response.json()
        gbif_df = pd.json_normalize(gbif_data['results'])
        gbif_df['search_name'] = name
        data.append(gbif_df)
    except Exception as e:
        print(f'failed to retrieve data for: {name} - {e}')
        continue
# %% Combine GBIF data into one geoDataFrame
gbif = pd.concat(data)
gbif = gpd.GeoDataFrame(gbif, crs="EPSG:4326", geometry=gpd.points_from_xy(
    gbif.decimalLongitude, gbif.decimalLatitude))
gbif
# %%
gbif_cols = ['search_name',
             'basisOfRecord',
             'kingdom',
             'decimalLongitude',
             'decimalLatitude',
             'coordinateUncertaintyInMeters',
             'eventDate',
             'recordedBy',
             'informationWithheld',
             'gbifID',
             'occurrenceID',
             'catalogNumber',
             'institutionCode',
             'identificationRemarks',
             'occurrenceRemarks',
             'locality',
             'habitat',
             'locationRemarks',
             'georeferenceRemarks',
             'geometry']

gbif = gbif[gbif_cols]
gbif
# %% get counts of each species in gbif
gbif_counts = gbif['search_name'].value_counts().reset_index()
gbif_counts.columns = ['search_name', 'count']
gbif_counts
# %%
