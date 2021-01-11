import os
from pymongo import MongoClient
from plpygis import Geometry
from datetime import datetime
from collections import defaultdict

client = MongoClient('localhost', 27017)
traj_db = client.traj_db
dsta_collection = traj_db.dsta_collection

def xreadlines(f):
    while True:
        line = f.readline()
        if line !='':
            yield line
        else:
            break

def range_enumerate(iter, start, end=None):
    if end is None:
        end = start
        start = 0
    for i, v in zip(range(start, end), iter):
        yield i, v


file_path = '../DSTA_dataset/'

# def append_mongodb(collection, dict):
#     for mmsi, traj in dict.items():
#         collection.update(
#             {"_id": mmsi}, 
#             {"$push": {"traj" : {"$each": traj} } }, 
#             upsert=True
#         )

def append_mongodb(collection, toks):
    # schema: time_utc|geom|mmsi_hash|speed|heading|flag_hash|breadth_extream|draught_large|grt|dwt|ship_type_code|updt_dt
    lng = Geometry(toks[1]).shapely.x
    lat = Geometry(toks[1]).shapely.y
    timestamp = datetime.strptime(toks[0], r'%Y-%m-%d %H:%M:%S').timestamp()
    mmsi = toks[2]
    speed = toks[3]
    heading = toks[4]
    type_code = toks[10]

    # print(lat, lng, timestamp, mmsi, speed, heading, type_code)

    collection.insert({
        'lat': lat, 
        'lng': lng, 
        'timestamp': timestamp, 
        'tid': mmsi, 
        'speed': speed,
        'heading': heading, 
        'type': type_code
    })

#references: read data in May/June/July
ref_filenames = ['geo5_w_static_Aug2017.csv']
collection_names = ['dsta_May', 'dsta_Jun', 'dsta_Jul', 'dsta_Aug']
# ref_filenames = ['geo5_w_static_May2017.csv', 'geo5_w_static_Jun2017.csv', 'geo5_w_static_Jul2017.csv', 'geo5_w_static_Aug2017.csv']

traj_dict = defaultdict(list)

cnt = 0
for filename_, collection_name in zip(ref_filenames, collection_names):
    filename = os.path.join(file_path, filename_)
    print('filename=', filename)
    with open(filename) as f:
        for i, line in enumerate(xreadlines(f)):
            cnt += 1
            if cnt < 841978262 - 740112496 + 1:
                continue
            # print(i, line.split('|'))
            toks = line.split('|')
            if i>=1:
                append_mongodb(dsta_collection, toks)
            if i%100000==0:
                print('i=', i)

#tata cable: Aug 1st to Aug 15th
