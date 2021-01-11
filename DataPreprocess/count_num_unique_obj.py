import os


def xreadlines(f):
    while True:
        line = f.readline()
        if line !='':
            yield line
        else:
            break

file_path = '../DSTA_dataset/'
filenames = ['geo5_w_static_May2017.csv', 'geo5_w_static_Jun2017.csv', 'geo5_w_static_Jul2017.csv', 'geo5_w_static_Aug2017.csv', 'geo5_w_static_Sep2017.csv']



obj_dict = {}

cnt = 0
for filename_ in filenames:
    filename = os.path.join(file_path, filename_)
    print('filename=', filename)
    with open(filename) as f:
        for i, line in enumerate(xreadlines(f)):
            cnt += 1
            if cnt %1000000==0:
                print(cnt)
            toks = line.split('|')
            mmsi = toks[2]
            obj_dict[mmsi] = True

    print('len=', len(obj_dict))

#tata cable: Aug 1st to Aug 15th
