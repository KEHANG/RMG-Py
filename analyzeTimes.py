import numpy as np
import json

def loadTimesFromJSON(filename):
    data = []
    with open(filename, 'r') as f_in:
        for line in f_in:
            data.append(json.loads(line))
    return data

def reformJSONData(data):

    # get keys in data
    firstRecord = data[0]
    keys = firstRecord.keys()

    # initialize reformulated data
    data_ref = {}
    for key in keys:
        data_ref[key] = [firstRecord[key]]

    # packaging reformulated data
    for record in data[1:]:
        for key in keys:
            data_ref[key].append(record[key])

    return data_ref

def reject_outliers(data, m = 2.):
    data = np.array(data)
    d = abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return s<m

def reject_outliers_forMultipleKeys(data_ref, m = 2.):

    keys = data_ref.keys()
    # reject outliers
    truthIndices = np.array([True]*len(data_ref[keys[0]]))
    for key in keys:
        truthIndices *= reject_outliers(data_ref[key], m)

    data_ref_remain = {}
    for key in keys:
        data_ref_remain[key] = np.array(data_ref[key])[truthIndices]
        print '***************'
        print '***' + key + '***'
        print '***************'
        print 'Filtered data: ' + str(data_ref_remain[key])
        print 'Mean: '+str(np.mean(data_ref_remain[key]))
        print 'Median: '+str(np.median(data_ref_remain[key]))
        print 'St. Dev: '+str(np.std(data_ref_remain[key]))
    return data_ref_remain

data = loadTimesFromJSON('times.json')
data_ref = reformJSONData(data)
data_ref_remain = reject_outliers_forMultipleKeys(data_ref)


