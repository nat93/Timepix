# Robin Landstrom (per.robin.landstrom@cern.ch) 
import glob

def getFesaFiles():
    PATH = "/home/petruccs-local/scripts/robin/data/"
    FILE_BASE = "CpFMMeasurements_"
    FILE_IDS = [str(id) for id in range (483,500)] 
    SUFFIX = ".csv"    
    return (PATH + FILE_BASE + id + SUFFIX for id in FILE_IDS)

def getLocalFiles():
    PATH = "/home/petruccs-local/scripts/robin/data/" #Events from the same time copied from PCEN35019
    FILE_BASE = "measurement-"    
    SUFFIX = ".dat"    
    return glob.glob(PATH + FILE_BASE + "*" + SUFFIX)

def chunks(l, n):
    """Yield successive n-sized chunks from l, used to parse flattened 2d arrays from rda logger"""
    for i in range(0,len(l), n):
        yield l[i:i+n]

#Used for parsing rda logger files
data_arrays_2d = { # Column name and index (from 0)
    'amplitude' : 2,
    'baseline' : 3,
    'charge' : 5,
    'fallingEdgeTime' : 6,
    'peakTime' : 7,
    'peaks' : 8,
    'risingEdgeTime' : 9    
}

#Used for parsing rda logger files
data_arrays_1d = {
    'channelNames' : 4,
    'timestamps': 10
}

#Used for parsing local files
column_order = ["id", "channelNames", "eventId", "baseline", "amplitude", "peakTime", "charge", "risingEdgeTime", "fallingEdgeTime", "peaks"]

#Return a list of lists of channel dicts from a row from a rda logger file
def parse_row(row):
    number_of_samples = len(row[data_arrays_1d['timestamps']].split("|"))
    number_of_channels = 32

    result = [[{} for _ in range(number_of_channels)] for _ in range(number_of_samples)]

    #Populate item with 2D arrays
    for key, val in data_arrays_2d.items():
        #print("Processing " + key + " on column " + str(val))
        generator = chunks(row[val].split("|"), number_of_samples)
        for channel_index, channel_data in enumerate(generator):
            #print("    channel index:" + str(channel_index))
            for sample_index, d in enumerate(channel_data):
                result[sample_index][channel_index][key] = float(d)

    #Add timestamps
    timestamps = row[data_arrays_1d['timestamps']].split("|")
    for item, timestamp in zip(result, timestamps):
        for channel in item:
            channel['timestamp'] = int(timestamp)

    #Add channel names
    channelNames = row[data_arrays_1d['channelNames']].split("|")
    for item in result:
        for channel, name in zip(item, channelNames):
            channel['channelName'] = name
    
    return result            

def parseFesaFiles():
    result = []
    for filename in getFesaFiles():
        f = open(filename, "r")
        #print("Parsing file: " + filename)
    
        headers = f.readline().split(",") # First line in the file is the header
#        print(headers)

        for line in f: # Each line in the file is a snapshot of the AcquisitionProperty of the FESA class taken every second, normally this contains 42-45 samples
            data = line.split(",")
            timestamps = data[data_arrays_1d['timestamps']].split("|") #timestamps is the last field but the lines in the file ends with a ',' therefor timestamp is the next to last field (-2)
            number_of_samples = len(timestamps) #To determine the number of samples on this row, look at the number of timestamps
            result += parse_row(data)
    return result


# Yield a list of lines related to one event in the measurement files
def event_lines(f):
    f.readline()
    f.readline() #Skip two header lines    
    header = f.readline()
    res = []
    for line in f.readlines():
        if (line.startswith("===")):
            res[:0] = [header]
            header = line
            #print(res)
            yield res
            res = []
            continue
        res += [line]
    res[:0] = [header]
    yield res	

#Return a list of lists of channel dicts from a measurement file
def parse(data_file):
    events = []
    f = open(data_file, "r")
    for event_data in event_lines(f):
        event_header = event_data[0].split(" ")
        ts = event_header[5]
        event = []
        for channel_line in event_data[1:]:
            channel_line = channel_line.split(" ")
            ch = {}
            for col, val in zip(column_order, channel_line):
                if col in data_arrays_2d:
                    ch[col] = float(val)
                else:
                    ch[col] = val
            ch["timestamp"] = int(ts) * 1000 #File ts in us, fesa ts in ns
            event += ch,            
        events += event,
    return events

def parseLocalFiles():
    result = []
    for filename in getLocalFiles():
        #print("Parsing file: " + filename)
        result += parse(filename)
    return result

def allignEvents(a, b):
    while 1:
        a_val = a.next()
        b_val = b.next()

        #result = []
        skip_a = 0
        skip_b = 0
        a_ahead = False
        b_ahead = False
        while a_val[0]["timestamp"] < b_val[0]["timestamp"]:
            a_val = a.next()
            skip_a += 1
            a_ahead = True

        while b_val[0]["timestamp"] < a_val[0]["timestamp"]:
            b_val = b.next()
            skip_b += 1
            b_ahead = True

        if(a_ahead):
            print("Skipped", skip_a, "values from a to sync timestamps")
        if(b_ahead):
            print("Skipped", skip_b, "values from b to sync timestamps")

        yield (a_val, b_val)
        

sorted_local = sorted(parseLocalFiles(), key=lambda x:x[0]["timestamp"])
sorted_fesa = sorted(parseFesaFiles(), key=lambda x:x[0]["timestamp"])

data = list(allignEvents(iter(sorted_local), iter(sorted_fesa)))
#data now contains a list of event tuples,
#the first item in the tuple contains the values from local files, the second from the rda logger
#each event is a list of channels
#each channel is a dictionary of values

allowed_error = 0.0005
print("Validating", len(data), "timestamp alligned events, allowed error", allowed_error)
for field in data_arrays_2d: #								Local file values - RDA logged values 
    print(field, "verified", all(all(map(lambda x:abs(x[0][channel][field] - x[1][channel][field]) < allowed_error, data)) for channel in range(3)))
	
allowed_error = 0.00001
print("Validating", len(data), "timestamp alligned events, allowed error", allowed_error)
for field in data_arrays_2d: #								Local file values - RDA logged values 
    print(field, "verified", all(all(map(lambda x:abs(x[0][channel][field] - x[1][channel][field]) < allowed_error, data)) for channel in range(3)))

