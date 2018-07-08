#########################################
#										#
# Stefano Petrucci - s.petrucci@cern.ch	#
#										#
#########################################

import csv
import numpy
import glob
import numpy as n
import re
import time
import datetime
import sys

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

channel_names = []
timestamps = []
files_dictionary = {}
matrix_dict = {}
folder = []

for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/Detectors Data/ScintillatorsMeasurements_*.csv'), key=numericalSort):
    folder.append(infile)

offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000

j = 0
k = 0
#for i in range(0, len(folder)):
for i in range(1380, 1394):
	with open(folder[i]) as csvfile:
#		print folder[i]
		reader = csv.DictReader(csvfile)
#		if i == 0:
#		if i == 800:
		header = reader.fieldnames
		header.remove('')
#		print header
		for row in reader:
			if channel_names == []:
				channel_names = row[header[4]].split('|')
				# removing empty string from the list otherwise len(channel_names) is wrong
				channel_names = filter(None, channel_names)
#				channel_names[:] = (value for value in channel_names if value != 'n.c.' )
			time_vector = list(map(int, row[header[10]].split('|')))
			for l in range(0, len(time_vector)):
				time_vector[l] = str(int(time_vector[l])/10**6 - int(offset))
			timestamps.extend(time_vector)
			for items in header:
				if items != 'Time' and items != 'acqTimestamp' and items != 'timestamps' and items != 'channelNames':
					temp_values = (list(map(float, row[items].split('|'))))
					matrix_dict[items] = n.reshape(temp_values,(len(channel_names), len(time_vector)))
					for name, meas in zip(channel_names, matrix_dict[items]):
						if name != 'n.c.':
							files_dictionary[name] = open('/home/petruccs-local/scripts/'+sys.argv[2]+'/wc_meas/'+name+'_'+items+'.dat', 'a')
#							if j == 0 and i == 0 and k == 0:
#								files_dictionary[name].write("timestamps,"+items+"\n")
							j = 1
							files_dictionary[name].write('\n'.join('%s %s' % x for x in zip(time_vector, meas)))
							files_dictionary[name].write('\n')
							files_dictionary[name].close()
						j = 0
			k = 1
