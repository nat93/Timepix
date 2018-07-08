#########################################
#                                       #
# Stefano Petrucci - s.petrucci@cern.ch #
#                                       #
#########################################

#################################################################
#                                                               #
# USAGE: python this-script.py folder_date real_date			#
#                                                               #
# DATE FORMAT: yyyy_mm_dd                                       #
#                                                               #
#################################################################

#################################################################
#																#
# Write a description                                           #
#																#
#################################################################

import csv 
import glob
import re
import time
import datetime
import sys 
import numpy

# enlarging the file size of the single line otherwise sometime is impossible to read files
csv.field_size_limit(sys.maxsize)

numbers = re.compile(r'(\d+)')

# function for sort the list of file inside a folder in the right way
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

files = []
j = 0
old_time = 0
events_chip0 = 0
events_chip1 = 0

#for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/Detectors Data/Medipix*.csv'), key=numericalSort):
for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/Medipix/Medipix*.csv'), key=numericalSort):
    files.append(infile)

offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple())) * 1000 # in ms

for i in range(0, len(files)):
	with open(files[i]) as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
# TIME CONVERSION
			timestring = row['Time']
			time_ms = timestring.split(':')[-1]
			time_no_ms = re.sub(':\d+$', '', timestring)
			timestamp = int(time.mktime(datetime.datetime.strptime(time_no_ms, '%d %b %y %H:%M:%S').timetuple())) * 1000 + int(time_ms) - offset # in ms
# STARTING SOMETHING
			if j == 0:
				old_time == int(round(timestamp/1000))
			j += 1
			if int(round(timestamp/1000)) == old_time:
				if int(row['Chip']) == 0:
					events_chip0 += 1
				if int(row['Chip']) == 1:
					events_chip1 += 1
			if int(round(timestamp/1000)) != old_time:
				print int(old_time), events_chip0, events_chip1
				old_time = round(timestamp/1000)
				events_chip0 = 0
				events_chip1 = 0
				if int(row['Chip']) == 0:
					events_chip0 = 1
				if int(row['Chip']) == 1:
					events_chip1 = 1
	j = 0
