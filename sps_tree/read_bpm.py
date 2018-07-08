import csv
import glob
import os
import re
import ntpath
import sys
import time
import datetime

numbers = re.compile(r'(\d+)')

# function for sort the list of file inside a folder in the right way
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

files = []
timestamps = []
positions = []
offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000
bpm_list = ['SPS_BPMCOL_BOX1_HOR_UP', 'SPS_BPMCOL_BOX1_HOR_DOWN']
j = 1
old_ts = 0
i = 0
timestamp = 0
meas = 0

for bpm in bpm_list:
	for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/BI instrumentation/'+bpm+'_*.csv'), key=numericalSort):
		files.append(infile)
	with open(sys.argv[2]+'_'+bpm.lower()+'.dat', 'w') as output_file:
		print 'Opening bpm: ', sys.argv[2]+'_'+bpm.lower()+'.dat'
		for i in range(0, len(files)):
			with open(files[i]) as csvfile:
				reader = csv.DictReader(csvfile)
				for row in reader:
					timestamp = (int(row['acqTimestamp'])/1000000 - offset)/1000
					if i == 0 or timestamp < 0:
						old_ts = timestamp
						i = 1
					if timestamp == old_ts:
						meas += float(row['horLinearPos'])
						j += 1;
					if timestamp == old_ts + 1:
						old_ts = timestamp
						j = 1
						output_file.write(str(timestamp)+' '+str(meas/j)+'\n')
						meas = float(row['horLinearPos'])
	files = []


#					timestamps.append(str(int(row['acqTimestamp'])/1000000 - offset))
#					print str((int(row['acqTimestamp'])/1000000 - offset)/1000), j
#					positions.append(row['horLinearPos'])
#					if old_ts == (int(row['acqTimestamp'])/1000000 - offset)/1000:	
#						j += 1
#					else:
#						j = 0
#		for time, meas in zip(timestamps, positions):
#			output_file.write(time+' '+meas+'\n')
