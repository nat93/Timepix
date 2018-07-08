#########################################
#										#
# Stefano Petrucci - s.petrucci@cern.ch	#
#										#
#########################################

#######################################################
#                                                     #
# USAGE: python this-script.py yyyy_mm_dd yyyy_mm_dd  #
#													  #
#######################################################

import csv
import glob
import re
import time
import datetime
import sys

# enlarging the file size of the single line otherwise sometime is impossible to read files
csv.field_size_limit(sys.maxsize)

numbers = re.compile(r'(\d+)')

# function for sort the list of file inside a folder in the right way
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

path = '/home/anatochi/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/Motors Data/'
#path = '/home/petruccs-local/sps/MD_'+sys.argv[1]+'/Motors Data/'
files = []				# list of the files to open
gonio_name = ['TECS.51797', 'TECS.51652']
offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000

############ To add the control for argv[1]. Cfr the csv reader -> root file

for gonio in gonio_name:
	print 'Reading gonio: ', gonio
	for infile in sorted(glob.glob(path + gonio+'_*.csv'), key=numericalSort):
	    files.append(infile)
	filename = './'+sys.argv[1]+'/'+sys.argv[2]+'_'+gonio.lower().replace('.', '_')+'.dat'
	#filename = sys.argv[2]+'_'+gonio.lower().replace('.', '_')+'.dat'
	with open(filename, 'w') as output_file:
		for i in range(0, len(files)):
			with open(files[i]) as csvfile:
				#print files[i]
				reader = csv.DictReader(csvfile)
				for row in reader:
					output_file.write(str(float((int(row['acqStamp'])*1E-6 - offset)/1E3))+' '+row['axisAPosition_ctrl']+' '+row['axisAPosition_lvdt']+' '+row['axisBPosition_ctrl']+' '+row['axisBPosition_lvdt']+' '+row['crystalAAngle_ctrl']+' '+row['crystalAAngle_lvdt']+' '+row['crystalBAngle_ctrl']+' '+row['crystalBAngle_lvdt']+'\n')
	files = []
