#########################################
#										#
# Stefano Petrucci - s.petrucci@cern.ch	#
#										#
#########################################

#############################################################
#															#
# USAGE: python this-script.py bct_number folder_date date	#
#															#
# Date format: yy_mm_dd										#
#															#
#############################################################

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

files = []				# list of the files to open

############ To add the control for argv[1]. Cfr the csv reader -> root file

for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/Motors Data/TACW.51998_*.csv'), key=numericalSort):
    files.append(infile)

offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000
filename = sys.argv[2]+'_tacw_51998.dat'

with open(filename, 'w') as output_file:
	for i in range(0, len(files)):
		with open(files[i]) as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				output_file.write(str((int(row['acqStamp'])*1E-6 - offset)/1E3)+' '+row['position']+'\n')
