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
offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000
#motors_name = ['TACW.51998'] #TAC
#motors_name = ['TCSM.51934'] #LHC collimator
motors_name = ['BSHV.51793', 'TAL.52196.LIN1', 'TAL.52196.LIN2', 'TAL.52196.LIN3', 'TCXHW.51651.H1', 'TCXHW.51651.H2', 'TQCD.51794', 'TQCD.51991', 'XRPH.51937.H1', 'XRPH.51937.H2', 'XRPH.52202.H1', 'XRPH.52202.H2', 'TCPC.51795_BP', 'TCPC.51795_LIN', 'TCPC.51795_ROT', 'TQCD.201271']

############ To add the control for argv[1]. Cfr the csv reader -> root file

for motors in motors_name:
	print 'Reading: ', motors
	for infile in sorted(glob.glob(path + motors+'_*.csv'), key=numericalSort):
		files.append(infile)
	filename = './'+sys.argv[1]+'/'+sys.argv[2]+'_'+motors.lower().replace('.', '_')+'.dat'
	with open(filename, 'w') as output_file:
		for i in range(0, len(files)):
			with open(files[i]) as csvfile:
				reader = csv.DictReader(csvfile)
				for row in reader:
					#output_file.write(str(float(int(row['acqStamp'])*1E-6 - offset)/1000)+' '+row['position']+'\n') #TAC
					#output_file.write(str(float(int(row['acqStamp'])*1E-6 - offset)/1000)+' '+row['lvdt_left_downstream']+' '+row['lvdt_left_upstream']+' '+row['lvdt_right_downstream']+' '+row['lvdt_right_upstream']+'\n') #LHC collimator
					output_file.write(str(float(int(row['acqStamp'])*1E-6 - offset)/1000)+' '+row['controllerPosition']+' '+row['lvdtPosition']+' '+row['resolverPosition']+'\n')
	files = []











