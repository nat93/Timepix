#########################################
#										#
# Stefano Petrucci - s.petrucci@cern.ch	#
#										#
#########################################

#################################################################
#																#
# USAGE: python this-script.py folder_date real_date bct_number #
#																#
# DATE FORMAT: yyyy_mm_dd										#
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

files = []				# list of the files to open
timestamp = 0			# the initial timestamp
list_timestamp = []		# list with timestamps
ts_exponent = []		# list with timestamps unit exponent
intensity_exp = 0		# exponent of the intensity (10**intensity_exponent)
intensity = []			# list with intensities
number_of_meas = 0		# number of measurements inside every meas vector
len_array = 0			# new lenght of the array, computed to be multiple of 200 (see the code)

############ To add the control for argv[1]. Cfr the csv reader -> root file

for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/BI instrumentation/BI_'+sys.argv[3]+'_BCTDC_*.csv'), key=numericalSort):
    files.append(infile)

offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000

with open(sys.argv[2]+'_bi_'+sys.argv[3]+'_bctdc.dat', 'w') as output_file:
	for i in range(0, len(files)):
		with open(files[i]) as csvfile:
			reader = csv.DictReader(csvfile)
			for row in reader:
				number_of_meas = int(row['nbOfMeas'])
				ts_exponent = int(float(row['measStamp_unitExponent']))
				intensity_exp = int(float(row['totalIntensity_unitExponent']))
				timestamp = int(row['acqTimestamp'])
				if number_of_meas != 0 and (row['SPS_USER'] == 'SPS.USER.LHCMD1' or row['SPS_USER'] == 'SPS.USER.LHCMD2'):
					len_array = int(number_of_meas/200) * 200
					list_timestamp = numpy.array(list(map(float, row['measStamp'].split('|'))))*10**(ts_exponent + 9)
					intensity = numpy.array(list(map(float, row['totalIntensity'].split('|'))))*10**intensity_exp
					intensity_mean = numpy.mean(numpy.reshape(numpy.resize(intensity, len_array), (-1, 200)), axis = 1)
					intensity_std = numpy.std(numpy.reshape(numpy.resize(intensity, len_array), (-1, 200)), axis = 1)
					time = numpy.reshape(numpy.resize(list_timestamp, len_array), (-1, 200))
					for l in range(0, len_array/200):
						output_file.write(str(int(((time[l][-1] + timestamp)*1E-6 - offset)*1E-3))+' '+str(intensity_mean[l])+' '+str(intensity_std[l])+'\n')
				else:
					for l in range(0, int(number_of_meas/200)):
						output_file.write(str(int(((timestamp + l)*1E-6 - offset)*1E-3))+' '+'0.0'+' '+'0.0'+'\n')
