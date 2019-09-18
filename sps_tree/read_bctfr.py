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

#################################################################
#																#
# Script for reading the intensity of the most pupulated bunch.	#
# The measurements are averaged out on 60 seconds. This is		#
# possible averaging on 20 measurements (the measurements rate	#
# is 1/50 ms). For doing that the array with timestamps and		#
# with the measurements is reshaped (-1,20). For the average of	#
# the measurements is performed an average on the entire array.	#
# For the timestamps, is taken the last element of the			#
# 20-values array.												#
#																#
#################################################################

import csv
import glob
import re
import time
import datetime
import sys
import numpy
import math

# enlarging the file size of the single line otherwise sometime is impossible to read files
csv.field_size_limit(sys.maxsize)

numbers = re.compile(r'(\d+)')

# function for sort the list of file inside a folder in the right way
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

files = []
number_of_meas = 0
bunch_exp = 0
total_intensity_exp = 0
timestamp = 0				# to change the name: it is cycle stamp
array_len = 0				# Variable used for resizing the array if it is not multiple of "sample"
fill_delay = 1015			# delay in ms between the injection time and the start of the measStamp (for LHCMD1 and LHCMD2)
array_3bunch = []			# NOT USED?
user = 'SPS.USER.LHCMD1'

for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/BI instrumentation/BI_'+sys.argv[3]+'_BCTFR_*.csv'), key=numericalSort):
    files.append(infile)

offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000000000

sample = 20

if int(sys.argv[3]) == 31450:
	sample = 20
if int(sys.argv[3]) == 51895:
	sample = 20

for i in range(0, len(files)):
	with open(files[i]) as csvfile:
		print files[i]
		reader = csv.DictReader(csvfile)
		for row in reader:
# Reading the parameters and setting the arrays to zero
			number_of_meas = int(row['nbOfMeas'])
			bunch_exp = int(row['gateIntensity_unitExponent'])
			total_intensity_exp = int(row['totalIntensity_unitExponent'])
# TIMESTAMP IN ns
			timestamp = int(row['cycleStamp']) + (fill_delay * 1000000)
# Stamps in ms
			stamps = numpy.zeros(number_of_meas)
# Selecting only data from the right SPS cycle
			if row['cycleName'] == user and int(row['nbOfMeas']) != 0:
#			if (row['cycleName'] == 'SPS.USER.LHCMD1' or row['cycleName'] == 'SPS.USER.LHCMD2') and (int(row['nbOfMeas']) != 0):
#			if (row['cycleName'] == 'SPS.USER.LHCION3') and int(row['nbOfMeas']) != 0 and int(row['cycleStamp']) != 0:
				print row['cycleName'], 'measStamp = ', row['measStamp'], len(row['measStamp'])
				array_len = (int(number_of_meas) / sample) * sample
				stamps = numpy.array(list(map(int, row['measStamp'].split('|')))) * 1000000 + timestamp
				bunch_intensity_fast = numpy.array(list(map(float, row['bunchIntensityFast'].split('|'))))*10**bunch_exp
				total_intensity_fast = numpy.array(list(map(float, row['totalIntensityFast'].split('|'))))*10**total_intensity_exp
# Creating the matrix with the intensity of every bunch
				matrix_intensity_fast = numpy.reshape(bunch_intensity_fast, (number_of_meas, -1))
# Selectin the array with the most populated bunch, resizing it with array_len (deleting the measurements in excess),
# reshaping for averaging on "sample" elements, averaging the new array
				mean_array = numpy.mean(numpy.reshape(numpy.resize(numpy.amax(matrix_intensity_fast, axis = 1), array_len), (-1, sample)), axis = 1)
				std_array = numpy.std(numpy.reshape(numpy.resize(numpy.amax(matrix_intensity_fast, axis = 1), array_len), (-1, sample)), axis = 1) / math.sqrt(sample)
				mean_intensity_fast = numpy.mean(numpy.reshape(numpy.resize(total_intensity_fast, array_len), (-1, sample)), axis = 1)
#				for jj in range (0, len(numpy.amax(matrix_intensity_fast, axis = 1))):
#					if numpy.argmax(matrix_intensity_fast, axis = 1)[jj] + 2 < len(numpy.amax(matrix_intensity_fast, axis = 0)):
#						array_3bunch[jj] = matrix_intensity_fast[jj][numpy.argmax(matrix_intensity_fast, axis = 1)[jj] + 2] + matrix_intensity_fast[jj][numpy.argmax(matrix_intensity_fast, axis = 1)[jj] + 1] + matrix_intensity_fast[jj][numpy.argmax(matrix_intensity_fast, axis = 1)[jj]]
#						print 'ok'
#					if numpy.argmax(matrix_intensity_fast, axis = 1)[jj] + 1 < len(numpy.amax(matrix_intensity_fast, axis = 0)):
#						array_3bunch[jj] = matrix_intensity_fast[jj][numpy.argmax(matrix_intensity_fast, axis = 1)[jj] + 1] + matrix_intensity_fast[jj][numpy.argmax(matrix_intensity_fast, axis = 1)[jj]]
#						print 'ok'
#					if numpy.argmax(matrix_intensity_fast, axis = 1)[jj] < len(numpy.amax(matrix_intensity_fast, axis = 0)):
#						array_3bunch[jj] =  matrix_intensity_fast[jj][numpy.argmax(matrix_intensity_fast, axis = 1)[jj]]
#						print 'ok'

				for l in range(0, array_len/sample):
					print str(((int(numpy.reshape(numpy.resize(stamps, array_len), (-1, sample))[l][-1]) - offset)/1000000000)), mean_array[l], std_array[l], mean_intensity_fast[l]
# The other measurements/SPS cycle
			elif int(row['nbOfMeas']) != 0:
				bunch_intensity_fast = numpy.zeros(int(int(row['measStamp'].split('|')[-1])/1000))
				for k in range(0, int(int(row['measStamp'].split('|')[-1])/1000)):
					stamps[k] = ((timestamp + (int(k) * 1000000000) - offset )/1000000000)
				for time, meas in zip(stamps, bunch_intensity_fast):
					print str(int(time)), meas, 0.0, 0.0
