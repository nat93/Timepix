#################################################
#												#
# Simone Montesano - simone.montesano@cern.ch	#
#												#
# Stefano Petrucci - s.petrucci@cern.ch			#
#												#
#################################################

import csv
import numpy
import glob
import re
import time
import datetime
import sys

numbers = re.compile(r'(\d+)')

def numericalSort(value):
	parts = numbers.split(value)
	parts[1::2] = map(int, parts[1::2])
	return parts

counternames = []
timestamps = []
countervalues = {}
counts = {}
folder = []
offset = int(time.mktime(datetime.datetime.strptime(sys.argv[1], '%Y_%m_%d').timetuple()))*1000

for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/Detectors Data/Wavecatcher_*.csv'), key=numericalSort):
	folder.append(infile)

counternames = ['WC_A+RF', 'WC_AA+RF', 'WC_AD+I', 'WC_AD+I+RF', 'WC_AD+RF', 'WC_AH+AI', 'WC_AH+AI+RF', 'WC_AH+RF', 'WC_AI+RF', 'WC_C+D', 'WC_C+D+RF', 'WC_C+RF', 'WC_D+RF', 'WC_E+N', 'WC_E+N+RF', 'WC_E+RF', 'WC_G+M', 'WC_G+M+RF', 'WC_G+RF', 'WC_H+RF', 'WC_I+RF', 'WC_J+K', 'WC_J+K+RF', 'WC_J+RF', 'WC_K+RF', 'WC_M+RF', 'WC_N+RF', 'WC_RF', 'WC_Z+H+RF', 'WC_Z+RF']

for items in counternames:
	countervalues[items] = []

for i in range(0, len(folder)):
	with open(folder[i]) as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			time_vector = list(map(int, row['aqnStamp'].split('|')))
			for l in range(0, len(time_vector)):
				time_vector[l] = str(int(int(time_vector[l])/10**6 - int(offset))/1000)
			for items in counternames:
				counts[items] = list(map(int, row[items].split('|')))
				countervalues[items] = sum(counts[items])
			for items in counternames:
				out_file = open('/home/petruccs-local/scripts/'+sys.argv[1]+'/'+sys.argv[1]+'_'+items+"_counters"+".dat", "a")
				out_file.write(time_vector[-1]+' '+str(countervalues[items])+'\n')
				out_file.close()
			countervalues[items] = []
