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
folder = []
offset = int(time.mktime(datetime.datetime.strptime(sys.argv[2], '%Y_%m_%d').timetuple()))*1000

for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_'+sys.argv[1]+'/Detectors Data/ScintillatorsCounters_*.csv'), key=numericalSort):
	folder.append(infile)

# fare un file per ogni canale e levare canali n.c. frev_N con N > 1, etc

#j = 0
for i in range(0, len(folder)):
	with open(folder[i]) as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			if counternames == []:
				counternames = row['counterNames'].split('|')
				counternames = filter(None, counternames)
				for name in counternames:
					countervalues[name] = []
			time_vector = list(map(int, row['timestamps'].split('|')))
			for l in range(0, len(time_vector)):
				time_vector[l] = str(int(int(time_vector[l])/10**6 - int(offset))/1000)
			num_values = len(time_vector);
			num_counters = len(counternames)
			timestamps.extend(time_vector)
			value_vector = numpy.array(list(map(float, row['counterValues'].split('|'))))
			if num_values * num_counters !=len(value_vector):
				print("ERROR!! Dimension of counters*values matrix are wrong!")
			value_matrix = numpy.reshape(value_vector,(num_counters, num_values))
			for name, values in zip(counternames, value_matrix):
#				countervalues[name].extend(values)
				countervalues[name].append(sum(values))
			for items in counternames:
				if items != 'n.c.' and items != 'n.c':
					out_file = open('/home/petruccs-local/scripts/'+sys.argv[2]+'/'+sys.argv[2]+'_'+items+"_counters"+".dat", "a")
#					if j == 0 and i == 0:
#						out_file.write("#timestamps,counterValues\n")
#						j = 1
					out_file.write('\n'.join('%s %s' % x for x in zip(time_vector, countervalues[items])))
					out_file.write('\n')
					out_file.close()
#				j = 0
				countervalues[items] = []
