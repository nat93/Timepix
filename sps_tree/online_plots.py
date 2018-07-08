import glob
import csv
import re
import time
import datetime
import sys

numbers = re.compile(r'(\d+)')

def numericalSort(value):
	parts = numbers.split(value)
	parts[1::2] = map(int, parts[1::2])
	return parts

files = []
measures = []
timestamp = 0
avg_measures = 0
offset = int(time.mktime(datetime.datetime.strptime('2016_10_18', '%Y_%m_%d').timetuple()))
print offset
mylist = []

for infile in sorted(glob.glob('/media/dfs/Experiments/UA9/Data_SPS_runs/MD_2016_10_18/Detectors Data/Scintillators_*.csv'), key=numericalSort):
    files.append(infile)
for i in range(0, len(files)):
	with open(files[i]) as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			measures = row['SCI_F'].split('|')
			timestamp = long(row['aqnStamp'].split('|')[0])
			mylist.append(timestamp)
			for meas in measures:
				avg_measures += int(meas)
			print long(timestamp*1E-9 - offset), int(avg_measures) / len(measures)
			avg_measures = 0

#temp = long(0)
#j = 0
##kk = 0
##mylist.sort()
##diff = mylist[150] - mylist[149]
##
##for k in range(len(mylist)-1):
##	print mylist[k]
###	if mylist[k+1] - mylist[k] > diff + 274951625:
###		kk += 1
####		print mylist[k+1] - mylist[k]
###print kk, len(mylist)
#
#for element in mylist:
#	if long(temp) > long(element):
#		j += 1
#		print long(temp), long(element)
#	else: 
#		temp = long(element)
#print j, len(mylist), str(float((j * 100) / len(mylist)))+"%"
