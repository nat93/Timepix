#########################################
#                                       #
# Stefano Petrucci - s.petrucci@cern.ch #
#                                       #
#########################################

import sys

data = open(sys.argv[1], 'r')
measurements = 0
j = 1
avg = 200

for lines in data:
	measurements += float(lines.split(' ')[1])
	j += 1
	if j == avg:
		print str(int(float(lines.split(' ')[0])/1000))+' '+str(measurements/avg)
		j = 0
		measurements = 0
