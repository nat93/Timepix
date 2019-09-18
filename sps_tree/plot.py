# Stefano Petrucci (s.petrucci@cern.ch
#
# Script for plotting the wavecatcher .dat files

import matplotlib.pyplot as plt
import argparse
import glob
import time
import sys
from sys import stdout
from time import sleep

# Header file
# ChannelNb	ChannelName	EventID	MeasuredBaseline[Volts]	Amplitude[Volts]	PeakPosition[Cell]	Charge[pC]	LeadingEdgeTime[ns]	TrailingEdgeTime[ns]
# 0			1			2		3						4					5					6			7					8						9
# ChannelNb	ChannelName	EventID	MeasuredBaseline[Volts]	Amplitude[Volts]	PeakTime[ns]		Charge[pC]	LeadingEdgeTime[ns]	TrailingEdgeTime[ns]	Peaks[count]
#
# Special operation value (special_op variable)
# Trailing edge time - leading edge time	=> 1
# Trailing edge time - frev					=> 2 to be implemented



#################### START FUNCTION DECLARATION ##################

def get_channel_name(fileName):
	g = open(fileName, "r")
	count = 0
	dictionary = {}
	for line in g:
		if count == 4:
			break
		if line.startswith("==="):
			count += 1
		if not line.startswith("==="):
			dictionary[line.split(" ")[1]] = line.split(" ")[0]
	return dictionary
	g.close()

#################### END FUNCTION DECLARATION ####################

############### PARSING THE COMMAND LINE VARIABLES ###############

parser = argparse.ArgumentParser()
#parser.add_argument('string', help="first string is the folder's path, second is the measure to plot", nargs=1)
parser.add_argument('-p', help='full path to the folder to analyze', action='store', dest='path', type=str, required=True)
parser.add_argument('-a', help='action to execute i.e. amplitude or leading_edge_time, etc.', action='store', dest='action', type=str, required=True)
parser.add_argument('-c', help='channel name, it is case sensitive and for now there is no protection against wrong parameter', action='store', dest='channel_name', type=str, required=True)
parser.add_argument('-s', help='save the histogram to a pdf file', action='store_true', dest='save_to_file', default=False)
parser.add_argument('-t', help='threshold value in mV', action='store', dest='threshold', type=float, required=False)
args = parser.parse_args()

############# END PARSING THE COMMAND LINE VARIABLES #############

folder = glob.glob(args.path+"measurement*.dat")

op = -1
special_op = -1
threshold = -0.020 # Volts

if args.threshold:
	threshold = args.threshold

channel_dict = get_channel_name(folder[0])
progress = 0

data = []
f_rev = 0.
f = [0 for x in range(0, len(folder))]
total_entries = 0

# parsing of the operations
if args.action == "base_line":
	op = 3
if args.action == "amplitude":
	op = 4
if args.action == "peak_time":
	op = 5
if args.action == "charge":
	op = 6 
if args.action == "leading_edge_time":
	op = 7
if args.action == "trailing_edge_time":
	op = 8
if args.action == "leading_trailing_diff":
	special_op = 1
if args.action == "trailing_no_frev":
	special_op = 2 
# configuration for other parameters
# time-window for time plots
min_time = -3000
max_time = 3000

# Checking op and special_op value (only if there is no operation/special operation)
if op == -1 & special_op == -1:
	sys.exit("ERROR: wrong action")

print("Progress:")
for i in range(0,len(folder)):
	f[i] = open(folder[i], "r")
	progress = i*100/(len(folder)-1)
	stdout.write("\r%d" %progress)
	stdout.flush()
	for line in f[i]:
		if line.startswith("==="):
			continue
		if line.startswith("2 frev1"):
			f_rev = float(line.split(" ")[8])
		if line.startswith(channel_dict[args.channel_name]+" "+args.channel_name):
			total_entries = total_entries + 1
			if total_entries > 10000:
				break
			if special_op < 0:
				if float(line.split(" ")[4]) < threshold:
					data += ( float(line.split(" ")[op])),
			if special_op == 1:
				if float(line.split(" ")[4]) < threshold:
					# remove this if if you don't need trailing - leading
					if not(float(line.split(" ")[7]) == 0.0 and float(line.split(" ")[8])!= 0.0):
						data += ( float( line.split(" ")[7]) - float(line.split(" ")[8] ) ),
			if special_op == 2:
				if float(line.split(" ")[4]) < threshold:
					data += ( float( line.split(" ")[8]) - f_rev ),
	f[i].close()
print("\nNumber of entries: "+str(len(data)))

# HISTOGRAM
#plt.hist(data, bins = 1024, range=[0, 3000])
plt.hist(data, bins = 1024)
plt.title(args.action+" "+" - CH "+args.channel_name+" - thr = "+str(threshold)+" V")
plt.yscale('log', nonposy='clip')
plt.ylabel('Counts')
if op == 4 or op == 3: 
	plt.xlabel('Volts')
if op == 6: 
	plt.xlabel('pC')
if op == 7 or op == 8 or op == 5 or special_op > 0:
	plt.xlabel('Time (ns)')

if args.save_to_file == True:
	plt.savefig(args.action+"_"+args.channel_name+".pdf", format='pdf')
	print("Saved histogram in the current directory as "+args.action+"_"+args.channel_name+".pdf")
else:
	plt.show()
