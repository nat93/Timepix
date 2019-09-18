#!/usr/bin/env python

#########################################
#										#
# Stefano Petrucci - s.petrucci@cern.ch	#
#										#
#########################################

import ROOT 
import sys 
import os

def parse_CSV_file_with_TTree_ReadStream(tree_name, afile):
	ROOT.gROOT.SetBatch()
#	# The mapping dictionary defines the proper branch names and types given a header name.  
#	header_mapping_dictionary = { 
#		'timestamps'	: ('t'		, long),
#		'bct4'			: ('bct4'	, float),
#	}   
#
#	type_mapping_dictionary = { 
#		str		: 'C',
#		long	: 'L',
#		int		: 'I',
#		float	: 'F'
#	}
#
#	header_row = open(afile).readline().strip().split(' ')
#	print header_row
#	print ':'.join([header_mapping_dictionary[row][0]+'/'+type_mapping_dictionary[header_mapping_dictionary[row][1]] for row in header_row])
#
#	branch_descriptor = ':'.join([header_mapping_dictionary[row][0]+'/'+type_mapping_dictionary[header_mapping_dictionary[row][1]] for row in header_row])
#	print branch_descriptor
	
	output_ROOT_file_name = os.path.splitext(afile)[0] + '.root'
	output_file = ROOT.TFile(output_ROOT_file_name, 'recreate')
	print "Outputting %s -> %s" % (afile, output_ROOT_file_name)

	output_tree = ROOT.TTree(tree_name, tree_name)
	print afile
	output_tree.ReadFile(afile)

	output_file.cd()
	output_tree.Write()

######################### END OF FUNCTION #########################

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print "Usage: %s file_to_parse.dat" % sys.argv[0]
		sys.exit(1)
	parse_CSV_file_with_TTree_ReadStream("example_tree1", sys.argv[1])
