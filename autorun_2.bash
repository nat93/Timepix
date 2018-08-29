#!/bin/bash

##----------------------##
## SPS MD auto analysis ##
##----------------------##

data_dir=/home/anatochi/Medipix/SPS_DATA/TEST/TIMEPIX;

#1

cd $data_dir;
nFiles=$(ls -1q *.root | wc -l);

echo "";
echo "Number of files in the directory: "$nFiles;
echo "";

cd -;

#2

make convert_common;
./convert_common $data_dir/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/TEST.root

#3

make analysis_common;
./analysis_common /home/anatochi/Medipix/ROOT_FILES/TEST.root /home/anatochi/Medipix/ROOT_FILES/TEST_HISTO.root 3

# devRP0E=0 #G02-W0108    FITpix 0384
# devRP3I=1 #K09-W0255    FITpix 0393
# devRP3E=2 #C08-W0255    FITpix 0399
# devRP1I=3 #F04-W0108    FITpix 0409
# devRP0I=4 #I02-W0108    FITpix 0415



