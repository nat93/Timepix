#!/bin/bash

##----------------------##
## SPS MD auto analysis ##
##----------------------##

# Good channeling: Wed Oct 18 2017 01:23:00 | Final time: Wed Oct 18 2017 01:28:00

data_dir_in=/home/andrii/Medipix/SPS_DATA/MD_2017_10_17/TIMEPIX/; # Input directory name
data_dir_out=/home/andrii/Medipix/ROOT_FILES/; # Output directory name
data_file_out_1=SPS_MD_2017_10_17.root; # Output file name with data
data_file_out_2=SPS_MD_2017_10_17_HISTO.root; # Output file name with histo
data_file_out_3=SPS_MD_2017_10_17_HISTO_old.root; # Output file name with histo
ut_min=1508305080000; # Minimum UnixTime for data selection --> Tuesday, October 18, 2017 07:38:00 GMT+02:00 DST
ut_max=1508306400000; # Maximum UnixTime for data selection --> Tuesday, October 18, 2017 08:00:00 GMT+02:00 DST
start_event=630643;
chipid=2;

make convert_v5_2;
make analysis_v5_2;

nfiles_old=$(ls -1 $data_dir_in | wc -l); # Initial number of files in the directory
ls -v $data_dir_in*.root > filelist2.dat; # Initial list of the files in the directory

time ./convert_v5_2 filelist2.dat  $ut_min $ut_max $data_dir_out$data_file_out_1 $start_event; # Convert existing files in the directory
time ./analysis_v5_2 $data_dir_out$data_file_out_1 $data_dir_out$data_file_out_2 $chipid $ut_min $ut_max # Analysis with plotting histograms
cp $data_dir_out$data_file_out_2 $data_dir_out$data_file_out_3;

while true; do
    nfiles_new=$(ls -1 $data_dir_in | wc -l);
    if [ "$nfiles_new" -ne "$nfiles_old" ]
    then
        nfiles_old=$nfiles_new;
        ls -v $data_dir_in*.root > filelist2.dat;
        time ./convert_v5_2 filelist2.dat  $ut_min $ut_max $data_dir_out$data_file_out_1 $start_event;
        time ./analysis_v5_2 $data_dir_out$data_file_out_1 $data_dir_out$data_file_out_2 $chipid $ut_min $ut_max
        cp $data_dir_out$data_file_out_2 $data_dir_out$data_file_out_3;
    else
        echo "... Wating for the new file ...";
        sleep 10 #Wait 30 seconds
    fi
done




