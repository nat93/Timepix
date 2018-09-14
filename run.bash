#!/bin/bash

##---------------------------------------------------------##
##
##---------------------------------------------------------##
#data_dir=/home/anatochi/Medipix/H8_DATA/2018_04_11_pions;
#data_dir=/home/anatochi/Medipix/SPS_DATA/TEST_2018_08_07;
#data_dir=/home/anatochi/Medipix/SPS_DATA/MD_2018_08_15;
#data_dir=/home/anatochi/Medipix/SPS_DATA/MD_2017_10_17/TIMEPIX;
#data_dir=/home/anatochi/Medipix/SPS_DATA/MD_2018_06_18;
data_dir=/home/anatochi/Medipix/H8_DATA/2018_09_12_pions;

for runrunID in $(seq 1 1 26)
do

    #1

    ## ASCII FORMAT ##

    cd $data_dir/RUN_$runrunID;
    nFiles=$(ls -1q *.dat | wc -l);

    echo "";
    echo "Number of files in the directory: "$nFiles;
    echo "";

    cd -;

    make clean; make ascii2root_common;
    for runid in $(seq 1 1 $nFiles)
    do
            ./ascii2root_common $data_dir/RUN_$runrunID/Medipix_$runid.dat $data_dir/RUN_$runrunID/Medipix_$runid.root
    done

    ## CSV FORMAT ##

    #cd $data_dir/RUN_$runrunID;
    #nFiles=$(ls -1q *.csv | wc -l);

    #echo "";
    #echo "Number of files in the directory: "$nFiles;
    #echo "";

    #cd -;
    #let "nFiles -= 1";

    #make clean; make csv2root_common;
    #for runid in $(seq 0 1 $nFiles)
    #do
    #        echo "";
    #        du -sh $data_dir/RUN_$runrunID/Medipix_$runid.csv;
    #        ./csv2root_common $data_dir/RUN_$runrunID/Medipix_$runid.csv $data_dir/RUN_$runrunID/Medipix_$runid.root
    #done

    #2

    make clean; make convert_common;
    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_RUN_$runrunID.root
    #./convert_common $data_dir/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_RUN_$runrunID.root
    #./convert_common $data_dir/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/TEST_2018_08_07_RUN_$runrunID.root
    #./convert_common $data_dir/RUN_$runrunID/Medipix_ 507 508 /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_RUN_$runrunID.root
    #./convert_common $data_dir/RUN_$runrunID/Medipix_ 232 261 /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_TACW_SCAN_RUN_$runrunID.root
    #./convert_common $data_dir/RUN_$runrunID/Medipix_ 970 973 /home/anatochi/Medipix/ROOT_FILES/MD_2017_10_17_CRY2CHCRY4AMP1_RUN_$runrunID.root

    #3

    make clean; make analysis_common;
    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_$runrunID.root 0
    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_$runrunID.root 1
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_TACW_SCAN_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_TACW_SCAN_HISTO_RP1I_RUN_$runrunID.root 3
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_TACW_SCAN_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_TACW_SCAN_HISTO_RP0I_RUN_$runrunID.root 4
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2017_10_17_CRY2CH_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2017_10_17_CRY2CH_HISTO_RP1I_RUN_$runrunID.root 2
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2017_10_17_CRY2CHCRY4CHP7_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2017_10_17_CRY2CHCRY4CHP7_HISTO_RP1I_RUN_$runrunID.root 2
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2017_10_17_CRY2CHCRY4AMP1_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2017_10_17_CRY2CHCRY4AMP1_HISTO_RP1I_RUN_$runrunID.root 2

    #4

    #make clean; make analysis_trk;
    #./analysis_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_CLUSTERINFO_RUN_$runrunID.root

    #5

    #make clean; make trackreco_trk;
    #./trackreco_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_CLUSTERINFO_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_TRACKINFO_RUN_MODE_2_RUN_$runrunID.root 2

done

# From SPS:
#   devRP0E=0 #G02-W0108    FITpix 0384
#   devRP3I=1 #K09-W0255    FITpix 0393
#   devRP3E=2 #C08-W0255    FITpix 0399
#   devRP1I=3 #F04-W0108    FITpix 0409
#   devRP0I=4 #I02-W0108    FITpix 0415


# git commit -am "third commit"
# git push origin develop

