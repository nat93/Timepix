#!/bin/bash

##---------------------------------------------------------##
##
##---------------------------------------------------------##
data_dir=/media/andrii/F492773C92770302/MedipixData/SPS_DATA/MD_2018_09_17;
#data_dir=/media/andrii/F492773C92770302/MedipixData/H8_DATA/2018_09_12_pions;
#data_dir=/home/anatochi/Medipix/H8_DATA/2018_09_12_pions;
#data_dir=/home/anatochi/Medipix/SPS_DATA/MD_2018_06_18;
#data_dir=/home/anatochi/Medipix/SPS_DATA/MD_2018_08_15;
#data_dir=/media/andrii/F492773C92770302/MedipixData/SPS_DATA/MD_2018_10_24;

for runrunID in $(seq 8 1 8)
do

    #1

    ## ASCII FORMAT ##

    cd $data_dir/RUN_$runrunID;
    nFiles=$(ls -1q *.dat | wc -l);

    echo "";
    echo "Number of files in the directory: "$nFiles;
    echo "";

    cd -;

#    make clean; make ascii2root_common;
#    for runid in $(seq 1 1 $nFiles)
#    do
#        if [ -e $data_dir/RUN_$runrunID/Medipix_$runid.dat ]
#        then
#            ./ascii2root_common $data_dir/RUN_$runrunID/Medipix_$runid.dat $data_dir/RUN_$runrunID/Medipix_$runid.root
#        else
#            echo "WARNING!!! The $data_dir/RUN_$runrunID/Medipix_$runid.dat cannot be find"
#        fi
#    done

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

    make convert_common;
    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 530 559 /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L19_RUN_$runrunID.root
#    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 3044 3094 /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_LCP61_CH_RUN.root
#    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 323 330 /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_STF123_AM_RUN.root

#    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_RUN_$runrunID.root
#    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 345 360 /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_M1_RUN_$runrunID.root
#    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 362 363 /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_M2_RUN_$runrunID.root
#    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 250 285 /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_M3_RUN_$runrunID.root
#    ./convert_common $data_dir/RUN_$runrunID/Medipix_ 1 $nFiles /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_10_24_RUN_$runrunID.root

    #3

#    make analysis_common;
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_10_24_RUN_$runrunID.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_10_24_HISTO_RP0E_RUN_$runrunID.root 0
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_10_24_RUN_$runrunID.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_10_24_HISTO_RP1I_RUN_$runrunID.root 1
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_10_24_RUN_$runrunID.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_10_24_HISTO_RP3E_RUN_$runrunID.root 2

#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_STF123_AL_RUN.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_STF123_AL_RUN_CHIP_0.root 0
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_STF123_AL_RUN.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_STF123_AL_RUN_CHIP_1.root 1
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_LCP61_CH_RUN.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_LCP61_CH_RUN_CHIP_0.root 0
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_LCP61_CH_RUN.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_LCP61_CH_RUN_CHIP_1.root 1
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_ACP85_AM_RUN.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_ACP85_AM_RUN_CHIP_0.root 0
#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_ACP85_AM_RUN.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/H8_2018_09_12_ACP85_AM_RUN_CHIP_1.root 1

#    ./analysis_common /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L2_RUN_$runrunID.root /media/andrii/F492773C92770302/MedipixData/ROOT_FILES/MD_2018_09_17_L2_HISTO_RP0I_RUN_$runrunID.root 3
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_09_17_K0_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_09_17_K0_HISTO_RP0I_RUN_$runrunID.root 3
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip0_RUN_$runrunID.root 0
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_09_12_pions_HISTO_Chip1_RUN_$runrunID.root 1
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_M1_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_HISTO_RP1I_M1_RUN_$runrunID.root 3
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_M1_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_HISTO_RP0I_M1_RUN_$runrunID.root 4
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_M2_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_HISTO_RP1I_M2_RUN_$runrunID.root 3
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_M2_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_HISTO_RP0I_M2_RUN_$runrunID.root 4
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_M3_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP1I_M3_RUN_$runrunID.root 3
#    ./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_M3_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_08_15_HISTO_RP0I_M3_RUN_$runrunID.root 4

    #4

    #make clean; make analysis_trk;
    #./analysis_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_CLUSTERINFO_RUN_$runrunID.root

    #5

    #make clean; make trackreco_trk;
    #./trackreco_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_CLUSTERINFO_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_TRACKINFO_RUN_MODE_2_RUN_$runrunID.root 2

done

# From SPS:
#   devRP0E=0(0) #G02-W0108    FITpix 0384
#   devRP3I=X(1) #K09-W0255    FITpix XXXX
#   devRP3E=2(2) #C08-W0255    FITpix 0399
#   devRP1I=1(3) #F04-W0108    FITpix 0393
#   devRP0I=3(4) #I02-W0108    FITpix 0415


# git commit -am "third commit"
# git push origin develop

