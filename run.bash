#!/bin/bash

#make clean;
#make;

#./csv2root ./medipix_h8_data/Medipix_465.dat ./medipix_h8_data/Medipix_465.root
#./ascii2root ./medipix_h8_data/Medipix_1.dat ./medipix_h8_data/Medipix_1.root
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip87_1ms.dat ./medipix_h8_data/Timepix_Chip87_1ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip87_2ms.dat ./medipix_h8_data/Timepix_Chip87_2ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip87_5ms.dat ./medipix_h8_data/Timepix_Chip87_5ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip87_10ms.dat ./medipix_h8_data/Timepix_Chip87_10ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip87_20ms.dat ./medipix_h8_data/Timepix_Chip87_20ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip87_100ms.dat ./medipix_h8_data/Timepix_Chip87_100ms.root;

#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip89_1ms.dat ./medipix_h8_data/Timepix_Chip89_1ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip89_2ms.dat ./medipix_h8_data/Timepix_Chip89_2ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip89_5ms.dat ./medipix_h8_data/Timepix_Chip89_5ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip89_10ms.dat ./medipix_h8_data/Timepix_Chip89_10ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip89_20ms.dat ./medipix_h8_data/Timepix_Chip89_20ms.root;
#./ascii2root_v2 ./medipix_h8_data/Timepix_Chip89_100ms.dat ./medipix_h8_data/Timepix_Chip89_100ms.root;

#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6001.dat ./medipix_h8_data/tracks_rec/Medipix_6001.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6002.dat ./medipix_h8_data/tracks_rec/Medipix_6002.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6003.dat ./medipix_h8_data/tracks_rec/Medipix_6003.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6004.dat ./medipix_h8_data/tracks_rec/Medipix_6004.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6005.dat ./medipix_h8_data/tracks_rec/Medipix_6005.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6006.dat ./medipix_h8_data/tracks_rec/Medipix_6006.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6007.dat ./medipix_h8_data/tracks_rec/Medipix_6007.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6008.dat ./medipix_h8_data/tracks_rec/Medipix_6008.root
#./ascii2root ./medipix_h8_data/tracks_rec/Medipix_6009.dat ./medipix_h8_data/tracks_rec/Medipix_6009.root

#make clean;
#make plot_quad;
#./plot ./medipix_h8_data/Medipix_7000.dat OUTPUT_FRAME.root -1 0 255 0 255
#./plot_quad ./medipix_h8_data/Medipix_Deflected.dat OUTPUT_FRAME_Deflected.root -1 175 400 200 400
#./plot_quad ./medipix_h8_data/Medipix_NotDeflected.dat OUTPUT_FRAME_NotDeflected.root -1 150 500 200 400

#while true; do
#        echo ''
#        date
#        echo ''
#        ./plot ./medipix_h8_data/Medipix_7000.dat OUTPUT_FRAME.root -1 0 255 0 255
#        echo ''
#        echo 'Waiting 3 seconds ...'
#        sleep 3;
#done

#for runid in $(seq 0 1 537)
#do
#        ./csv2root /home/anatochi/dfs/Experiments/UA9/Data_SPS_runs/MD_2017_07_04/Detectors\ Data/Medipix_$runid.csv ./sps_data/2017_07_04/Medipix_$runid.root
#done

#make clean;
#make analysis2;
#./analysis2 OUTPUT_COMMON.root ./sps_data/aligned_data.root OUTPUT.root

#for runid in $(seq 1 1 5)
#do
#        ./ascii2root_v3 /home/anatochi/Medipix/TEST_2017_07_11/3MHz/Medipix_$runid.dat /home/anatochi/Medipix/TEST_2017_07_11/3MHz/Medipix_$runid.root
#done

#./convert OUTPUT_MD_EXTRACTION.root;
#./convert_v3 OUTPUT_MD_EXTRACTION_LOCAL_PC.root;

##----------------------##
## Test-beam H8 TOA TMP ##
##----------------------##

#1
#runrunID=30
#nFiles=115

#make clean; make ascii2root_v5;
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_1Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_1Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_5Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_5Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_10Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_10Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_15Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_15Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_20Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_20Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_25Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_25Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_50Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_50Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_100Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_100Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_1000Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_1000Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_10000Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_10000Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_100000Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_100000Hz.root
#./ascii2root_v5 ../H8_DATA/2017_12_04_ions/Medipix_1000000Hz.dat ../H8_DATA/2017_12_04_ions/Medipix_1000000Hz.root

#for runid in $(seq 86 1 $nFiles)
#do
        #./ascii2root_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_27_PL08_TOA/DAT/Medipix_$runid.dat /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_27_PL08_TOA/ROOT/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_28_PL08_MPX/DAT/Medipix_$runid.dat /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_28_PL08_MPX/ROOT/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_24_QMP54/DAT/Medipix_$runid.dat /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_24_QMP54/ROOT/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_QMP54/DAT/Medipix_$runid.dat /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_QMP54/ROOT/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_26_QMP54/DAT/Medipix_$runid.dat /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_26_QMP54/ROOT/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/Medipix_$runid.dat /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/Medipix_$runid.root
        #./ascii2root_v5 ./sps_data/SPS_MD_2017_09_18_LOCAL_PC/Medipix_$runid.dat ./sps_data/SPS_MD_2017_09_18_LOCAL_PC/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/SPS_DATA/MD_2017_09_18/TIMEPIX/Medipix_$runid.dat /home/anatochi/Medipix/SPS_DATA/MD_2017_09_18/TIMEPIX/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/SPS_DATA/MD_2017_09_27/TIMEPIX/TEST_1/Medipix_$runid.dat /home/anatochi/Medipix/SPS_DATA/MD_2017_09_27/TIMEPIX/TEST_1/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/SPS_DATA/MD_2017_09_27/TIMEPIX/Medipix_$runid.dat /home/anatochi/Medipix/SPS_DATA/MD_2017_09_27/TIMEPIX/Medipix_$runid.root
        #./ascii2root_v5 ../SPS_DATA/MD_2017_09_27/TIMEPIX/Medipix_$runid.dat ../SPS_DATA/MD_2017_09_27/TIMEPIX/Medipix_$runid.root
        #./ascii2root_v5 ../H8_DATA/2017_12_ions/Medipix_$runid.dat ../H8_DATA/2017_12_ions/Medipix_$runid.root

        #./ascii2root_v5 ../H8_DATA/2018_04_07_protons/difuser/Medipix_$runid.dat ../H8_DATA/2018_04_07_protons/difuser/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/H8_DATA/2018_04_11_pions/RUN_1/Medipix_$runid.dat /home/anatochi/Medipix/H8_DATA/2018_04_11_pions/RUN_1/Medipix_$runid.root
        #./ascii2root_v5 ../H8_DATA/2018_04_15_pions/RUN_$runrunID/Medipix_$runid.dat ../H8_DATA/2018_04_15_pions/RUN_$runrunID/Medipix_$runid.root
        #./ascii2root_v5 /home/anatochi/Medipix/SPS_DATA/MD_2018_06_18/Medipix_$runid.dat /home/anatochi/Medipix/SPS_DATA/MD_2018_06_18/Medipix_$runid.root
#done

#make clean; make ascii2root_v4;
#./ascii2root_v4 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile1MHz.dat /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile1MHz.root 1 11810
#./ascii2root_v4 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile3MHz.dat /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile3MHz.root 3 3936
#./ascii2root_v4 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile6MHz.dat /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile6MHz.root 6 1968
#./ascii2root_v4 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile12MHz.dat /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile12MHz.root 12 984
#./ascii2root_v4 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile24MHz.dat /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile24MHz.root 24 492
#./ascii2root_v4 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile48MHz.dat /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile48MHz.root 48 246
#./ascii2root_v4 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile96MHz.dat /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile96MHz.root 96 123


#2

#make clean; make convert_v5;
#./convert_v5 OUTPUT_H8_PIONS_PL08_TOA.root
#./convert_v5 OUTPUT_H8_PIONS_PL08_MPX.root
#./convert_v5 OUTPUT_H8_PIONS_QMP54_24_TOA.root
#./convert_v5 OUTPUT_H8_PIONS_QMP54_25_TOA.root
#./convert_v5 OUTPUT_H8_PIONS_QMP54_26_TOA.root
#./convert_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/TMP_CpFM_RUN2.root
#./convert_v5 ./sps_data/SPS_MD_2017_09_18_LOCAL_PC/SPS_MD_2017_09_18.root
#./convert_v5 /home/anatochi/Medipix/SPS_DATA/MD_2017_09_18/TIMEPIX/Medipix_ 1 628 /home/anatochi/Medipix/ROOT_FILES/SPS_MD_2017_09_18.root
#./convert_v5 /home/anatochi/Medipix/SPS_DATA/MD_2017_09_27/TIMEPIX/TEST_1/Medipix_ 1 8 /home/anatochi/Medipix/ROOT_FILES/SPS_MD_2017_09_27.root
#./convert_v5 /home/anatochi/Medipix/SPS_DATA/MD_2017_09_27/TIMEPIX/Medipix_ 90 133 /home/anatochi/Medipix/ROOT_FILES/SPS_MD_2017_09_27.root
#./convert_v5 ../SPS_DATA/MD_2017_09_18/TIMEPIX/Medipix_ 1 628 ../ROOT_FILES/SPS_MD_2017_09_18.root
#./convert_v5 ../SPS_DATA/MD_2017_10_17/TIMEPIX/Medipix_ 1 1320 ../ROOT_FILES/SPS_MD_2017_10_17.root
#./convert_v5 ../H8_DATA/2017_12_ions/Medipix_ 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_05.root
#./convert_v5 ../H8_DATA/2017_12_ions/Medipix_ 2 2 ../ROOT_FILES/H8_Test_Beam_2017_12_05_second.root
#./convert_v5 /home/anatochi/Medipix/SPS_DATA/MD_2017_09_27/TIMEPIX/Medipix_ 1 328 /home/anatochi/Medipix/ROOT_FILES/SPS_MD_2017_09_27.root

#./convert_v5 /home/anatochi/Medipix/H8_DATA/2018_04_07_protons/difuser/Medipix_ 1 3 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_07_protons_DIFUSER.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/2018_04_07_protons/blade/Medipix_ 1 3 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_07_protons_BLADE.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/2018_04_11_pions/RUN_1/Medipix_ 1 37 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_1.root
#./convert_v5 ../H8_DATA/2018_04_15_pions/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_$runrunID.root

#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_1Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_5Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_10Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_15Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_15Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_20Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_25Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_25Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_50Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_100Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_100Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_1000Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_10000Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10000Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_100000Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_100000Hz.root
#./convert_v5 ../H8_DATA/2017_12_04_ions/Medipix_1000000Hz 1 1 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz.root

#./convert_v5 /home/anatochi/Medipix/H8_DATA/2017_08_pions/H8_2017_08_27_PL08_TOA/ROOT/Medipix_ 1 285 ../ROOT_FILES/OUTPUT_H8_PIONS_PL08_TOA.root

#./convert_v5 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile1MHz 1 1 ../ROOT_FILES/OutputFile1MHz_HISTO.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile3MHz 1 1 ../ROOT_FILES/OutputFile3MHz_HISTO.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile6MHz 1 1 ../ROOT_FILES/OutputFile6MHz_HISTO.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile12MHz 1 1 ../ROOT_FILES/OutputFile12MHz_HISTO.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile24MHz 1 1 ../ROOT_FILES/OutputFile24MHz_HISTO.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile48MHz 1 1 ../ROOT_FILES/OutputFile48MHz_HISTO.root
#./convert_v5 /home/anatochi/Medipix/H8_DATA/medipix_h8_data/OutputFile96MHz 1 1 ../ROOT_FILES/OutputFile96MHz_HISTO.root
#./convert_v5 /home/anatochi/Medipix/SPS_DATA/MD_2018_06_18/Medipix_ 86 115 ../ROOT_FILES/SPS_MD_2018_06_18.root

#3

#make clean; make analysis_v5;
#./analysis_v5 ../ROOT_FILES/OUTPUT_H8_PIONS_PL08_TOA.root ../ROOT_FILES/OUTPUT_H8_PIONS_PL08_TOA_HISTO.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_1.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_1_HISTO.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_0_$runrunID.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_1_$runrunID.root 1

#make analysis_v5;
#./analysis_v5 ../ROOT_FILES/SPS_MD_2018_06_18.root ../ROOT_FILES/SPS_MD_2018_06_18_HISTO_RP0INT.root 3
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_2.root ../ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_2_HISTO.root 0

#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_1.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_1.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_2.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_2.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_3.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_3.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_4.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_4.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_5.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_5.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_6.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_6.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_7.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_7.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_8.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_8.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_9.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_9.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_10.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_10.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_11.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_11.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_12.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_12.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_13.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_13.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_14.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_14.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_15.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_15.root 0

#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_16.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_16_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_16.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_16_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_17.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_17_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_17.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_17_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_18.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_18_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_18.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_18_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_19.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_19_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_19.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_19_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_20.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_20_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_20.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_20_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_21.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_21_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_21.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_21_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_22.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_22_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_22.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_22_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_23.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_23_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_23.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_23_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_24.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_24_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_24.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_24_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_25.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_25_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_25.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_25_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_26.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_26_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_26.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_26_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_27.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_27_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_27.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_27_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_28.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_28_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_28.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_28_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_29.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_29_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_29.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_29_1.root 1
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_30.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_30_0.root 0
#./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_30.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_15_pions_RUN_HISTO_30_1.root 1


#./analysis_v5 OUTPUT_H8_PIONS_PL08_MPX.root OUTPUT_H8_PIONS_PL08_MPX_HISTO.root 0
#./analysis_v5 OUTPUT_H8_PIONS_QMP54_24_TOA.root OUTPUT_H8_PIONS_QMP54_24_TOA_HISTO.root 0
#./analysis_v5 OUTPUT_H8_PIONS_QMP54_25_TOA.root OUTPUT_H8_PIONS_QMP54_25_TOA_HISTO.root 0
#./analysis_v5 OUTPUT_H8_PIONS_QMP54_26_TOA.root OUTPUT_H8_PIONS_QMP54_26_TOA_HISTO.root 0
#./analysis_v5 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/TMP_CpFM_RUN2.root /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/TMP_RUN2.root 0
#./analysis_v5 ./sps_data/SPS_MD_2017_09_18_LOCAL_PC/SPS_MD_2017_09_18.root ./sps_data/SPS_MD_2017_09_18_LOCAL_PC/SPS_MD_2017_09_18_HISTO.root 1
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_18.root ../ROOT_FILES/SPS_MD_2017_09_18_HISTO_TMPINT.root 0
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_18.root ../ROOT_FILES/SPS_MD_2017_09_18_HISTO_TMPEXT.root 1
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP1INT_dump.root 2
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP3EXT_tmp.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_05_second.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_HISTO_CHIP89_tmp_second.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_05_second.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_HISTO_CHIP87_tmp_second.root 1

#(1) 15:45:50 -- 15:50:33
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP3INT_1.root 0
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP3EXT_1.root 1
#(2) 17:58:47 -- 18:00:45
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP3INT_2.root 0
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP3EXT_2.root 1
#(3) 18:01:51 -- 18:04:29
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP3INT_3.root 0
#./analysis_v5 ../ROOT_FILES/SPS_MD_2017_09_27.root ../ROOT_FILES/SPS_MD_2017_09_27_HISTO_RP3EXT_3.root 1

#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_15Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_15Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_15Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_15Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_25Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_25Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_25Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_25Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_100Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_100Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_100Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_100Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10000Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10000Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_100000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_100000Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_100000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_100000Hz_HISTO_CHIP87.root 1
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP89.root 0
#./analysis_v5 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP87.root 1

#Tracking H8
#make clean; make analysis_tracking;
#./analysis_tracking /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_1.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_1_HISTO.root


#4 data from TMP and CPFM_COBRA
#make clean; make analysis_v7;
#./analysis_v7 ./H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/TMP_RUN2.root ./H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/COBRA_RUN2.root ./H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_2/COBRA_TMP_RUN2.root

#5 For TMPX MOTORS BLM BPM from  SPS MD
#make clean; make analysis_v8;
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_09_18.root ../ROOT_FILES/SPS_MD_2017_09_18_HISTO_TMPEXT_BMP_BLM_MTR_2_tmp.root 1 ../SPS_DATA/MD_2017_09_18/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_09_18/BLM/TIMBER_DATA_BLM5.root ../SPS_DATA/MD_2017_09_18/BPM/TIMBER_DATA_BPM_BV518H.root


#(1) 11:48 -- 12:04
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_1.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508233680000 1508234640000
#(2) 12:35 -- 13:10
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_2.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508236500000 1508238600000
#(3) 13:10 -- 13:30
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_3.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508238600000 1508239800000
#(4) 13:25 -- 13:50
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_4.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508239500000 1508241000000
#(5) 13:50 -- 14:00
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_5.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508241000000 1508241600000
#(6) 14:00 -- 14:15
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_6.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508241600000 1508242500000
#(7) 18:20 -- 18:32
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_7.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508257200000 1508257920000
#(8) 18:34 -- 18:40
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_8.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508258040000 1508258400000
#(9) 18:23:32 -- 18:23:46
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_9.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508257412000 1508257426000
#(10) 18:37:45 -- 18:37:59
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_10.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508258265000 1508258279000
#(11) 18:47 -- 18:56
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_11.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508258820000 1508259360000
#(12) 18:51:00 -- 18:51:40
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_12.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508259060000 1508259100000
#(13) 18:47 -- 18:49
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_13.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508258820000 1508258940000
#(14) 18:51 -- 18:52
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_14.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508259060000 1508259120000
#(15) 18:54 -- 18:56
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_15.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508259240000 1508259360000
#(16) 19:05 -- 19:10
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_16.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508259900000 1508260200000
#(17) 19:15 -- 19:30
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_17.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508260500000 1508261400000
#(18) 19:37 -- 19:47
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_18.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508261820000 1508262420000
#(19) 19:42 -- 19:45
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_19.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508262120000 1508262300000
#(20) 21:40 -- 21:51
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_20.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508269200000 1508269860000
#(21) 21:46 -- 21:50
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_21.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508269560000 1508269800000
#(22) 22:21 -- 22:25
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_22.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508271660000 1508271900000
#(23) 22:42:20 -- 22:45:00
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_23.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508272940000 1508273100000
#(24) 22:16:20 -- 22:19:00
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_24.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508271380000 1508271540000
#(25) 22:45 -- 22:49
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_25.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508273100000 1508273340000
#(26) 23:06 -- 23:09
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_26.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508274360000 1508274540000
#(27) 23:09 -- 23:12
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_27.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508274540000 1508274720000
#(28) 23:20 -- 23:24
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_28.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508275200000 1508275440000
#(29) 23:24 -- 23:28
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_29.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508275440000 1508275680000
#(30) 23:36 -- 23:39
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_30.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508276160000 1508276340000
#(31) 23:39 -- 23:43
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_31.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM1.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508276340000 1508276580000
#(32) 21:49:40 -- 21:50:00
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_32.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508269780000 1508269800000
#(33) 01:14 -- 01:23
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_33.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508282040000 1508282580000
#(34) 18:51 -- 18:53
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_34.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508259060000 1508259180000
#(35) 01:19 -- 01:22
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_35.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508282340000 1508282520000
#(36) 01:21:10 -- 01:21:20
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_36.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508282470000 1508282480000
#(37) 01:24 -- 01:29
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_37.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508282640000 1508282940000
#(38) 02:06 -- 02:10
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_38.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508285160000 1508285400000
#(39) 02:21 -- 02:24
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_39.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508286060000 1508286240000
#(40) 02:34 -- 02:37
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_40.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508286840000 1508287020000
#(41) 02:47 -- 02:49
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_41.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508287620000 1508287740000
#(42) 02:59 -- 03:03
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_42.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508288340000 1508288580000
#(43) 03:13 -- 03:15
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_43.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508289180000 1508289300000
#(44) 03:26 -- 03:31
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_44.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508289960000 1508290260000
#(45) 03:45:10 -- 03:48
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_45.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508291110000 1508291280000
#(46) 04:07 -- 04:09
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_46.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508292420000 1508292540000
#(47) 21:49:40 -- 21:50:00
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_47.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508269780000 1508269800000
#(48) 01:55 -- 01:59
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_48.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM3.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508284500000 1508284740000
#(49) 02:41 -- 02:42:40
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_49.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM6.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508287260000 1508287360000
#(50) 19:43 -- 19:45
#./analysis_v8 ../ROOT_FILES/SPS_MD_2017_10_17.root ../ROOT_FILES/SPS_MD_2017_10_17_HISTO_50.root 2 ../SPS_DATA/MD_2017_10_17/MOTORS/MOTORS.root ../SPS_DATA/MD_2017_10_17/BLM/TIMBER_DATA_BLM4.root ../SPS_DATA/MD_2017_10_17/BPM/TIMBER_DATA_BPM_BV518H.root 1508262180000 1508262300000

#6 data from TMP#89 and TMP#87
#make clean; make analysis_v9;
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_05_HISTO_CHIP89_tmp_second.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_HISTO_CHIP87_tmp_second.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_HISTO_CHIP89_CHIP87_tmp_second.root
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_05_1Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_1Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_1Hz_HISTO_CHIP89_CHIP87_withcut.root
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_05_5Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_5Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_5Hz_HISTO_CHIP89_CHIP87_withcut.root 1
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_05_10Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_10Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_05_10Hz_HISTO_CHIP89_CHIP87_withoutcut.root 0

#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP89_CHIP87_without.root 0
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1Hz_HISTO_CHIP89_CHIP87_with.root 1
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP89_CHIP87_without.root 0
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_5Hz_HISTO_CHIP89_CHIP87_with.root 1
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP89_CHIP87_without.root 0
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_10Hz_HISTO_CHIP89_CHIP87_with.root 1
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP89_CHIP87_without.root 0
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_20Hz_HISTO_CHIP89_CHIP87_with.root 1
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP89_CHIP87_without.root 0
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_50Hz_HISTO_CHIP89_CHIP87_with.root 1
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP89_CHIP87_without.root 0
#./analysis_v9 ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP89.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP87.root ../ROOT_FILES/H8_Test_Beam_2017_12_04_1000000Hz_HISTO_CHIP89_CHIP87_with.root 1

##---------------------------##
## Test-beam H8 OnlyTOA TMP3 ##
##---------------------------##

#1

#make clean; make ascii2root_v6;
#./ascii2root_v6 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_24_TCP67v2_TMP3/Timepix3Data.txt /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_24_TCP67v2_TMP3/Timepix3Data.root
#./ascii2root_v6 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_3/Timepix_CpFM.txt /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_3/Timepix_CpFM.root
#./ascii2root_v6 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_0.txt /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_0.root
#./ascii2root_v6 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_1.txt /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_1.root
#./ascii2root_v6 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_2.txt /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_2.root
#./ascii2root_v6 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_3.txt /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_3.root
#./ascii2root_v6 /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_4.txt /home/anatochi/Medipix/H8_2017_08_pions/H8_2017_08_25_Timepix_CpFM/RUN_4/Timepix_CpFM_4.root

#2
#make clean; make analysis_v6;
#./analysis_v6 ./H8_2017_08_pions/H8_2017_08_24_TCP67v2_TMP3/Timepix3Data.root OUTPUT_H8_PIONS_TCP67v2_TMP3_HISTO.root

#TH1D* h_5_cut = h_5->Clone("h_5_cut");for(int i = 90; i < 150; i++){h_5_cut->SetBinContent(i,0);};
#TF1* f_gaus = new TF1("f_gaus","gaus(0)",0,256);f_gaus->SetParameter(0,4.64602e+02);f_gaus->SetParameter(1,4.53658e+02);f_gaus->SetParameter(2,1.54396e+02);TH1D* h_bg = new TH1D("h_bg","h_bg",256,0,256);for(int i = 0; i < 256; i++){h_bg->SetBinContent(i+1,f_gaus->Eval(i));};TH1D* hh = h_5->Clone("hh");hh->Add(h_bg,-1);for(int i = 0; i < 256; i++){if(hh->GetBinContent(i+1) < 0){hh->SetBinContent(i+1,0);}}
#hh->Integral(219.7-3.0*10.98,219.7+3.0*10.98)


##----------------------##
## Test-beam H8 TRACKER ##
##----------------------##

## 1 ##

#make clean; make ascii2root_trk;
#for runid in $(seq 1 1 127)
#do
#        ./ascii2root_trk /home/anatochi/Medipix/H8_DATA/2018_04_11_pions/RUN_7/Medipix_$runid.dat /home/anatochi/Medipix/H8_DATA/2018_04_11_pions/RUN_7/Medipix_$runid.root
#done

## 2 ##

#make clean; make convert_trk;
#./convert_trk /home/anatochi/Medipix/H8_DATA/2018_04_11_pions/RUN_7/Medipix_ 1 127 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_7.root

## 3 ##

#make clean; make analysis_trk;
#./analysis_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_7.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_7_CLUSTERINFO.root

## 4 ##
#make clean; make trackreco_trk;
#./trackreco_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_7_CLUSTERINFO.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_7_TRACKINFO_RUN_MODE_2.root 2


##-------------------------##
## Test-beam H8 128 Septum ##
##-------------------------##

#for runrunID in $(seq 1 1 66)
#do

    #1

#    cd ../H8_DATA/2018_05_16_pions/RUN_$runrunID;
#    nFiles=$(ls -1q *.dat | wc -l);

#    echo "";
#    echo "Number of files in the directory: "$nFiles;
#    echo "";

#    cd -;

#    make clean; make ascii2root_v5;

#    for runid in $(seq 1 1 $nFiles)
#    do
#            ./ascii2root_v5 ../H8_DATA/2018_05_16_pions/RUN_$runrunID/Medipix_$runid.dat ../H8_DATA/2018_05_16_pions/RUN_$runrunID/Medipix_$runid.root
#    done

    #2

#    make clean; make convert_v5;
#    ./convert_v5 ../H8_DATA/2018_05_16_pions/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_$runrunID.root

    #3

#    make clean; make analysis_v5;
#   ./analysis_v5 /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_05_16_pions_RUN_HISTO_0_$runrunID.root 0

#done

##---------------------------------------------------------##
##
##---------------------------------------------------------##
data_dir=/home/anatochi/Medipix/H8_DATA/2018_04_11_pions;
#data_dir=/home/anatochi/Medipix/SPS_DATA/MD_2018_06_18;

for runrunID in $(seq 3 1 3)
do

    #1

    cd $data_dir/RUN_$runrunID;
    nFiles=$(ls -1q *.dat | wc -l);

    echo "";
    echo "Number of files in the directory: "$nFiles;
    echo "";

    cd -;

    #make clean; make ascii2root_common;
    #for runid in $(seq 1 1 $nFiles)
    #do
    #        ./ascii2root_common $data_dir/RUN_$runrunID/Medipix_$runid.dat $data_dir/RUN_$runrunID/Medipix_$runid.root
    #done

    #2

    #make clean; make convert_common;
    #./convert_common $data_dir/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_$runrunID.root
    #./convert_common $data_dir/RUN_$runrunID/Medipix_ 1 $nFiles /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_RUN_$runrunID.root

    #3

    #make clean; make analysis_common;
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_HISTO_Chip2_RUN_$runrunID.root 2
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_HISTO_RP1I_RUN_$runrunID.root 3
    #./analysis_common /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/MD_2018_06_18_HISTO_RP0I_RUN_$runrunID.root 4

    #4

    #make clean; make analysis_trk;
    #./analysis_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_CLUSTERINFO_RUN_TEMP_$runrunID.root

    #5

    make clean; make trackreco_trk;
    ./trackreco_trk /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_CLUSTERINFO_RUN_$runrunID.root /home/anatochi/Medipix/ROOT_FILES/H8_Test_Beam_2018_04_11_pions_TRACKINFO_RUN_MODE_2_RUN_$runrunID.root 2

done

# From SPS:
#   devRP0E=0 #G02-W0108    FITpix 0384
#   devRP3I=1 #K09-W0255    FITpix 0393
#   devRP3E=2 #C08-W0255    FITpix 0399
#   devRP1I=3 #F04-W0108    FITpix 0409
#   devRP0I=4 #I02-W0108    FITpix 0415


# git commit -am "third commit"
# git push origin develop
