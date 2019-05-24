#PART I: DATA ACQUISITION

1. Run executable file: Pixelman.exe
2. For Python plug-in use the following DAQ files:

SPS PC:                         pcen33603.cern.ch
Location:                       C:\Users\Public\Documents\Pixelman_2.2.4_64\scripts\ua9
FileName:                       MD_acquisition_RPall_LastVersion_4 (date modified 27/11/2018)
Configuration files directory:  C:\Users\Public\Documents\Pixelman_2.2.4_64\configs

All necessary parameters are indicated in the scripts for each detector.
3. Data collection for the SPS:
→ On the pcen33603.cern.ch run COPY_DATA_FILES.py (date modified 27/11/2018) script, which is looking for a new file in the output directory and makes a copy to a backup folder and CERN DFS.

#PART II: DATA ANALYSIS

→ Mount a remote CERN DFS directory on the user's PC. On user's PC run autorun_1.bash, it will take a new file from the mounted CERN DFS folder and convert ASCII/CSV raw data file to a ROOT files on the user's PC by means of the ascii2root_common/csv2root_common routine. Or just make a copy of the data acquisition files from the CERN DFS to the user's PC, and run this script.
→ Run autorun_2.bash script, which is doing a merging of all ROOT files into a single ROOT files by means of the convert_common routine. The last step is a running of the analysis_common
 routine for the final plot of all necessary histograms and cluster analysis (CA). To activate CA use a boolean variable _cluster_analysis in the analysis_common routine.
→ Other way to do is to use a run.bash script, which will do the analysis automatically.
→ For particle tracking at the H8: use analysis_trk & trackreco_trk routines.


Medipix Collaboration Picture:

![alt text](https://github.com/nat93/Timepix/blob/master/logo2.png)

The beam profile:

![alt text](https://github.com/nat93/Timepix/blob/master/picture_beam_profile.PNG)

HISTORICAL INFO:

>>>>>>>>>>SPS MD 18-19/06/2018

RUN_0
For the article, good channeling

RUN_1
03:10 - 3:13 19/06/2018
MPX mode, 1 sec, 48 MHz, 100 frames

RUN_2
03:15 - 3:19 19/06/2018
ToA mode, 0.0246 sec, 0.48 MHz, 500 frames

RUN_3
03:20 - 3:24 19/06/2018
ToA mode, 0.00246 sec, 4.8 MHz, 500 frames

RUN_4
03:25 - 3:30 19/06/2018
ToA mode, 0.000246 sec, 48 MHz, 500 frames

RUN_5
03:31 - 3:35 19/06/2018
ToA mode, 0.000123 sec, 96 MHz, 500 frames

RUN_6
03:36 - 3:40 19/06/2018
ToT mode, 0.000246 sec, 48 MHz, 500 frames

RUN_7
20:01 - 6:18 18-19/06/2018
MPX mode, 0.1 sec, 48 MHz ASCII format

RUN_8
08:05 - 15:57 18/06/2018
MPX mode, 0.1/1.0 sec, 48 MHz, new CSV format

>>>>>>>>>> SPS MD 15/08/2018

13:15 Medipix_288.dat, dTHL = 100, Bias (RP1I) = 0 V

13:15 Medipix_289.dat, dTHL = 100, Bias (RP1I) = 10 V

13:20 Medipix_290.dat, dTHL = 100, Bias (RP1I) = 20 V

13:22 Medipix_291.dat, dTHL = 100, Bias (RP1I) = 40 V

13:23 Medipix_292.dat, dTHL = 100, Bias (RP1I) = 60 V

13:25 Medipix_293.dat, dTHL = 100, Bias (RP1I) = 80 V

13:26 Medipix_294.dat, dTHL = 100, Bias (RP1I) = 100 V

13:28 Medipix_295.dat, dTHL = 100, Bias (RP1I) = 40 V, Bias (all the rest) = 0 V

13:30 Medipix_296.dat, dTHL = 100, Bias (RP1I) = 40 V, Bias (all the rest) = 10 V

13:31 Medipix_297.dat, dTHL = 100, Bias (RP1I) = 40 V, Bias (all the rest) = 20 V

13:33 Medipix_298.dat, dTHL = 100, Bias (RP1I) = 40 V, Bias (all the rest) = 40 V

13:35 Medipix_299.dat, dTHL = 100, Bias (RP1I) = 40 V, Bias (all the rest) = 60 V

13:36 Medipix_300.dat, dTHL = 100, Bias (RP1I) = 40 V, Bias (all the rest) = 80 V

13:38 Medipix_301.dat, dTHL = 100, Bias (RP1I) = 40 V, Bias (all the rest) = 100 V

13:41 Medipix_302.dat, dTHL (RP1I) = 100, dTHL (all the rest) = 40, Bias = 40 V

13:46 Medipix_303.dat, dTHL (RP1I) = 40, dTHL (all the rest) = 100, Bias = 40 V


