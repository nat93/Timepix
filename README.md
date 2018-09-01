
1. To convert each ASCII file to a sinle ROOT file (For H8 and SPS):          ascii2root_common

2. To merge all single ROOT files to a single ROOT file (For H8 and SPS):     convert_common

3. Frame analysis for a single chip (For H8 and SPS):                         analysis_common

4. Tracking with more than one plane (For H8 and only for 4 planes)

    4.1 To extract the info about each cluster:                             analysis_trk

    4.2 To performe a track reconstruction:                                 trackreco_trk

Medipix Collaboration Picture:

![alt text](https://github.com/nat93/Timepix/blob/master/logo2.png)

To run the acquisition and data analysis scripts during SPS MD:

--> On the PCEN33603 PC:
Switch on the Pixelman and run the Python scripting utilit with DAQ script MD_acquisition_RPall_LastVersion.
Via WINDOWS CMD: ---$ python COPY_DATA_FILES.py --> it will check the new file and copy it to the DFS and to the backup folder.

--> On the user PC:
Mount the DFS directory with output files.
Launch data analysis scripts: (1) autorun_1.bash (continuously) (2) autorun_2.bash (for defined files)

For SPS MD data analysis use sps_md_analysis.C to calculate the ration between CH and DCH particles.

The neutron beam profile from nTOF experiment:

![alt text](https://github.com/nat93/Timepix/blob/master/nTOF_UPSTREAM_QUADPIX_2018_08_30.PNG)

SPS MD 18-19/06/2018

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

SPS MD 15/08/2018

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


