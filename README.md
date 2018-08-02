
1. To convert each ASCII file to a sinle ROOT file (For H8 and SPS):          ascii2root_common

2. To merge all single ROOT files to a single ROOT file (For H8 and SPS):     convert_common

3. Frame analysis for a single chip (For H8 and SPS):                         analysis_common

4. Tracking with more than one plane (For H8 and only for 4 planes)

    4.1 To extract the info about each cluster:                             analysis_trk

    4.2 To performe a track reconstruction:                                 trackreco_trk


To run the acquisition and data analysis scripts during SPS MD:

--> On the PCEN33603 PC:
Switch on the Pixelman and run the Python scripting utilit with DAQ script MD_acquisition_RPall_LastVersion.
Via WINDOWS CMD: ---$ python COPY_DATA_FILES.py --> it will check the new file and copy it to the DFS and to the backup folder.

--> On the user PC:
Mount the DFS directory with output files.
Launch data analysis scripts: (1) autorun_1.bash (continuously) (2) autorun_2.bash (for defined files)
