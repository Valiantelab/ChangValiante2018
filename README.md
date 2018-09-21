# Epileptiform-Event-Detection-Algorithm (detection.m)

Description: 
High speed detection of ictal events (SLEs), as Chang et al., 2018 would mark them. This simple script detects all the spikes in the time series provided (.abf file) and groups spikes that are within 10 sec of each other as one event. These events are then classified based on their spiking characfteristics (duration, rate, intensity, and amplitude). Events can be classified as ictal event (SLE), interictal event (IIE), interictal spike (IIS), or an artifact.

For demonstration: 
open exampleScript.m (which is exactly the same script as detection.m), but with on screen instructions provided. Provide the .abf file, "exampleFile.abf".

Input options:
This script gives end user the option to set the threshold for epileptiform spike detection as some multiple of the time series's sigma (Default is 3.9x, but use up to 10x in noisier data sets). End user can also set the threshold for artifacts in a similar manner (default is 70x; future versions will detect artifacts by their width). There is also an option to request figures of detected SLEs and all events (for troubleshooting purposes); select 1 as the option (see figure below). Otherwise, click "OK" and use default options to proceed. The script will then request for a time series (LFP recording) to analyze. Select "exampleFile" as indicated by on-screen instuctions.

Output Files:
The script outputs an excel sheet and powerpoint of detect SLEs and figures demonstrating how the events were segmented based on their features (amplitude, frequency, intensity, and duration). These output files will be placed into the MatLab working folder. The excel sheet reports all the detected events and their features, organized by category of event in different tabs. 



Thank you for choosing to use the Epileptiform Event Detector for your research needs. 

Refer to Chang et al., 2018. Neurobiology of Disease.

Authors: Michael Chang (michael.chang@live.ca), Christopher Lucasius, Fu-der (Fred) Chen, Liam Long, Thomas Lordello, Vitaly Topekha, Taufik A. Valiante.
 

Copyright 2018, Valiante Lab 
