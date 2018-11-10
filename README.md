# Epileptiform-Event-Detector(detectionAlgorithm.m)

Requirements:
- MatLab R2015a (or later)
- Microsoft Office

Quick start (For demonstration):
1) Check MatLab's working directory is set to ChangValiante2018.
2) Open detectionAlgorithm.m
3) Run the script
4) Specify Detection Parameters, or click "OK" and use default settings 
5) Select .abf you want to analyze, or select 13226009(exampleFile).abf provided

Description:
High speed detection of ictal events (SLEs), as Chang et al., 2018 would mark them. This simple script detects all the spikes in the time series provided (.abf file, channel 1) and groups spikes that are within 10 sec of each other as one event. These events are then classified based on their spiking characteristics (duration, rate, intensity, and amplitude). Events can be classified as ictal event (SLE), interictal event (IIE), interictal spike (IIS), or an artifact.

Input options (detailed):
This script gives end user the option to set the threshold for epileptiform spike detection as some multiple of the time series's sigma (Default is 3.9x, but use up to 10x in noisier data sets or in vivo recordings). End user can also set the threshold for artifacts in a similar manner (default is 70x; future versions will detect artifacts by their width). There is also an option to request figures of detected SLEs and other events (for troubleshooting purposes). There is also a option. The script will then request for a time series (LFP recording) to analyze. Select the .abf file you want to analyze as indicated by on-screen instuctions.

Recommendations for Threshold setting:
Use the recommended epileptiform spike threshold: 3.9 x sigma, for in vitro recordings. However, use a higher threshold: 10 x sigma, for noisy data, results where multiple seizure-like events are grouped together, or in vivo recordings.  

Output Files:
The script outputs an excel sheet and powerpoint of detect SLEs and figures demonstrating how the events were segmented based on their features (amplitude, frequency, intensity, and duration). These output files will be placed into the MatLab working folder. The excel sheet reports all the detected events and their features, organized by category of event in different tabs.

Thank you for choosing to use the Valiante Lab's Epileptiform Event Detector for your research needs.

Refer to Chang et al., 2018. Neurobiology of Disease.
Link: https://www.sciencedirect.com/science/article/pii/S0969996117302255

Authors: Michael Chang (michael.chang@live.ca), Christopher Lucasius, Liam Long and Taufik A. Valiante.

Acknowledgements: Fu-der (Fred) Chen, Thomas Lordello, Vitaly Topekha, Kramay Patel, Adam Gierlach, and Gerard O'Leary.  

Copyright 2018, Valiante Lab 
