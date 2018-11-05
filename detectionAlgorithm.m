%Program: Epileptiform Activity Detector
%Author: Michael Chang (michael.chang@live.ca)
%Copyright (c) 2018, Valiante Lab
%Version 8.1 (In Vivo)

%Description: Standard stage 1 to detect epileptiform events. THe script
%can process a single file that you select or a collection of files from a
%directory that you have provided. Additionally, this script will save the
%work space as a .mat file which contains the time points of the seizure
%onset and offset and classificaitons, as well as details regarding how the
%files were detected. These variables will be saved as a struct (.mat
%file).

%% Stage 1: Detect Epileptiform Events
%clear all (reset)
close all
clear all
clc

%Manually set File Directory
inputdir = 'F:\ictal segments\young #4';

%GUI to set thresholds
%Settings, request for user input on threshold
titleInput = 'Specify Detection Thresholds';
prompt1 = 'Epileptiform Spike Threshold: average + (3.9 x Sigma)';
prompt2 = 'Artifact Threshold: average + (70 x Sigma)';
prompt3 = 'Figure: Yes (1) or No (0)';
prompt4 = 'Stimulus channel (enter 0 if none):';
prompt5 = 'Troubleshooting: plot SLEs(1), IIEs(2), IISs(3), Artifacts (4), Review(5), all(6), None(0):';
prompt6 = 'To analyze multiple files in folder, provide File Directory:';
prompt = {prompt1, prompt2, prompt3, prompt4, prompt5, prompt6};
dims = [1 70];
definput = {'10', '70', '0', '2', '0', ''};

opts = 'on';    %allow end user to resize the GUI window
InputGUI = (inputdlg(prompt,titleInput,dims,definput, opts));  %GUI to collect End User Inputs
userInput = str2double(InputGUI(1:5)); %convert inputs into numbers

if (InputGUI(6)=="")
    %Load .abf file (raw data), analyze single file
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', inputdir);%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
    [spikes, events, SLE, details, artifactSpikes] = detectionInVivo(FileName, userInput, x, samplingInterval, metadata);
    save(sprintf('%s.mat', FileName(1:end-4)))  %Save Workspace  
else
    % Analyze all files in folder, multiple files
    PathName = char(InputGUI(6));
    S = dir(fullfile(PathName,'*.abf')); 

    for k = 1:numel(S)
        clear IIS SLE_final events fnm FileName x samplingInterval metadata %clear all the previous data analyzed
        fnm = fullfile(PathName,S(k).name);
        FileName = S(k).name;
        [x,samplingInterval,metadata]=abfload(fnm);
        [spikes, events, SLE, details, artifactSpikes] = detectionInVivo(FileName, userInput, x, samplingInterval, metadata);               
        save(sprintf('%s.mat', FileName(1:end-4)))  %Save Workspace    
    end
end

fprintf(1,'\nThank you for choosing to use Michaels Epileptiform Activity Detector.\n')

   


