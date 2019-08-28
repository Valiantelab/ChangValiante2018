%Program: Epileptiform Event Detector 
%Corresponding Author: Michael Chang (michael.chang@uhnresearch.ca) 
%Copyright (c) 2018, Valiante Lab
%Version 8.6 

% For quick start: i) Run the script, ii) click OK, iii) select the 
% .abf file to analyze; for demo select "13226009(exampleFile).abf" 

%% Stage 1: Detect Epileptiform Events
%clear all (reset)
close all
clear all
clc

%Add all subfolders in working directory to the path.
addpath(genpath(pwd));  

%Manually set File Directory
inputdir = 'C:\Users\Michael\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\Control Data\1) Control (VGAT-ChR2, light-triggered)\1) abf files';

%GUI to set thresholds
%Settings, request for user input on threshold
titleInput = 'Specify Detection Parameters';
prompt1 = 'Epileptiform Spike Threshold: average + (3.9 x Sigma)';
prompt2 = 'Artifact Threshold: average + (70 x Sigma)';
prompt3 = 'Figure: Yes (1) or No (0)';
prompt4 = 'Stimulus channel (enter 0 if none):';
prompt5 = 'Troubleshooting: plot SLEs(1), IIEs(2), IISs(3), Artifacts (4), Review(5), all(6), None(0):';
prompt6 = 'To analyze multiple files, provide the folder directory (leave blank to select individual files):';
prompt = {prompt1, prompt2, prompt3, prompt4, prompt5, prompt6};
dims = [1 70];
definput = {'4', '100', '1', '2', '0', ''};

opts = 'on';    %allow end user to resize the GUI window
InputGUI = (inputdlg(prompt,titleInput,dims,definput, opts));  %GUI to collect End User Inputs
userInput = str2double(InputGUI(1:5)); %convert inputs into numbers

if (InputGUI(6)=="")
    %Load .abf file (raw data), analyze single file
    [FileName,PathName] = uigetfile ('*.abf','pick .abf file', inputdir);%Choose abf file
    [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
    [spikes, events, SLE, details] = detectionInVitro4AP(FileName, userInput, x, samplingInterval, metadata);
else
    % Analyze all files in folder, multiple files
    PathName = char(InputGUI(6));
    S = dir(fullfile(PathName,'*.abf'));

    for k = 1:numel(S)
        clear IIS SLE_final events fnm FileName x samplingInterval metadata %clear all the previous data analyzed
        fnm = fullfile(PathName,S(k).name);
        FileName = S(k).name;
        [x,samplingInterval,metadata]=abfload(fnm);
        [spikes, events, SLE, details] = detectionInVitro4AP(FileName, userInput, x, samplingInterval, metadata);
        %Collect the average intensity ratio for SLEs
        %indexSLE = events(:,7) == 1;
        %intensity{k} = events(indexSLE,18);                   
    end
end

fprintf(1,'\nA summary of the detection results can be found in the current working folder: %s\n', pwd)
fprintf(1,'\nThank you for choosing to use Chang & Valiante (2018) Epileptiform Event Detector.\n')

   


