%Program: Epileptiform Event Analyzer (aka Stage 2)
%Corresponding Author: Michael Chang (michael.chang@uhnresearch.ca) 
%Copyright (c) 2018, Valiante Lab
%Version FrequencyAnalysis V1.0 

%Description: Standard stage 2 to analyze your epileptiform events. The script
%can process a single .mat file that you select or a collection of .mat files from a
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

%Add all subfolders in working directory to the path.
addpath(genpath(pwd));  

%Manually set File Directory
inputdir = 'C:\Users\Michael\OneDrive - University of Toronto\3) Manuscript III (Nature)\Section 2\B. HEPES (VGAT-ChR2)';

%GUI to set thresholds
%Settings, request for user input on threshold
titleInput = 'Specify Spectrogram Parameters';
prompt1 = 'Window Size (sec)';
prompt2 = 'Overlap between windows (sec)';
prompt3 = 'Figure: Yes (1) or No (0)';
prompt4 = 'Test Group:';
prompt5 = 'Unique Title';
prompt6 = 'To analyze multiple files in folder, provide File Directory:';
prompt = {prompt1, prompt2, prompt3, prompt4, prompt5, prompt6};
dims = [1 70];
definput = {'1.1', '.1', '1', 'HEPES-buffered ACSF (VGAT-ChR2)', 'result_manuscriptIII.xlsx', ''};


opts = 'on';    %allow end user to resize the GUI window
InputGUI = (inputdlg(prompt,titleInput,dims,definput, opts));  %GUI to collect End User Inputs
% userInput = str2double(InputGUI(1:3)); %convert inputs into numbers

if (InputGUI(6)=="")
    %Load .mat file (raw data), analyze single file    
    [FileName,PathName] = uigetfile ('*.mat','pick .mat file to load Workspace', inputdir);%Choose file    
    [epileptiformEvent, frequencyContentAnalysis] = frequencyAnalysis(InputGUI, PathName, FileName); 
        
%     [x,samplingInterval,metadata]=abfload([PathName FileName]); %Load the file name with x holding the channel data(10,000 sampling frequency) -> Convert index to time value by dividing 10k
%     [spikes, events, SLE, details, artifactSpikes] = detectionInVitro4AP(FileName, userInput, x, samplingInterval, metadata);
    save(sprintf('%s(stage2).mat', FileName(1:end-4)))  %Save Workspace  
else
    % Analyze all .mat files in folder, multiple files
    PathName = char(InputGUI(6));
    S = dir(fullfile(PathName,'*.mat'));

    for k = 1:numel(S)
        clear FileName epileptiformEvent frequencyContentAnalysis %clear all the previous data analyzed
%         fnm = fullfile(PathName,S(k).name);
        FileName = S(k).name;
        [epileptiformEvent, frequencyContentAnalysis] = frequencyAnalysis(InputGUI, PathName, FileName); 
%         [x,samplingInterval,metadata]=abfload(fnm);
%         [spikes, events, SLE, details, artifactSpikes] = detectionInVitro4AP(FileName, userInput, x, samplingInterval, metadata);
        save(sprintf('%s(stage2).mat', FileName(1:end-4)))  %Save Workspace 
        %Collect the average intensity ratio for SLEs
        %indexSLE = events(:,7) == 1;
        %intensity{k} = events(indexSLE,18);                   
    end
end

fprintf(1,'\nA summary of the detection results can be found in the current working folder: %s\n', pwd)
fprintf(1,'\nThank you for choosing to use Chang & Valiante (2018) Epileptiform Event Analyzer.\n')
