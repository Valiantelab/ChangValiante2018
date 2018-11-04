%% Stage 1: Import .Mat Files (workspace)
%clear all (reset)
close all
clear all
clc

%Manually set File Directory
inputdir = 'C:\Users\Michael\OneDrive - University of Toronto\8) Seizure Detection Program\Workspace\Nov2_2018\light-triggered';

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
definput = {'3.9', '70', '0', '2', '0', ''};

opts = 'on';    %allow end user to resize the GUI window
InputGUI = (inputdlg(prompt,titleInput,dims,definput, opts));  %GUI to collect End User Inputs
userInput = str2double(InputGUI(1:5)); %convert inputs into numbers

if (InputGUI(6)=="")    
    %Load .abf file (raw data), analyze single file
    [FileName,PathName] = uigetfile ('*.mat','pick .mat file to load Workspace', inputdir);%Choose file    
    fnm = fullfile(PathName,FileName);
    myVars = {'spikes', 'events', 'SLE', 'artifactSpikes', 'details', 'samplingInterval', 'x', 'metadata'};
    load(sprintf('%s', fnm), myVars{:})  
    
else
    % Analyze all files in folder, multiple files
    PathName = char(InputGUI(6));
    S = dir(fullfile(PathName,'*.mat'));

    for k = 1:numel(S)
        fnm = fullfile(PathName,S(k).name);
        FileName = S(k).name;
        myVars = {'spikes', 'events', 'SLE', 'artifactSpikes', 'details', 'samplingInterval', 'x', 'metadata'};
        load(sprintf('%s', fnm), myVars{:})
        
        detected.spikes = spikes;
        detected.events = events;
        detected.SLE = SLE;
        detected.artifactSpikes = artifactSpikes;
        detected.details = details;
        detected.samplingInterval = samplingInterval;
        detected.x = x;
        detected.metadata = metadata;
        
        save(sprintf('%s.mat', FileName(1:8)), 'detected')  %Save Workspace    
        
    end
end