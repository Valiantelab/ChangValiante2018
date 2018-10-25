function [label,classification] = decipher (events, i)
%decipher is a function to convert my number classification system into
%normal language. Author: Michael Chang, Liam Long, and Thomas Lordello
%   This function turns the numbers into words for the purposes of
%   labelling the epileptiform events that are plotted for troubleshooting
%   and demonstration purposes. The numbers are classification codes for
%   the type of epileptiform event and if tonic-phase is present.

%If default values are not provided
if nargin <2
    i = 1;
end

%% Label the epiletiform event (for plotting); convert numbers to words
switch events(i,7) 
    case 0
        label = 'Review Event';
    case 1 
        label = 'SLE';
    case 2
        label = 'IIE';
    case 3
        label = 'IIS';
    case 4
        label = 'artifact';
    case 1.5
        label = 'Questionable SLE';
    otherwise                                    
        label = 'Error: review the function decipher.m to see if you missed a class of epileptiform events';
end



%% classification of the epileptiform event based on its tonic-phase         
switch events(i,13) 
    case -1
        classification = 'IIS (duration <3 s)';
    case 0
        classification = 'No tonic phase';
    case 1
        classification = 'Tonic-Clonic SLE';
    case 2
        classification = 'Tonic-only SLE';
    otherwise
        classification = 'Error: review the function decipher.m; perhaps you discovered a new ictal phase ';
end   


