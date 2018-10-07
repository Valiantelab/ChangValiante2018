function [label,classification] = decipher (events, i)
%decipher is a function to convert my number classification system into
%normal language
%   This function turns the numbers into words for the purposes of
%   labelling the epileptiform events that are plotted for troubleshooting
%   and demonstration purposes. The numbers are classification codes for
%   the type of epileptiform event and if tonic-phase is present.

%% Label the epiletiform event (for plotting); convert numbers to words
if events(i,7) == 1 
    label = 'SLE';
else if events(i,7) == 2
        label = 'IIE';
    else if events(i,7) == 3
            label = 'IIS';
        else if events(i,7) == 4
                label = 'artifact';
            else if events (i,7) == 1.5
                label = 'Questionable SLE';
                else if events (i,7) == 2.5
                        label = 'Your hypothesis failed: IIS w/ tonic';
                    else                                    
                    label = 'error; review the function decipher.m to see if you missed a class of epileptiform events';
                    end                                
                end                                                    
            end
        end
    end
end  

%% classification of the epileptiform event based on its tonic-phase         
if events(i,15) == 0
classification = 'No tonic phase';
else if events(i,15) == 1
        classification = 'Tonic-Clonic SLE';
    else if events(i,15) == 2
            classification = 'Tonic-only SLE';
        else
            classification = 'IIS (duration <3 s)';
        end
    end
end   


