%% Struct
%File Details
details.FileNameInput = FileName;
details.inputdir = inputdir;
details.FileNameOutput = excelFileName;
details.frequency = frequency;
%User Input and Hard-coded values
details.spikeThreshold = (userInput(1));
details.distanceSpike = distanceSpike;
details.artifactThreshold = userInput(2);
details.distanceArtifact = distanceArtifact;
details.minSLEduration = minSLEduration;
%Detect events (epileptiform/artifacts) | Absolute Values
details.minPeakHeightAbs = minPeakHeight; 
details.minPeakDistanceAbs = minPeakDistance;
details.minArtifactHeightAbs = minArtifactHeight;
details.minArtifactDistanceAbs = minArtifactDistance;
%SLECrawler 
details.durationOnsetBaseline = durationOnsetBaseline;
details.durationOffsetBaseline = durationOffsetBaseline;
details.calculateMeanOffsetBaseline = calculateMeanOffsetBaseline;

%detection results
details.IISsDetected = numel(IIS(:,1));
details.eventsDetected = numel(events(:,1));
details.SLEsDetected = numel (SLE_final(:,1);

%LED details
if userInput(4)>0
    details.stimulusChannel = userInput(4);
    details.onsetDelay = onsetDelay;
    details.offsetDelay = offsetDelay;
end


