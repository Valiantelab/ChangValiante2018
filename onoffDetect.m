function [onloc,offloc]=onoffDetect(rawData,startTime,endTime,sampleTime)


    startIndex=ceil(startTime/1e-4);
    endIndex=floor(endTime/1e-4);

%     [rawData,sampleTime,h]=abfload(filename);
    time=1:size(rawData,1);
    time=time'.*sampleTime;

    %bandpass filter CF1: low cutoff, CF2:high cutoff (in Hz)
    cf1=1;
    cf2=300;

    filteredData=bandpassfilter(rawData,cf1,cf2,sampleTime);
    croppedData=filteredData(startIndex:endIndex,1);
    croppedTime=time(startIndex:endIndex);



    [onlocs,offlocs,seizureSeg]=onoffTimeDetection(croppedData,sampleTime);
    seizureSeg=seizureSeg'+startTime-1;
    figure
    plot(croppedTime,croppedData)
    hold on 
    plot(onlocs+startTime-1,zeros(1,size(onlocs,2)),'o')
    plot(offlocs+startTime-1,zeros(1,size(offlocs,2)),'o')
    onlocs=(onlocs+startTime-1)';
    offlocs=(offlocs+startTime-1)';
    
    
end