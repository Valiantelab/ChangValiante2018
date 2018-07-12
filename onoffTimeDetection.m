function [onlocs,offlocs,seizureSeg]=onoffTimeDetection(data,sampleTime)
    Fs=1/sampleTime;
    
    %sliding window and sliding steps of the variance calculation
    slidingsteps=0.5;
    slidingwindow=3;
    %find the min within 60s
    baselinecal=60;
    
    timeextension=round(baselinecal./slidingsteps);
    
    %calculate the vairance along the moving window
    varTime=windowVar(data,Fs,slidingwindow,slidingsteps);
    sizeVarTime=size(varTime,2);
    
    %detect peaks above certain threshold
    [~,peaklocs] =findpeaks(varTime(2,:),'MinPeakHeight',0.1e-3);

    index=1;
    seizure=0;
    onlocs=[];
    offlocs=[];
    seizureSeg=[];
    minthresholdcoeff=0;
    while index<size(peaklocs,2)
        %location of the on time
        templocs= peaklocs(index);
   
        
        sortEndIndex=templocs+timeextension;   
        minthresholdcoeff=0.25;
        
        %cropped the variance after the on time for latter calculation and
        %sort the variance vlaue from low to high
        if (sortEndIndex)<sizeVarTime
            sortedData=sort(varTime(2,templocs:sortEndIndex));
        else
            sortedData=sort(varTime(2,templocs:end));
        
        end
        
        % calculate the median of the last 25 percent of the lowest data.
        % Use this vlaue as the threshold for detecting off time
        minthreshold=median(sortedData(1,1:round(size(sortedData,2).*minthresholdcoeff)));
        offindex=find(varTime(2,templocs:end)<minthreshold,1,'first')+templocs;
        
        
        
        peakvalue=varTime(2,templocs);
        %the condition for detecting off time for seizure event is
        %different due to the large amplitude of the signal
        if peakvalue>0.002 && ((peakvalue/minthreshold)<30)
            seizure=1;
            sortEndIndex=templocs+timeextension*3;
            minthresholdcoeff=0.25;
            
                    
            if (sortEndIndex)<sizeVarTime
                sortedData=sort(varTime(2,templocs:sortEndIndex));
            else
                sortedData=sort(varTime(2,templocs:end));

            end
            
            minthreshold=median(sortedData(1,1:round(size(sortedData,2).*minthresholdcoeff)));
            
            while (offindex-templocs)<timeextension
                offindex=find(varTime(2,offindex+1:end)<minthreshold,1,'first')+offindex;
            end
        end
        
        
        
        %peak value of the on time should be three time larger than off
        %time or else its not a burst or a seizure(continue looking for off
        %time)
        if peakvalue>(minthreshold*3)
            if(size(offindex,2)==0)
                break
            end
            %append on and off time in th array
            onlocs=[onlocs templocs];
            offlocs=[offlocs offindex];
            % set the index of the next peak
            index=find(peaklocs>offindex,1,'first');
            
            if seizure
                seizureSeg=[seizureSeg [templocs;offindex]];
      
            end
            
        else
            index=index+1;
        end
        seizure=0;
        
    end
    onlocs=onlocs.*slidingsteps;
    offlocs=offlocs.*slidingsteps;
    seizureSeg=seizureSeg.*slidingsteps;
    figure;
    plot(varTime(1,:),varTime(2,:))
    hold on 
    plot(onlocs,zeros(1,size(onlocs,2)),'o')
    plot(offlocs,zeros(1,size(offlocs,2)),'o')
    
    
end