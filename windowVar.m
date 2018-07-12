function varTime=windowVar(data,Fs,winSize,stepSize)
    dataLen=size(data,1);
    duration=dataLen;
    winSize=winSize*Fs;
    stepSize=stepSize*Fs;
    numStep=round((duration-winSize)/stepSize);
    
    varTime=zeros(2,numStep);
    for i=1:numStep
        startIndex=round(stepSize*(i-1)+1);
        endIndex=round(stepSize*(i-1)+winSize);
        %return time and value
        varTime(1,i)=(endIndex+startIndex)/(2*Fs);
        varTime(2,i)=var(data(startIndex:endIndex));
    end
    
end