function filteredData=bandpassfilter(MEAdata,cf1,cf2,sampleTime)
    Fs = 1/(sampleTime);
         
    bpFilt = designfilt('bandpassiir', ...       % Response type
       'StopbandFrequency1',0.9*cf1, ...    % Frequency constraints
       'PassbandFrequency1',cf1, ...
       'PassbandFrequency2',0.9*cf2, ...
       'StopbandFrequency2',cf2, ...
       'StopbandAttenuation1',40, ...   % Magnitude constraints
       'PassbandRipple',1,...
       'DesignMethod','cheby2', ...         % Design method
       'MatchExactly','passband', ...
       'SampleRate',Fs);
   
    fvtool(bpFilt)
    for i=1:size(MEAdata,2)
        %zero phase delay filter
        filteredData(:,i)=filter(bpFilt,MEAdata(:,1));
    end

end