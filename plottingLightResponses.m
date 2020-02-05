
%organize groups index
indexControl = events(:,4)==1;
indexTest = events(:,4)==2;
indexPosttest = events(:,4)==3;
indexWashout = events(:,4) == 4;    %optional, may not be present

%% Analyze Light-triggered response in tissue
    LightResponse_window = 1;    %seconds
    
        
    %Pre allocate matrix to store values
    lightResponse = zeros (numel(P.range(:,1)),6);
        
    for i = 1:numel(P.range(:,1))        
        %% Find indices (for light-triggered responses in the tissue)
        if P.range(i,1)+(LightResponse_window*frequency) < numel(LFP)        
            lightResponse_index = (P.range(i,1):P.range(i,1)+(LightResponse_window*frequency));
        else
            continue
        end
        
        %Store indices for later analysis
        lightResponse_indices{i} = lightResponse_index;
        
        
        %% Analysis
        %Calculate max deflection from starting point when light is delievered
        lightResponseCentered = LFP(lightResponse_indices{i}) - LFP(lightResponse_indices{i}(1));        %Centered the light response so starting point is 0
        [maxDeflectionValue, maxDeflectionLocation] = max(abs(lightResponseCentered)); %Find the maximum distance from the starting point
        maxDeflection = lightResponseCentered(maxDeflectionLocation); %will account for the deflection is positive or negative 
        
        %Calculate power of Light-triggered Response
        totalEnergy = sum(powerFeature(lightResponse_indices{i})); %This is the energy of the event; totalPower is a misnomer
        power = totalEnergy /LightResponse_window;  %This is the power of the event; you previously considered it as the average power/minute

        %Calculate peak-to-peak amplitude of Light-triggered Response
        eventVectorLFP = LFP_centered(lightResponse_indices{i});
        p2pAmplitude = max(eventVectorLFP) - min (eventVectorLFP);

        %% Store Output
        lightResponse(i,1) = P.range(i,1)/frequency;    %store the onset time 
        lightResponse(i,2) = lightResponse(i,1)+ LightResponse_window;   %offset of spike window
        lightResponse(i,3) = maxDeflection;  %Max deflection (mV)               
        %lightResponse(i,4) = 
        lightResponse(i,5) = power; %Power (mV^2/s)
        lightResponse(i,6) = p2pAmplitude; %Amplitude (mV)    
    end    

%     figure
%     plot(t(lightResponse_indices{i}), LFP(lightResponse_indices{i}))

    %% finding the light pulse of interest | Automated or User Specified        
    %control period - start
    if userInput_time_spec(2) > -0.1        
        startControl = userInput_time_spec(2);
    else
        startControl= events(indexControl,1)-0.15;  %move the start time up earlier so you definitely include the light pulse
    end
    
    %control period - end
    if userInput_time_spec(3) > -0.1        
        endControl = userInput_time_spec(3);
    else
        endControl = startControl(end);
    end
    
    %Locate light pulses during control period specified
    controlCondition = lightResponse(:,1) > startControl(1) & lightResponse(:,1) < endControl;   %logical index of all the light responses during control condition
    lightResponse(controlCondition, 4) = 1;
        
    %Test period - start   
    if userInput_time_spec(5) > -0.1
        startTest = userInput_time_spec(5);
    else
        startTest= events(indexTest,1)-0.15;
    end
    
    %Test period - end
    if userInput_time_spec(6) > -0.1
        endTest = userInput_time_spec(6);
    else
        endTest = startTest(end);
    end

    %Locate light pulses during test period specified
    testCondition = lightResponse(:,1) > startTest(1) & lightResponse(:,1) < endTest;
    lightResponse(testCondition, 4) = 2;
    
    %Posttest period - start
    if userInput_time_spec(8) > -0.1
        startPosttest = userInput_time_spec(8);
    else
        startPosttest = events(indexPosttest,1) - 0.15; 
    end
    
    %Posttest period - end
    if userInput_time_spec(9) > -0.1
        endPosttest = userInput_time_spec(9);
    else
        endPosttest = startPosttest(end);
    end
    
    %Locate light pulses during posttest period specified 
    posttestCondition = lightResponse(:,1) > startPosttest(1) & lightResponse(:,1) < endPosttest;
    lightResponse(posttestCondition, 4) = 3;

    %Optional: Washout Period 
    if ~isempty(InputGUI_time_spec{10})
        % Washout - Start
        if userInput_time_spec(11) > -0.1
            startWashout = userInput_time_spec(11);
        else
            startWashout = events(indexWashout,1) - 0.15;
        end
        
        %washout - End
        if userInput_time_spec(12) > -0.1
            endWashout = userInput_time_spec(12);
        else
            endWashout = startWashout(end);
        end
        
        %Locate light pulses during washout period specified
        washoutCondition = lightResponse(:,1) > startWashout(1) & lightResponse(:,1) < endWashout;
        lightResponse(washoutCondition, 4) = 4;
    end           
    
    %% Remove the light pulses during ictal events    
    %Find index of all ictal events of interest (that are labeled within a group, excel sheet)
    indexIctalEvents = events(:,4)>0;
    
    startIctalEvent= events(indexIctalEvents,1)-0.15;  %move the start time up earlier so you definitely capture some time before the light pulse
    endIctalEvent = events(indexIctalEvents,2) + 10;   %move the end time back, so you don't capture post-ictal bursts
    
    %remove spikes that trigger an ictal eventt or occur during the ictal events of interest
    for i=1:sum(indexIctalEvents)
        for j = 1:numel(lightResponse(:,1))
            if lightResponse(j,1) > startIctalEvent(i) & lightResponse(j,1) < endIctalEvent(i)
                lightResponse(j,4) = 0.1;
            end
        end
    end         

    
    %% Sort light-response into groups 
    index_LR_Control = lightResponse(:,4)==1;
    index_LR_Test = lightResponse(:,4)==2;
    index_LR_Posttest = lightResponse(:,4)==3;
    index_LR_Washout = lightResponse(:,4)==4;
 
%     %Amplitude (mV)
%     feature = 3;
% 
%     %Form Groups
%     amplitudeControl = lightResponse(lightResponse(:,4)==1,feature);
%     amplitudeTest = lightResponse(lightResponse(:,4)==2,feature);
%     amplitudePosttest = lightResponse(lightResponse(:,4)==3,feature);
%     amplitudeWashout = lightResponse(lightResponse(:,4)==4,feature);

%     %Test for Normality
%     [resultsAmplitude(1,:)] = stage3Analysis (amplitudeControl, 'parametric', 'figure');
%     [resultsAmplitude(2,:)] = stage3Analysis (amplitudeTest, 'parametric', 'figure');
%     [resultsAmplitude(3,:)] = stage3Analysis (amplitudePosttest, 'parametric', 'figure');
%     [resultsAmplitude(4,:)] = stage3Analysis (amplitudeWashout, 'parametric', 'figure');

%     %durationMatrix
%     LR_AmplitudeMatrix(1:numel(amplitudeControl),1) = amplitudeControl;
%     LR_AmplitudeMatrix(1:numel(amplitudeTest),2) = amplitudeTest;
%     LR_AmplitudeMatrix(1:numel(amplitudePosttest),3) = amplitudePosttest;
%     LR_AmplitudeMatrix(1:numel(amplitudeWashout),4) = amplitudeWashout;
%     LR_AmplitudeMatrix(LR_AmplitudeMatrix==0) = NaN;

%     %Comparisons - Kruskal Wallis
%     p = kruskalwallis(LR_AmplitudeMatrix);
%     title('Boxplot: duration of ictal events from different treatment groups')
%     xlabel ('Treatment Group')
%     ylabel ('Duration (s)')

    
    %Prep data for figures
    LightResponse_window_indices = numel(lightResponse_index);   %Number of indices in window

    %Preallocate matrix to store LFP, snipits of light response
    LFP_LightResponse = zeros(numel(lightResponse_indices),LightResponse_window_indices);
    
    %Store all the snippits of the LFP light response
    for i = 1:numel(lightResponse_indices)
        LFP_LightResponse(i,:) = LFP(lightResponse_indices{i})-LFP(lightResponse_indices{i}(1));  %Light Response centered        
    end
    
    %Form average waveforms of light response
    avg_control_LightResponse = mean(LFP_LightResponse(index_LR_Control,:));    %Control Condition
    avg_test_LightResponse = mean(LFP_LightResponse(index_LR_Test,:));  %Test condition
    avg_posttest_LightResponse = mean(LFP_LightResponse(index_LR_Posttest,:));  %Posttest condition
    avg_washout_LightResponse = mean(LFP_LightResponse(index_LR_Washout,:));    %Washout condition
    
    %% Plot Figures of average waveform
    figure('name',sprintf('%s', FileName))
    subplot(1,4,1)
    plot(t(1:LightResponse_window_indices),avg_control_LightResponse, 'k')
    hold on
    plot(t(1:LightResponse_window_indices), LED(lightResponse_indices{i}), 'b')
    
    if ~isempty(InputGUI_time_spec{1}) 
        title (sprintf('%s, n=%d', InputGUI_time_spec{1}, sum(index_LR_Control)))
    else
        title (sprintf('Control Condition, n=%d', sum(index_LR_Control)))
    end
    xlabel ('Time (s)')
    ylabel ('mV')

    subplot(1,4,2)
    plot(t(1:LightResponse_window_indices),avg_test_LightResponse, 'k')
    hold on
    plot(t(1:LightResponse_window_indices), LED(lightResponse_indices{i}), 'b')
    
    if ~isempty(InputGUI_time_spec{4})    
        title (sprintf('%s, n=%d', InputGUI_time_spec{4}, sum(index_LR_Test)))
    else
        title (sprintf('HEPES-buffered ACSF, n=%d', sum(index_LR_Test)))
    end
    xlabel ('Time (s)')
    ylabel ('mV')

    subplot(1,4,3)
    plot(t(1:LightResponse_window_indices),avg_posttest_LightResponse, 'k')
    hold on
    plot(t(1:LightResponse_window_indices), LED(lightResponse_indices{i}), 'b')
    
    if ~isempty(InputGUI_time_spec{7})
        title (sprintf('%s, n=%d', InputGUI_time_spec{7}, sum(index_LR_Posttest)))
    else
        title (sprintf('HEPES-buffered ACSF + BUM [50 uM], n=%d', sum(index_LR_Posttest)))
    end
    xlabel ('Time (s)')
    ylabel ('mV')

    %Optional: Washout out
    if ~isempty(InputGUI_time_spec{10})
        subplot(1,4,4)
        plot(t(1:LightResponse_window_indices),avg_washout_LightResponse, 'k')
        hold on
        plot(t(1:LightResponse_window_indices), LED(lightResponse_indices{i}), 'b')
        title (sprintf('%s, n=%d', InputGUI_time_spec{10}, sum(index_LR_Washout)))
        xlabel ('Time (s)')
        ylabel ('mV')
    end
    
    linkaxes;   %all same y-axis
    