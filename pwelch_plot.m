% Choose FFT size and calculate spectrum
function FREQ_ESTIMATE = pwelch_plot(xx,fsamp)

    xx = detrend(xx);
    
    p = nextpow2(numel(xx));
    Nfft = 2^p;

    padding = Nfft - numel(xx);
    padding_zeros = zeros (padding, 1);

    xx = [xx; padding_zeros];   %xx is now padded with zeros


    %pwelch function
    [Pxx,f] = pwelch(xx);
    
    %plot time signal samples
    subplot(2,1,1);
    plot(xx);
    ylim([-0.5 0.5]);
    title('Time Signal Data');

    % Plot frequency spectrum
    subplot(2,1,2);
    plot (f*fsamp/2/pi,Pxx)
    %loglog(f,Pxx);
    ylabel('PSD'); xlabel('Frequency (Hz)');
    grid on;
    xlim ([0 100])

    % Get frequency estimate (spectral peak)
    [~,loc] = max(Pxx);
    FREQ_ESTIMATE = f(loc)*fsamp/2/pi;
    title(['Frequency estimate = ',num2str(FREQ_ESTIMATE),' Hz']);
    
    pause(0.01)
   
end