function [chord] = Bennett_Murphy_Courtright(file)
    close all
    clc

    Fs=16000;
    [x,Fs]=audioread(file);
    x1=x.';
    x=x1;
    sound(x,Fs)
    
    tic
    % Pad zeros to size 2^15
    x = horzcat(x, zeros(1,16768));

    % Calculate the needed time and frequency vectors
    t = (0:(length(x)-1)) / Fs;

    %Implement moving average, M = 16
    decimator = horzcat(ones(1,16),zeros(1,length(x)-16));
    y = conv(x, decimator);
    t2 = (0:(length(y)-1)) / Fs;

    %Downsample, M = 16
    z = [];
    for i = 1:length(y)

        if mod(i,16) == 0
            z = horzcat(z,y(i));
        end
    end
    Fs = Fs /16;
    t3 = (0:(length(z)-1)) / Fs;
    numPts = length(z);
    
    % Do the Fourier Transform
    zFT = fft(z)/numPts;
    zFT_s = fftshift(zFT);

    % The All-Pole Method
    
    % Plotting the Pole Error Plot
    %{
    [err, p] = all_pole_error(z, Fs, 100);
    figure
    plot(p, err)
    xlabel('Number of Poles')
    ylabel('Error')
    title('All-Pole Method: Number of Poles vs Error')
    grid on
    %}
    
    p = 40;
    
    [a_k] = all_pole(z, p);
    [H_sqr, w] = spec_est(a_k, length(z));
    
    figure
    subplot(2,1,1)
    plot(linspace(0, Fs/2, floor(length(zFT_s)/2)), (abs(zFT_s(length(zFT_s)/2:end-1))/10))
    title('Discrete Fourier Transform Frequency Spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (V/Hz)')
    grid on
    subplot(2,1,2)
    plot(linspace(0, Fs/2, floor(length(H_sqr)/2)), abs(H_sqr(1:2:end-1)))
    title('All-Pole Frequency Spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (V/Hz)')
    grid on
    
    %
    %
    % Start of Chord Recognition Algorithm
    %
    %
    
    %convert to positive only indices
    half_length = round(length(zFT_s)/2);
    zPos = zeros(1,half_length);
    for i = 1:length(zPos)
        zPos(i) = zFT_s(i+half_length-1);
    end

    %create magnitude array for convenience
    zPosAbs = abs(zPos);
    %limit to less than 400 Hz - chosen from reference values given in table 1
    %of project assignment
    upper = length(zPosAbs);
    lower = round(upper * 4 / 5);
    for i = lower:upper
        zPosAbs(i) = 0;
    end    

    upper = round(length(zPosAbs)*.16);
    lower = 1;
    for i = lower:upper
        zPosAbs(i) = 0;
    end 
    
    %graph positive part details
    numPtsPos = length(zPosAbs);

    %Loop for simplifying findpeaks results - combines nearby peaks in
    %progressively smaller bins to smooth FFT output. Stops looping when no
    %further combinations are being made

    %copy to preserve original array
    zPosAbsCopy = zPosAbs;
    %findpeaks
    [pks, locs] = findpeaks(zPosAbsCopy);
    %size of array counters - old version used as condition to exit while loop
    filt_length_old = 0;
    filt_length = length(locs);
    %loop findpeaks with decreasing bin sizes
    while filt_length_old ~= filt_length
        %combining nearby peaks in bins - saves largest peak per bin
        filt_length_old = filt_length;
        for i = 1:filt_length-1
            if locs(i) + filt_length/2 > locs(i+1)
                if pks(i) >= pks(i+1)
                    pks(i+1) = 0;
                else
                    pks(i) = 0;
                end
            end
        end
        %reconstructing combined peak array for next findpeaks call
        zPosAbsCopy = zeros(1,length(zPosAbs));
        for i = 1:length(locs)
            zPosAbsCopy(locs(i)) = pks(i);
        end
        %call findpeaks on newly combined array
        [pks, locs] = findpeaks(zPosAbsCopy);
        filt_length = length(locs);
    end

    %save pks just in case
    pksCopy = pks;
    %insertion sort pks array from least to greatest
    %allows for removal of straggler peaks not eliminated in above loop
    key = 0;
    j = 0;
    for i = 2:length(pksCopy)
        key = pksCopy(i);
        j = i - 1;
        while(j >= 1 && pksCopy(j) > key)
            pksCopy(j+1) = pksCopy(j);
            j = j - 1;
        end
        pksCopy(j+1) = key;
    end

    %delete all peaks besides greatest 6 for reference frequency checking
    %finds value of 6th greatest peak and removes lower values
    %final combined FFT peaks stored in pksFinal and locsFinal
    valPeak6 = pksCopy(length(pksCopy)-5);
    pksFinal = zeros(1,6);
    locsFinal = zeros(1,6);
    c = 1;
    for i = 1:length(pks)
        if pks(i) >= valPeak6
            pksFinal(c) = pks(i);
            locsFinal(c) = locs(i);
            c = c + 1;
        end
    end
    %reconstruction of FFT after elmination down to 6 peaks
    zPosAbsCopy = zeros(1,length(zPosAbs));
    for i = 1:length(locsFinal)
        zPosAbsCopy(locsFinal(i)) = pksFinal(i);
    end

    %chord recognition part
    %input - processed signal - using final FFT of zFT_s here - frequency of
    %each of the 6 peaks stored in freqLocs array
    freqLocs = locsFinal*(Fs/(2*length(zPosAbs)));
    %reference arrays for peak positions
    cRef = [98 130.8 164.8 196 261.6 329.6];
    fRef = [87.31 130.8 174.6 220 261.6 349.2];
    dRef = [92.5 110 146.8 220 293.6 270];
    gRef = [98 123.5 146.8 196 246.9 392];
    %create 4 delta value trackers for each note
    closeC = 0;
    closeF = 0;
    closeD = 0;
    closeG = 0;
    %range for either side of ideal frequency vale - can be within +- delta
    delta = length(zPosAbs)*2/Fs;
    %for each reference peak record difference between processed signal and
    %reference frequency
    for x = 1:6
        %number of frequency values that are within delta number of Hz to ideal
        %frequency value given in reference arrays
        if abs(freqLocs(x) - cRef(x)) <= delta
            closeC = closeC + 1;
        end
        if abs(freqLocs(x) - fRef(x)) <= delta
            closeF = closeF + 1;
        end
        if abs(freqLocs(x) - dRef(x)) <= delta
            closeD = closeD + 1;
        end
        if abs(freqLocs(x) - gRef(x)) <= delta
            closeG = closeG + 1;
        end 
    end

    %selection value to keep track of min delta difference, starts at 1
    %1-C
    %2-F
    %3-D
    %4-G
    %currMax - tracks note letter with maximum hits within delta ranges of the
    %ideal frequency
    %starts equal to totDifC value

    selection = 1;
    currMax = closeC;
    if closeF > currMax
        selection = 2;
        currMax = closeF;
    end
    if closeD > currMax
        selection = 3;
        currMax = closeD;
    end
    if closeG > currMax
        selection = 4;
        currMax = closeG;
    end
    %outprints correct note based on selection value
    if selection == 1
        chord='C'
    elseif selection == 2
        chord='F'
    elseif selection == 3
        chord='D'
    else
        chord='G'
    end

    %
    %
    % End Chord Recogntion Algorithm
    %
    %
    toc
end

% Returns Normalized Pole vs Error Plot
function [err, p] = all_pole_error(x, Fs, max_poles)
    p = 1:max_poles;
    err = zeros(1, max_poles);
    
    % Correct Method for Pole Error Plot: E{(x[n]-sum(ak*x[n-k])^2}
    %{
    for i=p
        w_n = wgn(1, length(x), 1);
        [a] = all_pole(x, i);
        tmp = zeros(i, length(x));
        for j=1:i
            tmp(j,j:end) = x(1:end-j+1);
        end
        err(i) = mean((x-sum(a.*tmp)).^2);
    end
    %}
    
    % Error vs Poles calcuated by comparing spectrums
    for i=p
        [a] = all_pole(x, i);
        [H_sqr, w] = spec_est(a, length(x));
        x_f = fftshift(fft(x));
        err(i) = mean(abs(x_f(length(x)/2:end-1))-(abs(H_sqr(1:2:end-1))/100));
    end
    err = err/max(err); 
end

%Spectrum Estimation for the All-Pole Method
function [H_sqr, w] = spec_est(a, num_pts)
    w=linspace(0, pi, num_pts);
    i=1:num_pts;
    k=1:length(a);
    H_sqr = 1./(1-a(k)'*exp(-j*k'*w(i)));
end

%Generate All-Pole Constants a1, a2, ..., ap
function [a] = all_pole(x, p)
    [~, M] = corrmtx(x,p);
    R=M(1:p, 1:p);
    phi=M(2:end,1);
    a=inv(R)*phi;
end
