clear
clc

[a_c, a_d, a_f, a_g] = initialize_filters(40);

Fs = 16000;
file = 'c_1.wav';
[x, Fs] =audioread(file);
sound(x,Fs) %Sends audio to speakers

x = x.';

% Pad zeros to size 2^15
%x = horzcat(x, zeros(1,16768));
x1 = zeros(1, 16768);
x1(1:length(x)) = x;
x = x1;


% Calculate the needed time and frequency vectors
t = (0:(length(x)-1)) / Fs;

%Implement moving average, M = 16
decimator = horzcat(ones(1,16),zeros(1,length(x)-16));
y = conv(x, decimator);
t2 = (0:(length(y)-1)) / Fs;

%Downsample, M = 16
z = y(1:16:end);

Fs = Fs /16;
t3 = (0:(length(z)-1)) / Fs;
numPts = length(z);

test_c = spec_est(a_c, z);
test_d = spec_est(a_d, z);
test_f = spec_est(a_f, z);
test_g = spec_est(a_g, z);

figure
subplot(4,1,1)
plot(1:length(test_c), abs(test_c))
subplot(4,1,2)
plot(1:length(test_d), abs(test_d))
subplot(4,1,3)
plot(1:length(test_f), abs(test_f))
subplot(4,1,4)
plot(1:length(test_g), abs(test_g))



function [H_sqr] = spec_est(a, x)
    H_sqr = zeros(1, length(x));
    i=1;
    for w=linspace(0, pi, length(x))
        tmp = zeros(1, length(a));
        for k=1:length(a)
            tmp(k) = a(k)*exp(-j*w*k);
        end
        H_sqr(i) = 1/(1-sum(tmp));
        i=i+1;
    end
end

function [a_c_chord, a_d_chord, a_f_chord, a_g_chord] = initialize_filters(p)
    % C Chord
    Fs = 16000;
    file = 'c_1.wav';
    [x, Fs] =audioread(file);
    x = x.';

    % Pad zeros to size 2^15
    x1 = zeros(1, 16768);
    x1(1:length(x)) = x;
    x = x1;

    %Implement moving average, M = 16
    decimator = horzcat(ones(1,16),zeros(1,length(x)-16));
    y = conv(x, decimator);

    %Downsample, M = 16
    z = y(1:16:end);

    [~, M] = corrmtx(z,p);
    R=M(1:p, 1:p);
    phi=M(2:end,1);
    a_c_chord=inv(R)*phi;
    
    
    % D Chord
    Fs = 16000;
    file = 'd_1.wav';
    [x, Fs] =audioread(file);
    x = x.';

    % Pad zeros to size 2^15
    x1 = zeros(1, 16768);
    x1(1:length(x)) = x;
    x = x1;

    %Implement moving average, M = 16
    decimator = horzcat(ones(1,16),zeros(1,length(x)-16));
    y = conv(x, decimator);

    %Downsample, M = 16
    z = y(1:16:end);

    [~, M] = corrmtx(z,p);
    R=M(1:p, 1:p);
    phi=M(2:end,1);
    a_d_chord=inv(R)*phi;
    
    
    % F Chord
    Fs = 16000;
    file = 'f_1.wav';
    [x, Fs] =audioread(file);
    x = x.';

    % Pad zeros to size 2^15
    x1 = zeros(1, 16768);
    x1(1:length(x)) = x;
    x = x1;

    %Implement moving average, M = 16
    decimator = horzcat(ones(1,16),zeros(1,length(x)-16));
    y = conv(x, decimator);

    %Downsample, M = 16
    z = y(1:16:end);

    [~, M] = corrmtx(z,p);
    R=M(1:p, 1:p);
    phi=M(2:end,1);
    a_f_chord=inv(R)*phi;
    
    
    % G Chord
    Fs = 16000;
    file = 'g_1.wav';
    [x, Fs] =audioread(file);
    x = x.';

    % Pad zeros to size 2^15
    x1 = zeros(1, 16768);
    x1(1:length(x)) = x;
    x = x1;

    %Implement moving average, M = 16
    decimator = horzcat(ones(1,16),zeros(1,length(x)-16));
    y = conv(x, decimator);

    %Downsample, M = 16
    z = y(1:16:end);

    [~, M] = corrmtx(z,p);
    R=M(1:p, 1:p);
    phi=M(2:end,1);
    a_g_chord=inv(R)*phi;
end