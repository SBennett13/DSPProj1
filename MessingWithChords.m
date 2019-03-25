clear all
close all
clc

Fs = 16000;
file = 'd_3.wav';
[x,Fs] =audioread(file);
x1=x.';
x=x1;
sound(x,Fs)

% Pad zeros to size 2^15
x = horzcat(x, zeros(1,16768));

% Calculate the needed time and frequency vectors
t = (0:(length(x)-1)) / Fs;
numPts = length(x);
f = (-numPts/2 : numPts/2-1)*Fs/numPts;

% Do the Fourier Transform
xFT = fft(x)/numPts;
xFT_s = fftshift(xFT);

% Plot the raw result
figure
subplot(2,1,1)
plot(t, x, 'Linewidth', 1.5)
title('Raw f Chord :: Time Domain')
ylabel('Amplitude (V)')
xlabel('Time (s)')
grid on
subplot(2,1,2)
plot(f, abs(xFT_s), 'Linewidth', 1.5)
title('Raw f Chord :: Frequency Domain')
ylabel('Amplitude (V/Hz)')
xlabel('Frequency (Hz)')
grid on

%Implement moving average, M = 16
decimator = horzcat(ones(1,16),zeros(1,length(x)-16));
y = conv(x, decimator);
t2 = (0:(length(y)-1)) / Fs;
numPts = length(y);
f2 = (-numPts/2 : numPts/2-1)*Fs/numPts;

% Do the Fourier Transform
yFT = fft(y)/numPts;
yFT_s = fftshift(yFT);

% Plot the raw result
figure
subplot(2,1,1)
plot(t2, y, 'Linewidth', 1.5)
title('CIC Filtered f Chord :: Time Domain')
ylabel('Amplitude (V)')
xlabel('Time (s)')
grid on
subplot(2,1,2)
plot(f2, abs(yFT_s), 'Linewidth', 1.5)
title('CIC Filtered f Chord :: Frequency Domain')
ylabel('Amplitude (V/Hz)')
xlabel('Frequency (Hz)')
grid on

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
f3 = (-numPts/2 : numPts/2-1)*Fs/numPts;

% Do the Fourier Transform
zFT = fft(z)/numPts;
zFT_s = fftshift(zFT);

% Plot the raw result
figure
subplot(2,1,1)
plot(t3, z, 'Linewidth', 1.5)
title('Downsampled f Chord :: Time Domain')
ylabel('Amplitude (V)')
xlabel('Time (s)')
grid on
subplot(2,1,2)
plot(f3, abs(zFT_s), 'Linewidth', 1.5)
title('Downsampled f Chord :: Frequency Domain')
ylabel('Amplitude (V/Hz)')
xlabel('Frequency (Hz)')
grid on

%{
This code does a comparison of the individual harmonics
file = 'f_1.wav';
[x,Fs] =audioread(file);
sound(x,Fs)
x1=x.';
x=x1;
figure
subplot(2,1,1)
plot(t, x, 'Linewidth', 1.5)
title("f_1 :: Time Domain")
ylabel('Amplitude (V)')
xlabel('Time (s)')
grid on
subplot(2,1,2)
plot(f, abs(xFT_s), 'Linewidth', 1.5)
title('f_1 :: Frequency Domain')
ylabel('Amplitude (V/Hz)')
xlabel('Frequency (Hz)')
grid on

%Comparing another f
file = 'f_2.wav';
[x,Fs] =audioread(file);
sound(x,Fs)
x1=x.';
x=x1;
xFT = fft(x)/numPts;
xFT_s = fftshift(xFT);

figure
subplot(2,1,1)
plot(t, x, 'Linewidth', 1.5)
title("f_1 :: Time Domain")
ylabel('Amplitude (V)')
xlabel('Time (s)')
grid on
subplot(2,1,2)
plot(f, abs(xFT_s), 'Linewidth', 1.5)
title('f_1 :: Frequency Domain')
ylabel('Amplitude (V/Hz)')
xlabel('Frequency (Hz)')
grid on

%Comparing another f
file = 'f_3.wav';
[x,Fs] =audioread(file);
sound(x,Fs)
x1=x.';
x=x1;
xFT = fft(x)/numPts;
xFT_s = fftshift(xFT);

figure
subplot(2,1,1)
plot(t, x, 'Linewidth', 1.5)
title("f_1 :: Time Domain")
ylabel('Amplitude (V)')
xlabel('Time (s)')
grid on
subplot(2,1,2)
plot(f, abs(xFT_s), 'Linewidth', 1.5)
title('f_1 :: Frequency Domain')
ylabel('Amplitude (V/Hz)')
xlabel('Frequency (Hz)')
grid on
%}