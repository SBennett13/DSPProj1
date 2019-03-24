clear all
close all
clc

% Add all harmonics to a single vector
x1=0;
Fs = 16000;
file = 'f_1.wav';
[x,Fs] =audioread(file);
x1=x1 + x;
file = 'f_2.wav';
[x,Fs] =audioread(file);
x1=x1 + x;
file = 'f_3.wav';
[x,Fs] =audioread(file);
x1=x1 + x;
sound(x1,Fs)
x=x1.';

% Calculate the needed time and frequency vectors
t = (0:(length(x)-1)) / Fs;
numPts = 16000;
f = (-numPts/2 : numPts/2-1)*Fs/numPts;

% Do the Fourier Transform
xFT = fft(x)/numPts;
xFT_s = fftshift(xFT);

% Plot the raw result
figure
subplot(2,1,1)
plot(t, x, 'Linewidth', 1.5)
title("Raw f Chord :: Time Domain")
ylabel('Amplitude (V)')
xlabel('Time (s)')
grid on
subplot(2,1,2)
plot(f, abs(xFT_s), 'Linewidth', 1.5)
title('Raw f Chord :: Frequency Domain')
ylabel('Amplitude (V/Hz)')
xlabel('Frequency (Hz)')
grid on

%{
This code does a comparison of the individual harmonics
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