clear all
close all
clc

Fs = 16000;
file = 'f_1.wav';
[x,Fs] =audioread(file);
%sound(x,Fs)
x1=x.';
t = (0:(length(x1)-1)) / Fs;
numPts = 16000;
f = (-numPts/2 : numPts/2-1)*Fs/numPts;
xFT = fft(x1)/numPts;
xFT_s = fftshift(xFT);

figure
subplot(2,1,1)
plot(t, x1, 'Linewidth', 1.5)
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