
% file is the sound data sampled at 16kHz in the form of a .wav file
function [chord] = Bennett_Murphy_Courtright(file)

Fs = 16000;
[x,Fs] =audioread(file);
sound(x,Fs)

chord = 'We are dumb';
end

