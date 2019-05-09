function dB = dB20(lin)

% function dB = dB20(lin) computes the voltage/current based dB value of a 
% linear input vector 'lin'

dB = 20*log10(abs(lin));