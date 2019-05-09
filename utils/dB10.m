function dB = dB10(lin)

% function dB = dB10(lin) computes the power based dB value of a linear
% input vector 'lin'

dB = 10*log10(abs(lin));