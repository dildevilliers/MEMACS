% Test script for the FarField.shift function
close all
clear all

FF0 = FarField.readCSTffs('CircWG_origin');
FFd = FarField.readCSTffs('CircWG_shift');

delta = pnt3D(25e-3,-50e-3,125e-3);
FF0d = FF0.shift(delta);

[r1,r2,r3] = rms(FF0d - FFd)

