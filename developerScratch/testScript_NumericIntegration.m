% Script to test the numeric integration functions

close all
clear all

%% 1D tests
Nx = 101;
x = linspace(0,pi,Nx);
f = @(x) sin(x);
y = f(x);
It = integral1D(x,y,'trap');
Is = integral1D(x,y,'simp');
Iv = integral(f,min(x),max(x));

err_trap1D = abs(It - Iv)./abs(Iv)
err_simp1D = abs(Is - Iv)./abs(Iv)

%% 2D tests
r = 1;
Nx = 91;
x = linspace(0,2*pi,Nx);
Ny = 45;
y = linspace(0,pi,Ny);
[X,Y] = meshgrid(x,y);
Z = r.^2.*sin(Y);
Iv = 4*pi*r.^2;
It = integral2D(X,Y,Z,'trap');
Is = integral2D(X,Y,Z,'simp');

err_trap2D = abs(It - Iv)./abs(Iv)
err_simp2D = abs(Is - Iv)./abs(Iv)

