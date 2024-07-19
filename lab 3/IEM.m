%% Exercise 1, Lab 3: Improved Euler Method, ODE Solver

% t0 is the starting point of the interval
% tN is the ending point of the interval
% y0 is the initial condition of the ode
% h is the stepsize of the ode
% f is the ode given to solve

function [T, Y] = IEM(t0,tN,y0,h,f)

N = round((tN-t0)/h, 0); %take the space and divide by step size to produce N intervals, 0 means round to nearest integer
T = linspace(t0, tN, N); %make even spaces of size t based on first, last, and N interval
Y = zeros(1, N); %make a row matrix of zeros that is size N to store the solutions at each X value
Y(1) = y0; %set the first position of Y to the initial value

for i = 1:(N-1) %you start at 1 so stop one early
    LS = f(T(i), Y(i));
    RS = f( T(i) + h , Y(i) + h*LS );
    Y(i+1) = Y(i) + h*(0.5 * (LS+RS));
end
