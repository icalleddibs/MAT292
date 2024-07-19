%% Exercise 6, Lab 5: Second-Order ODE Solver

% t0 is the starting point of the interval
% tN is the ending point of the interval
% y0 is the initial condition of y(t)
% y1 is the initial condition of y'(t)
% h is the stepsize of the ode
% consider y'' + p(t)y' + q(t)y = g(t)
% where we pass f = y'' = -p(t)y' - q(t)y + g(t)

function [T, Y] = DE2(t0, tN, y0, y1, h, f)

N = round((tN-t0)/h, 0); %total space divided by step size to produce N intervals, 0 means round to nearest integer
T = linspace(t0, tN, N); %matrix of even spaces
Y = zeros(1, N); %should be the y solutions

Y(1) = y0;
Y(2)= y0 + (y1*h);

for i = 1:(N-2) %we are setting Y(i+2) so stop at N-2

    y_prime = (Y(i+1) - Y(i)) / h; %since we defined Y(2) and Y(1)
    y_dbprime = f(T(i+2), y_prime, Y(i+1));
    % f is already in form to solve y''
    % give it the next time step, new y_prime, and current y value
    % we set i+1 as the current because we start with i=1 and the first
    % position is at Y(1) since we cannot index into Y(0) in matlab array

    % we are given y'' = (y(i+2) - 2y(i+1) + y(i) / h^2
    % we are solving for y(i+2) which is the next step, so rearrange
    Y(i+2) = (y_dbprime * h^2) + 2*Y(i+1) - Y(i);
    
end
