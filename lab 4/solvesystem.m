%% Exercise 1, Lab 4: Improved Euler Method, ODE Solver

% t0 is the starting point of the interval
% tN is the ending point of the interval
% x0 is a vector for the initial condition, where x(t0)=x0
% h is the stepsize of the ode
% f and g are the odes given to solve

function [T, X] = solvesystem(t0, tN, x0, h, f, g)

N = round((tN-t0)/h, 0); %total space divided by step size to produce N intervals, 0 means round to nearest integer
T = linspace(t0, tN, N); %matrix of even spaces
X = zeros(2, N); %2 rows corresponding to the 2 odes (row 1 for x1', row 2 for x2' solutions)
X(1,1) = x0(1,1); %set the first position of row 1 to initial 1
X(2,1) = x0(2,1); %set the first position of row 2 to initial 2

for i = 1:(N-1) %you start at 1 so stop one early
    % left side point is f(T, x1, x2) at i
    x1_LS = f( T(i), X(1,i), X(2,i) );
    x2_LS = g( T(i), X(1,i), X(2,i) );

    % right side point is next time step and current point + h*LS (for both
    % x1 and x2 values this time)
    x1_RS = f( T(i)+h, X(1,i)+(h*x1_LS), X(2,i)+(h*x2_LS));
    x2_RS = g( T(i)+h, X(1,i)+(h*x1_LS), X(2,i)+(h*x2_LS));
    
    % next value x(i+1) is x(i) + h* (average of LS & RS which is LS+RS/2)
    X(1,i+1) = X(1,i) + h*(0.5 * (x1_LS+x1_RS));
    X(2,i+1) = X(2,i) + h*(0.5 * (x2_LS+x2_RS));

    % essentially the same thing as IEM from lab 3 but increment both values
end
