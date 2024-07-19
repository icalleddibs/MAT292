%% Exercise 4, Lab 3: Adapted Euler Method, ODE Solver

% t0 is the starting point of the interval
% tN is the ending point of the interval
% y0 is the initial condition of the ode
% h is the stepsize of the ode
% f is the ode given to solve

function [t, y] = AEM(t0,tN,y0,h,f)

N = round((tN-t0) / h, 0) + 1;
t = zeros(1, N);
y = zeros(1, N);
y(1) = y0;
i = 2;
to1 = 1e-8;

while (t(end) < tN) %while the last value of t is still less than final time
    while true
        Y = y(i-1) + (f(t(i-1), y(i-1))) * h; % euler step size of h
        Z0 = y(i-1) + (f(t(i-1), y(i-1))) * (h/2); % first euler step size of h/2
        Z = Z0 + (f(t(i-1) + (h/2), Z0)) * (h/2); % next  of h/2 using previous step
        D = Z - Y;

        if (abs(D) < to1)
            y(i) = Z + D; %successful so add new solution value
            t(i) = t(i-1) + h; %update the time value
            i = i + 1; %iterate to next step
            break;
        else
            h = 0.9 * h * min(max(to1/abs(D), 0.3), 2); %given updated step size
        end
    end
end
