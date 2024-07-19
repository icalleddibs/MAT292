%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.
%
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|

f_a = @(t,y) y*(tan(t)) + sin(t);
t0a = 0;
tNa = pi;
y0a = -0.5;
ha = 0.001;

soln_a = ode45(f_a, [t0a, tNa], y0a);

[ta, ya] = IEM(t0a, tNa, y0a, ha, f_a);
plot(soln_a.x, soln_a.y, ta, ya);

title('Graph of dy/dt = ytan(t) + sin(t), y(0) = -1/2');
xlabel('t');
ylabel('y');
legend('Numerical Solution', 'IEM Solution', 'Location','Best');

%%
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|

f_b = @(t,y) 1 / (y^2);
t0b = 1;
tNb = 10;
y0b = 1;
hb = 0.001;

soln_b = ode45(f_b, [t0b, tNb], y0b);

[tb, yb] = IEM(t0b, tNb, y0b, hb, f_b);
plot(soln_b.x, soln_b.y, tb, yb);

title('Graph of dy/dt = 1 / y^2 , y(1) = 1');
xlabel('t');
ylabel('y');
legend('Numerical Solution', 'IEM Solution', 'Location','Best');

%%
% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|

f_c = @(t,y) (1 - t*y)/2;
t0c = 0;
tNc = 10;
y0c = -1;
hc = 0.001;

soln_c = ode45(f_c, [t0c, tNc], y0c);

[tc, yc] = IEM(t0c, tNc, y0c, hc, f_c);
plot(soln_c.x, soln_c.y, tc, yc);

title('Graph of dy/dt = (1-ty)/2 , y(0) = -1');
xlabel('t');
ylabel('y');
legend('Numerical Solution', 'IEM Solution', 'Location','Best');

%%
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|

f_d = @(t,y) y^3 - t^2;
t0d = 0;
tNd = 1;
y0d = 1;
hd = 0.001;

soln_d = ode45(f_d, [t0d, tNd], y0d);

[td, yd] = IEM(t0d, tNd, y0d, hd, f_d);
plot(soln_d.x, soln_d.y, td, yd);

title('Graph of dy/dt = y^3 - t^2 , y(0) = 1');
xlabel('t');
ylabel('y');
legend('Numerical Solution', 'IEM Solution', 'Location','northwest');


% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.

% Part a has strange behaviour at approximately t= pi/2. At this point,
% from plotting in IODE, we can see that there is a vertial asymptote. The
% graph suddenly curves downwards to t = pi/2, shoots almost vertically
% upwards, and back down again before curving and rejoining the original
% path that we can see from the ode45 approximation

% Part b and c have pretty good approximations that follow the same path as
% the numerical solution from ode45

% Part d also has a pretty good solution that is essentially the same as
% the numerical solution. However, the solution itself has a vertical
% asymptote at t=5.066046e-01, which is what the error tells us. The
% function extends to infinity at this point.

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
% see euler.m file within iode files, call the euler function
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.5|.

f3 = @(t,y) 2*t*sqrt(1-(y^2));
t3 = linspace(0, 0.5);
y3 = euler(f3, 0, t3);

% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
%
% The general solution is y = sin(t^2) + C where at y(0) = 0 then c = 0
% The particular solution is y = sin(t^2)
% The solution at t=0.5 is y(0.5) = sin(0.5^2) = sin(0.25) = 0.2474

texact = linspace(0, 0.5);
yexact = sin(texact.^2);
plot(t3, y3, texact, yexact);
title("dy/dx = 2*t*sqrt(1-y^2 ), y(0) = 0");
ylabel("y");
xlabel("t");
legend('Euler Solution', 'Exact Solution', 'Location','northwest');

% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.
%
% The bound for En = (1/2 * (1 + M) * dt) * (exp(M*dt*n) - 1)
% M is the maximum of the function which always bounds the function f
% n is the amount of steps
% dt is the step size for t (our x-axis)
%
% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.

nd = 100;
dtd = (0.5/nd);
Md = 2; %from the maximum of f, df/dt, df/dy

% the actual error bound is:
En_d = (1/2 * (1 + Md) * dtd) * (exp(Md*dtd*nd) - 1);  %error: 0.0129
fprintf('Computed Error: %g\n', En_d);

% the actual error = abs(yn - y(tn)) so true value at n - calculated at n:
AE_d = abs(sin(0.5^2) - y3(nd)); %actual error: 0.0024
fprintf('Actual Error: %g\n', AE_d);

% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.

ne = 1000; 
dte = (0.5/ne); %reducing step size by a factor of 10
Me = 2;
t3e = linspace(0, 0.5, 1000);
y3e = euler(f3, 0, t3e);

% the actual error bound is:
En_e = (1/2 * (1 + Me) * dte) * (exp(Me*dte*ne) - 1); % error: 0.0013
fprintf('Computed Error: %g\n', En_e);

% the actual error = abs(yn - y(tn)) so true value at n - calculated at n:
AE_e = abs(sin(0.5^2) - y3e(ne)); %actual error: 0.00023636
fprintf('Actual Error: %g\n', AE_e);

% 0.0024 / 10 = 0.00024, consistent with the change in step size.
% the actual error & actual error bound is reduced by a factor of 10 
% when we reduce the step size by 10. it is proportional to the step size value.
% therefore it confirms the first order of euler's method.


%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.

% The purpose of the update formula is to reduce the step size of the
% original step size given to make sure the error is not larger than to1
% (1e-8). It will take 90% of the original step size in combination with
% the ratio between our set error(to1) and the calculated difference (D).
% If D is much larger than the error, then we will be reducing by a huge
% amount to try and ensure the next step is below the tolerance range.
% Otherwise, if D is only slightly larger than the error, it will only
% reduce by a little bit because the ratio will be close to 1.

%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.

f5 = @(t,y) 2*t*sqrt(1-(y^2));
t5 = linspace(0, 0.75, (0.75/0.025)); %there will be 30 steps if we use a stepsize of 0.025
y5 = euler(f5, 0, t5);

% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.

f_5 = @(t,y) 2*t*sqrt(1-(y^2));
t05 = 0;
tN5 = 0.75;
y05 = 0;
h5 = 0.025;

[AT5, AY5] = AEM(t05, tN5, y05, h5, f_5);

% (c) Plot both approximations together with the exact solution.

soln_exact_5y = sin(t5.^2);
soln_exact_5x = t5;

plot(t5, y5, AT5, AY5, soln_exact_5x, soln_exact_5y);
legend("Euler's Method Solution", 'AEM Solution', 'Exact Solution', 'Location','northwest');
title("dy/dx = 2*t*sqrt(1-y^2 ), y(0) = 0")
ylabel("y");
xlabel("t");

%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.
% 
% We know that Euler's method is further from the actual solution by
% lookign at the plot. The reason is because it calculated using less steps
% (we assigned it a step size so it followed through and only used 30
% steps). However, AEM is almost identical to the exact solution, and we
% can see from the size of the t and y matrices for AEM that they used 5647
% steps. It is much more precise because it took significantly more points
% as calculations and approximated over smaller step sizes. AEM adjusts the
% steps size, but the step size is fixed in EM. That means if there is a
% sudden large change in the slope, EM cannot accurately represent this,
% but AEM can adjust and use more points to demonstrate the change.
%
% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.

% euler approximation from iode
f6 = @(t,y) 2*t*sqrt(1-(y^2));
t6 = linspace(0, 1.5); % didn't specify step size
y6 = euler(f6, 0, t6);

% adaptive euler approximation
f_6 = @(t,y) 2*t*sqrt(1-(y^2));
t06 = 0;
tN6 = 1.5;
y06 = 0;
h6 = 0.025;

[AT6, AY6] = AEM(t06, tN6, y06, h6, f_6);

% exact solution
soln_exact_6y = sin(t6.^2);
soln_exact_6x = t6;

plot(t6, y6, AT6, AY6, soln_exact_6x, soln_exact_6y);
legend("Euler's Method Solution", 'AEM Solution', 'Exact Solution', 'Location','northwest');
title("dy/dx = 2*t*sqrt(1-y^2 ), y(0) = 0")
ylabel("y");
xlabel("t");

% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.
%
% In the command window, we are given warning: "imaginary parts of complex
% X and/or Y arguments ignored." Because our differential equation has a
% square root symbol, that is where our imaginary numbers come from when
% (y^2 > 1). This starts to happen around t = 1.25, which is also the peak
% of our exact solution (y = 1 in the sine wave, close to pi/2). Because the
% approximations cannot use imaginary numbers, the plotted solutions start
% to stray away from the exact solution, or have close to a 0 slope. They
% cannot use previous values because they cannot use imaginary numbers, so
% the next values after this point are identical or very close (approximately
% y=1, and they present as horizontal lines.
%
%
