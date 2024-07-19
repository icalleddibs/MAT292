%% Systems Lab: Systems of ODEs in MATLAB
%
% In this lab, you will write your own ODE system solver for the Heun method 
% (aka the Improved Euler method), and compare its results to those of |ode45|.
%
% You will also learn how to save images in MATLAB.
% 
% Opening the m-file lab4.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are four (4) exercises in this lab that are to be handed in on the
% due date of the lab.  Write your solutions in a separate file, including
% appropriate descriptions in each step. Save the m-files and the pdf-file 
% for Exercise 4 and submit them on Quercus.

%% Exercise 1
%
% Objective: Write your own ODE system solver using the Heun/Improved Euler 
% Method and compare it to |ode45|.
%
% Details: Consider the system of 2 ODEs:
%
% |x1'=f(t,x1,x2), x2'=g(t,x1,x2)|
%
% This m-file should be a function which accepts as variables 
% (t0,tN,x0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, h is the stepsize, and x0 is a vector 
% for the initial condition of the system of ODEs |x(t0)=x0|.
% Name the function solvesystem_<UTORid>.m (Substitute your UTORid for <UTORid>).
% You may also want to pass the functions into the ODE the way |ode45| 
% does (check MATLAB labs 2 and 3).
%
% Your m-file should return a row vector of times and a matrix of 
% approximate solution values (the first row has the approximation for |x1|
% and the second row has the approximation for |x2|).
% 
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

%% Exercise 2
%
% Objective: Compare Heun with an exact solution
%
% Details:  Consider the system of ODEs
%
% |x1' = x1/2 - 2*x2, x2' = 5*x1 - x2|
%
% with initial condition |x(0)=(1,1)|.
%
% Use your method from Exercise 1 to approximate the solution from |t=0|
% to |t=4*pi| with step size |h=0.05|.
%
% Compute the exact solution (by hand) and plot both phase portraits on 
% the same figure for comparison.
%
% Your submission should show the construction of the inline function, the
% use of your Heun's method to obtain the solution, a construction of the 
% exact solution, and a plot showing both.  In the comments, include the 
% exact solution.
%
% Label your axes and include a legend.

f1 = @(t,x1,x2) x1/2 - 2*x2;
g1 = @(t,x1,x2) 5*x1 - x2;

t01 = 0;
tN1 = 4*pi;
h1 = 0.05;

x01 = 2:1; %two rows, 1 column with values of 1 and 1
x01(1,1) = 1;
x01(2,1) = 1;

% IEM solution
[ta, ya] = solvesystem(t01, tN1, x01, h1, f1, g1);

% exact solution
t = linspace(0, 4*pi, 4*pi/0.05);

eg1 = (-0.25 - (sqrt(151)*1i)/4);
eg2 = (-0.25 + (sqrt(151)*1i)/4);
tens = ((10 * sqrt(151) * 1i)/151);
v1 = ((1/20)*(3-(sqrt(151)*1i)));
v2 = ((1/20)*(3+(sqrt(151)*1i)));
inv1 = (0.5 - (3*(sqrt(151)*1i/302)));
inv2 = (0.5 + (3*(sqrt(151)*1i/302)));

x1ex = ((exp(eg1*t))*((v1*tens) + (v1*inv1)) + (exp(eg2*t))*((v2*tens*(-1)) + (v2*inv2)));
x2ex = (exp(eg1*t))*(tens + inv1) + (exp(eg2*t))*((tens * (-1)) + inv2);

%this was solved by hand using matrix exponential, as taught in class.
%i am aware there is a solution using cos and sin which is shorter, but i
%thought i should do it as taught in class

%here is the whole solution written out instead of using placeholder
%variables:
%x1 = ((exp((-0.25 - (sqrt(151)*1i)/4)*t))*(((1/20)*(3-(sqrt(151)*1i))*
%((10 * sqrt(151) * 1i)/151)) + (((1/20)*(3-(sqrt(151)*1i)))*(0.5 - (3*(sqrt
%(151)*1i/302))))) + (exp((-0.25 + (sqrt(151)*1i)/4)*t))*((((1/20)*(3+(sqrt
%(151)*1i)))*((10 * sqrt(151) * 1i)/151)*(-1)) + (((1/20)*(3+(sqrt(151)*1i)
% ))*(0.5 + (3*(sqrt(151)*1i/302))))))

%x2 = (exp((-0.25 - (sqrt(151)*1i)/4)*t))*(((10 * sqrt(151) * 1i)/151) +
%(0.5 - (3*(sqrt(151)*1i/302)))) + (exp((-0.25 + (sqrt(151)*1i)/4)*t))*(((
%(10 * sqrt(151) * 1i)/151) * (-1)) + (0.5 + (3*(sqrt(151)*1i/302))))

%plot x1 against x2 instead of against time
plot(ya(1,:), ya(2,:), '+',  x1ex, x2ex, 'LineWidth', 1);
title({"IEM and Exact Solution of", "x1' = x1/2 - 2*x2, x2' = 5*x1 - x2, x(0)=(1,1)"});
legend("IEM Solution", "Exact Solution");
xlabel("x1");
ylabel("x2");

%% Exercise 3
%
% Objective: Compare your method with Euler's Method (from |iode|).
%
% Details: Use |iode| to plot the solution for the same problem with the same
% step size as on Exercise 2.
%
% Compare your solution on exercise 2, the exact solution from exercise 2
% and the approximation using Euler's method. Plot the solution for 
% Euler's method and make note of any differences.

t2 = linspace(0, 4*pi, 4*pi/0.05);
f2 = @(t,x) [x(1)/2 - 2*x(2); 5*x(1) - x(2)]; %make a matrix of the equations
x_euler = euler(f2, x01, t2); %2x1 matrix of solutions

plot(x_euler(1,:), x_euler(2,:), ya(1,:), ya(2,:), '+', x1ex, x2ex, 'LineWidth', 1);
title({"Euler, IEM, Exact Solutions of", " x1' = x1/2 - 2*x2, x2' = 5*x1 - x2, x(0)=(1,1)"});
legend("Euler Solution", "IEM Solution", "Exact Solution");
xlabel("x1");
ylabel("x2");

% The approximation using Euler's method for the same amount of points (and
% step size) as IEM is really bad. It is spiraling at a much slower rate 
% than the exact solution and IEM, as we can see from how it completes 6
% rotations but has barely converged towards the middle, versus how IEM and
% exact solution get very close to the middle (equilibrium point) after 6
% rotations. When I tested by decreasing the step size, I got much better
% results. This difference likely comes from using EM vs. IEM. We learned
% in class that the approximation for IEM produces much better solutions
% because of the steps it takes (trapezoid instead of rectangular). As
% shown before, the IEM and exact solution are extremely similar.

%% Saving Images in MATLAB
%
% To do the following exercises, you will need to know how to output 
% graphics from MATLAB.  Create a folder on your Desktop (or elsewhere) to 
% contain the files generated by these exercises.  Make this folder the 
% "Current Folder" in the left side of the main MATLAB window.  This will 
% ensure that the files output by MATLAB end up in the folder you created.
%
% To save an image of a phase portrait, use the following steps:
%
% 1. Get the phase portrait looking the way you want in the |iode| window.  
%
% 2. Leaving |iode| open, switch to the main MATLAB window.
%
% 3. Type the command |print -dpng -r300 'filename.png'| in the command 
% window.
%
% This command will create a PNG graphic called |filename.png| in the 
% current folder.  The |-dpng| option tells MATLAB to output the graphic in
% PNG format; MATLAB also allows output in other formats, such as BMP, EPS,
% PNG and SVG.  The |-r300| option tells MATLAB to set the resolution at 
% 300 dots per inch and can be adjusted if you wish.


%% Exercise 4
%
% Objective: Analyze phase portraits.
%
% Details: Compile the results of the following exercises into a single
% document (e.g. using a word processor) and export it to |PDF| for
% submission on Quercus. 
%
% For each of the first-order systems of ODEs 4.1 to 4.10 below, do the
% following exercises:
%
% (a) Generate a phase portrait for the system (centre the graph on the
% equilibrium point at (0,0)). Include a few trajectories.
%
% (b) Classify the equilibrium on asymptotic stability, and behaviour 
% (sink, source, saddle-point, spiral, center, proper node, improper node) 
% - check table 3.5.1 and figure 3.5.7.
% Classify also as for clockwise or counterclockwise movement, when relevant.
%
% (c) Compute the eigenvalues of the matrix (you do not need to show your
% calculations). Using the eigenvalues you computed, justify part (b).
%
% To avoid numerical error, you should use Runge-Kutta solver with a step
% size of |0.05|. Change the display parameters, if necessary, to best
% understand the phase portrait.
%
% ANSWERS INCLUDED ON THE PDF!
%
% 4.1. |dx/dt = [2 1; 1 3] x|
% b)i) asymptotic stability: unstable
% b)ii) behaviour: source
% b)iii) movement: N/A
% c) eigenvalue 1: (5+sqrt(5))/2
%    eigenvalue 2: (5-sqrt(5))/2
% justification: the eigenvalues are positive, real, and not the same, so
% there is an unstable node (source).
%
% 4.2. |dx/dt = [-2 -1; -1 -3] x|
% b)i) asymptotic stability: asymptotically stable
% b)ii) behaviour: sink
% b)iii) movement: N/A
% c) eigenvalue 1: (-5+sqrt(5))/2
%    eigenvalue 2: (-5-sqrt(5))/2
% justification: the eigenvalues negative, real, and not the same, so
% there is a stable node (sink).
%
% 4.3. |dx/dt = [-4 -6; 3 5] x|
% b)i) asymptotic stability: unstable
% b)ii) behaviour: saddle-point
% b)iii) movement: N/A
% c) eigenvalue 1: 2
%    eigenvalue 2: -1
% justification: one eigenvalue is positive, while the other is negative,
% both real, so there is an unstable saddle point.
%
% 4.4. |dx/dt = [4  6; -3 -5] x|
% b)i) asymptotic stability: unstable
% b)ii) behaviour: saddle-point
% b)iii) movement: N/A
% c) eigenvalue 1: 1 
%    eigenvalue 2: -2
% justification: one eigenvalue is positive, while the other is negative,
% both real, so there is an unstable saddle point.
%
% 4.5. |dx/dt = [0  -1; 1 -1] x|
% b)i) asymptotic stability: asymptotically stable
% b)ii) behaviour: spiral
% b)iii) movement: CCW
% c) eigenvalue 1: -1/2 + sqrt(3)i/2
%    eigenvalue 2: -1/2 - sqrt(3)i/2
% justification: the eigenvalues are complex, with negative real parts, so
% there is a stable inward spiral point
%
% 4.6. |dx/dt = [0  1; -1 1] x|
% b)i) asymptotic stability: unstable
% b)ii) behaviour: spiral
% b)iii) movement: CW
% c) eigenvalue 1: 1/2 + sqrt(3)i/2
%    eigenvalue 2: 1/2 - sqrt(3)i/2
% justification: the eigenvalues are complex, with positive real parts, so
% there is an unstable outward spiral point
%
% 4.7. |dx/dt = [2  8; -1 -2] x|
% b)i) asymptotic stability: stable
% b)ii) behaviour: centre
% b)iii) movement: CW
% c) eigenvalue 1: 2i
%    eigenvalue 2: -2i
% justification: the eigenvalues are complex, with no real parts, only
% imaginary parts with opposite signs, so it's a stable centre
%
% 4.8. |dx/dt = [-2 -8; 1 2] x|
% b)i) asymptotic stability: stable
% b)ii) behaviour: centre
% b)iii) movement: CCW
% c) eigenvalue 1: 2i
%    eigenvalue 2: -2i
% justification: the eigenvalues are complex, with no real parts, only
% imaginary parts with opposite signs, so it's a stable centre
%
% 4.9. |dx/dt = [-8 5; -13 8] x|
% b)i) asymptotic stability: stable
% b)ii) behaviour: centre
% b)iii) movement: CW
% c) eigenvalue 1: i
%    eigenvalue 2: -i
% justification: the eigenvalues are complex, with no real parts, only
% imaginary parts with opposite signs, so it's a stable centre
%
% 4.10. |dx/dt = [8 -5; 13 -8] x|
% b)i) asymptotic stability: stable
% b)ii) behaviour: centre
% b)iii) movement: CCW
% c) eigenvalue 1: i
%    eigenvalue 2: -i
% justification: the eigenvalues are complex, with no real parts, only
% imaginary parts with opposite signs, so it's a stable centre
%
%
%
