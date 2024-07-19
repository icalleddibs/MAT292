%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|.
%
% There are five (5) exercises in this lab that are to be handed in.  
% Write your solutions in a separate file, including appropriate descriptions 
% in each step.

%% Using symbolic variables to define functions
% 
% In this exercise we will use symbolic variables and functions.

syms t s x y

f = cos(t)
h = exp(2*x)


%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

%%
% By default it uses the variable |s| for the Laplace transform
% But we can specify which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

%% 
% We can also specify which variable to use to compute the Laplace
% transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

%% 
% We can also use inline functions with |laplace|. When using inline
% functions, we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

%% 
% MATLAB also has the routine |ilaplace| to compute the inverse Laplace
% transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

%% 
% If |laplace| cannot compute the Laplace transform, it returns an
% unevaluated call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% 
% But MATLAB "knows" that it is supposed to be a Laplace transform of a
% function. So if we compute the inverse Laplace transform, we obtain the
% original function

ilaplace(G)

%%
% The Laplace transform of a function is related to the Laplace transform 
% of its derivative:

syms g(t)
laplace(diff(g,t),t,s)


%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.
% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|
% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)| 
% 
% (in your answer, explain part (c) using comments).      
%
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.

syms t s 

%a) we will set the function f, and compute F using matlab
f = exp(2*t)*(t^3)
F = laplace(f)

%b) we will set the transform G, and compute g using matlab
G = ((s - 1)*(s - 2))/(s*(s + 2)*(s - 3))
g = ilaplace(G)

syms a t f(t)

%c) store the arbitrary function f(t) and use it in h(t) to find H
h(t) = exp(a*t)*f(t)
F = laplace(f(t))
H = laplace(h(t))

% based on displayed matlab results, 
% F is defined as laplace(f(t), t, s)
% H is defined as laplace(f(t), t, s-a)
% Therefore H is equal to F(s-a) because the variable in the s-domain is
% the only thing that has changed
% We know this is true based on our knowledge of laplace: multiplying a
% function by exp(at) in the time domain just means a shift of s in that
% domain by a.


%% Heaviside and Dirac functions
%
% These two functions are builtin to MATLAB: |heaviside| is the Heaviside
% function |u_0(t)| at |0|
%
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function (at |0|) is also defined with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)


%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% In your answer, explain the 'proof' using comments.

syms t s

a = 2; %gave |a| a value
u_a(t) = heaviside(t-a); %define u_a(t)
f(t) = exp(2*t)*(t^3); %took the same f(t) as question 1
g(t) = u_a(t)*f(t-a);

G = laplace(g(t))
F = laplace(f(t))

% by matlab results, we see that
% G = (6*exp(-2*s)) / (s - 2)^4
% F = 6/(s - 2)^4

% we can see that G can be rewritten as = F(s) * exp(-a*s)

% we know from theory that if we have a function g(t) that is defined using
% a heaviside with a shift of c, multiplied by a defined function f:
% u_c(t)*f(t-c), then the laplace should be exp(-c*t) * laplace(f(t)).
% This is the result we have found since here c = a = 2


%% Solving IVPs using Laplace transforms
%
% Consider the following IVP, |y''-3y = 5t| with the initial
% conditions |y(0)=1| and |y'(0)=2|.
% We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),1)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y


%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: Explain your steps using comments
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)|
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% * Is there an initial condition for which |y| remains bounded as |t| goes to infinity? If so, find it.

syms y(t) t Y s

ODE = diff(y(t),t,3)+2*diff(y(t),t,2)+diff(y(t),t,1)+2*y(t)+cos(t)==0;

L_ODE = laplace(ODE);

% Use the initial conditions

L_ODE = subs(L_ODE, y(0), 0);
L_ODE = subs(L_ODE, subs(diff(y(t), t), t, 0),0);
L_ODE = subs(L_ODE, subs(diff(y(t), t, 2), t, 0), 0);

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y);
Y=solve(L_ODE,Y);

y = ilaplace(Y)

ezplot(y,[0,10*pi])
title("Solution of the IVP");
xlabel("t")
ylabel("y")

% we can see that the behaviour of the ODE is to grow while oscillating, so
% there is no initial condition where |y| is bounded as t goes to infinity.
% we can see this from the graph because no matter what, it will behave
% like a growing oscillating function.

%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |y(0)=2 and y'(0)=1|
%
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%
% In your answer, explain your steps using comments.

syms y(t) t Y s;

% to get the piecewise functions at specific intervals we should use
% heaviside (it is t = 0, 2, 5)

u_0(t) = heaviside(t);
u_2(t) = heaviside(t-2);
u_5(t) = heaviside(t-5);

% graphed it on desmos to visualize the heaviside combinations.
% first part: 3*u_0 since it's just a line y=3 at t=0 to t=2
% second part: the ramp starts at t=2, so do (t-2)*u_2
% third part: the ramp stops at t=4 if g(t) = 5 so we want to start
% canceling out the ramp with a negative slope. do (-t+4)*u_5

g = @(t) 3*u_0(t)+ (t-2)*u_2(t) + (-t+4)*u_5(t);

%define the ode |y''+2y'+5y=g(t)|
ODE = diff(y(t),t,2) + 2*diff(y(t),t, 1) + 5*y(t) - g(t) == 0;
L_ODE = laplace(ODE);

%initial conditions
L_ODE = subs(L_ODE, y(0), 2);
L_ODE = subs(L_ODE, subs(diff(y(t), t), t, 0), 1);

%replace laplace with Y
L_ODE = subs(L_ODE, laplace(y(t), t, s), Y);
%solve for laplace
Y = solve(L_ODE, Y);

%solve for y using inverse laplace
y = ilaplace(Y)

%plot the solution
ezplot(y, [0, 12, 0, 2.25]);
title("Solution of the Piecewise ODE");
xlabel("t");
ylabel("y");

%verify the solution:
simplify(diff(y,t,2)+2*diff(y,t)+5*y-g)
%gives 3-3*heaviside(t) which is essentially 3 - 3 = 0


%% Exercise 5
%
% Objective: Use the Laplace transform to solve an integral equation
% 
% Verify that MATLAB knowns about the convolution theorem by explaining why the following transform is computed correctly.
syms t tau y(tau) s y(t)
I=int( exp(-2*(t-tau))* y(tau),tau,0,t)
laplace(I,t,s) % = laplace(y(t), t, s) / (s+2)

% I is defined as the integral from 0 to t of
% y(tau) * e^(-2*(t-tau)) * dtau
% which we can rewrite as the integral of f(t-tau) * g(tau) dtau, which is
% the definition of the convolution of (f*g)(t)

% based on our knowledge of convolution theorem
% this tells us that the laplace{f*g}(t) = laplace{f(t)} times laplace{g(t)}
% thus we expect laplace(I, t, s) to give us laplace{y(tau)} times
% laplace{exp(-2*(t-tau))}

% the integrals are being evaluated where tau goes from 0 to t, which means
% we can look at simply f(t) and g(t) by definition:
f = exp(-2*(t));
g = y(t);
%find the laplace
F = laplace(f)
G = laplace(g)
%multiply them
answer = F*G
answer - laplace(I,t,s) %answer is 0 so they are the same

% this gives us the same answer as taking the laplace of the integral I
% which confirms our understanding of how convolution relates to taking the
% laplace of functions within the integral definition.
% therefore, matlab understands convolution theorem, which states that
% L{(f*g)(t)} = L{f(t)} times L{g(t)}


