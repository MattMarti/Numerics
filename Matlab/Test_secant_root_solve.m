%% Test_secant_root_solve.m
% 
% Test case for the Secant Root Solver function.
% 
% Use Graphical technique, bisection method, false-position, fixed-point
% iteration, Netwon method, and secant method to find the first root of
%     f(x) = x*exp(x) - cos(x)
% 
% @author: Matt Marti
% @date: 2019-04-26

clear

% Define function
f = @(x) cos(x) - x.*exp(x);


%% Part C: Secant Method

% Parameters
a = 0; % Lower bound
b = 1; % Upper bound
errstop = 1e-12; % Stopping criteria
maxiter = 1000;

% Function call
x = secant_root_solve(f, a, b, maxiter, errstop);

% Check results
assert(abs(f(x)) < errstop, 'Results error not less than specified error');

