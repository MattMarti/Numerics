%% Test_cubicspline.m
% 
% Test case for the cubic spline function. Based on the driver script for
% the solution to AOE 4404 Assignment 5, problems 1 to 3. Cubic spline and
% Lagrange interpolation.
% 
% @author: Matt Marti
% @date: 2019-05-06

clear


%% Part 1
% Test the spline

% Given
f = @(x) sin(x);
df = @(x) cos(x);
xkvec = linspace(0, 10, 20);
fkvec = f(xkvec);
xinter = linspace(0, 10, 1000);
fslope = [ cos(xkvec(1)), cos(xkvec(end)) ]; % Clambed B.C.s

% Run function
cs = cubicspline(xkvec, fkvec, fslope);
[finter, dfinter] = cs.interp(xinter, true);

% 
% [finter, dfinter, akvec, bkvec, ckvec, dkvec, xstar] ...
%     = cubicspline(xkvec, fkvec, xinter, fslope);

% Test Function truth values
fitrue = f(xinter);
error = fitrue - finter;
maxerr = max(abs(error));
assert(maxerr < 2.5e-4, 'Spline error too high');

% Test Derivative truth values
dfitrue = df(xinter);
errord = dfitrue - dfinter;
maxerrd = max(abs(errord));
assert(maxerrd < 1.5e-3, 'Spline error too high');


%% Part 2
% Test spline function with multiple rows of data

% Given
f = @(x) [sin(x); cos(x)];
df = @(x) [cos(x); -sin(x)];
xkvec = linspace(0, 10, 20);
fkvec = f(xkvec);
xinter = linspace(0, 10, 1000);
fslope = [ cos(xkvec(1)), cos(xkvec(end)); -sin(xkvec(1)), -sin(xkvec(end)) ]; % Clambed B.C.s

% Run function
cs = cubicspline(xkvec, fkvec, fslope);
[finter, dfinter] = cs.interp(xinter, true);

% Test Function truth values
fitrue = f(xinter);
error = fitrue - finter;
maxerr = max(max(abs(error)));
assert(maxerr < 2.5e-4, 'Spline error too high');

% Test Derivative truth values
dfitrue = df(xinter);
errord = dfitrue - dfinter;
maxerrd = max(max(abs(errord)));
assert(maxerrd < 1.5e-3, 'Spline error too high');

