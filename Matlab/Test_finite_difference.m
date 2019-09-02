%% Test_finite_difference.m
% 
% Test case for the finite difference function.
% 
% @author: Matt Marti
% @date: 2019-05-06

clear


%% Test 1: Make sure all the orders work

% Function input
h = .1;
thist = (0:h:10)';
yhist = -9.81/2*thist.^2 + 50*thist + 200;

% Truth value
ydottruthhist = -9.81*thist + 50;

% Test Function call
errvec = [5e-1; 1e-11; 1e-11; 1e-11; 1e-11; 1e-11; 2e-11; 2.5e-11; 5e-11];
for n = 1:9
    [ydothist] = finite_difference(yhist, thist, n);
    errhist = ydothist - ydottruthhist;
    maxerr = max(abs(errhist));
    assert(maxerr < errvec(n), ...
        'Error for low order polynomial is not zero');
end


%% Test 2: Derivative of nonlinear function is approximate and all orders work
% Also test that multiple data entries work

% Function input
h = .1;
thist = (0:h:10)';
yhist = [sin(thist), cos(thist)];
n = 5;

% Truth value
ydottruthhist = [cos(thist), -sin(thist)];

% Function call
[ydothist] = finite_difference(yhist, thist, n);

% Test
% Test Function call
errvec = [5e-2; 5e-3; 2.5e-4; 2e-5; 2e-6; 5e-7; 2e-8; 2e-9; 1e-10; 1e-11; 2e-12];
for n = 1:11
    [ydothist] = finite_difference(yhist, thist, n);
    errhist = ydothist - ydottruthhist;
    maxerr = max(max(abs(errhist)));
    assert(maxerr < errvec(n), ...
        'Error for low order polynomial is not zero');
end

%% Test 3: Test data with uneven time vector values

% Function input
thist = zeros(1000, 1);
h = 0.1;
for ii = 2:length(thist)
    thist(ii) = thist(ii-1) + h*(mod(ii-1, 10)+1);
end
yhist = [5 + 6*thist + 4*thist.^2, 5 - 3*thist];
n = 5;

% Truth value
ydottruthhist = [6 + 8*thist, - 3*ones(size(thist))];

% Function call
[ydothist] = finite_difference(yhist, thist, n);

% Test
% Test Function call
errvec = [5e-8; 5e-8; 5e-8; 5e-8];
for n = 1:4
    [ydothist] = finite_difference(yhist, thist, n+1);
    errhist = ydothist - ydottruthhist;
    maxerr = max(max(abs(errhist)));
    assert(maxerr < errvec(n), ...
        'Error for low order polynomial is not zero');
end
