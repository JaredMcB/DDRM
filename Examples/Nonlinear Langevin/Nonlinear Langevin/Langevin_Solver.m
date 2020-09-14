%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file: Langevin equation solver
% author: Jared McBride (Mar 19, 2019)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forward Euler
m = 1;
lam = 2;

mu = 0;
sig = 10;

xo = [0;1];
start = 0;
stop = 10;
steps = 10e4;

F = @(x) [0 1; 0 -lam/m]*x;

plot(t,x(1,:));

% [x] = ForwardEuler(xo,F,start,stop,steps);

function [x] = RK23(xo,F,start,stop,steps)
t = linspace(start,stop,steps);
h = t(2) - t(1);

x = zeros(2,steps);
x(:,1) = xo;
for n = 1 : steps - 1
    nu = sig*randn(1) + mu;
    
    Y1 = x(:,n);
    Y2 = x(:,n) + h/2*;
    Y3 = 0;
    Y4 = 0;
end
end




function [x] = ForwardEuler(xo,F,start,stop,steps)
t = linspace(start,stop,steps);
h = t(2) - t(1);

x = zeros(2,steps);
x(:,1) = xo;
for n = 1 : steps - 1
    nu = sig*randn(1) + mu;
    x(:,n+1) = x(:,n) +h*(F(x(:,n)) + [0; nu]);
end

plot(t,x(1,:));
end