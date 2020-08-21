%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file: Langevin equation solver
% author: Jared McBride (Mar 19, 2019)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Euler
tic
mu = 0;
sig = 1;

xo = .1;
start = 0;
stop = 10000;
steps = 10e5+1;

V = @(x) .25*(x^2 - 1)^2;
Vp = @(x) -x*(x^2 - 1);

t = linspace(start,stop,steps);
h = t(2) - t(1);

x = zeros(1,steps);
x(1) = xo;
for n = 1 : steps - 1
    x(n+1) = x(n) + h*Vp(x(n)) + sqrt(h)*randn(1);
end

%
disc = round(stop/2);
z = x(disc:end);
t_z = t(disc:end);
steps = size(z,2);

PSI = Vp;

a = 0;
b = [0 0];

[y1,a1,b1] = gen_mod_reduc_L(z,PSI,a,b);
[y2,a2,b2] = local_mod_reduc_L(z,PSI,a,b);

[a1,b1]
[a2,b2]

subplot(3,1,1), plot(t,x);
subplot(3,1,2), plot(t_z,y1);
subplot(3,1,3), plot(t_z,y2);
toc

function [y,a,b] = local_mod_reduc_L(x,PSI,a,b)
    N = size(x,2);
    Ao = [a,b];
    obj_fun = @(A) Pre_obj_fun(x,PSI,A(1),[A(2),A(3)],N);
    options = optimoptions('lsqnonlin',...
                'Display','iter');
                %'UseParallel','true');
	A_sol = lsqnonlin(obj_fun,Ao,[],[],options);
    a = A_sol(1);
    b = A_sol(2:3);
    y = zeros(1,N);
    Y = zeros(1,N-1);
    
    y(1:2) = x(1:2);
    Y(1) = x(2);
    for n = 2 : N - 1
        Y(n) = -a*Y(n-1) + PSI(y(n-1))*b(1) + PSI(y(n))*b(2);
        y(n+1) = Y(n);
    end
end

function E = Pre_obj_fun(x,PSI,a,b,N)    
    y = zeros(1,N);
    Y = zeros(1,N-1);
    E = zeros(N-1,1);

    y(1:2) = x(1:2);
    Y(1) = x(2);
    for n = 2:N-1
        Y(n) = -a*Y(n-1) + PSI(x(n-1))*b(1) + PSI(x(n))*b(2);
        y(n+1) = Y(n);
        E(n+1) = y(n+1) - x(n+1);
    end
end