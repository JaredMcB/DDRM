function [x,a,b] = gen_mod_reduc_L(Z,PSI,a,b)

% Model reduction parameters
p = 1;      % Order of AR term
q = 1;      % Order of X term
r = floor(p/2);    

% Information from inputs
[d, nu] = size(PSI(Z(:,1)));   % Number of observered Variables
                               % and Number of basis functions
steps = size(Z,2);             % Number of steps


% Compute psi sequence
psi = zeros(d, nu, steps);
for i = 1:steps
    psi(:,:,i) = PSI(Z(:,i));
end

% Fit the reduced model
ERROR = @(Alpha) Error(Alpha,Z,d,r,p,q,steps,psi,nu)...
    + 10^8*stab_con(Alpha(1:p));

Alpha_init = zeros(1, p + nu*(q+1));
Alpha_init = [a,b];

options = optimoptions('lsqnonlin',...
                'Display','iter');
                %'UseParallel','true');
alpha_sol = lsqnonlin(ERROR,Alpha_init,[-1 -1 -1],[1 1 1],options);

e = ERROR(alpha_sol);
norm(e,2)

alpha = alpha_sol(1:p-r);
beta = alpha_sol(p-r+1:p);
b = alpha_sol(p+1:end);

a = get_a(alpha, beta);
b = reshape(b, nu, q+1);



% e = ERROR([a,b]);
% norm(e,2)
% disp(sprintf('Residual: %g', norm(e,2)));

% Run the reduced model
x = zeros(d,steps+1);
y = zeros(d,steps);
% We take as the initial conditions the first obj.p observations 
x(:,2:p+1) = Z(:,1:p);
y(:,1:p) = Z(:,1:p);

for n = p + 1 : steps
    psi_sum = zeros(d,1);
    for i = 1:q+1                   
        psi_sum = psi_sum + PSI(x(:,n-p+i-1))*b(:,i); 
    end
    % Here we add the autoregressive part. This step could also
    % be accomplished by the 'cascade.' Also, a(1) =
    % a_p-1,...,a(p) = a_0, so I flip it soa s to get
    % a_p-1 y_n-1 + ... + a_0 y_n-p as required

    y(:,n) = psi_sum - y(:,n-p:n-1)*a;
    x(:,n+1) = y(:,n) + .0*mvnrnd(zeros(1,d),eye(d))';
end

x = x(:,2:end);
end


% Form Energy function. This will be our objective function.
function error = Error(Alpha,Z,d,r,p,q,steps,psi,nu)
    alpha = Alpha(1:p-r);
    beta = Alpha(p-r+1:p);
    b = Alpha(p+1:end);
    b = reshape(b,nu,q+1);
    a = get_a(alpha,beta);
    Y = zeros(d,steps);
    X = zeros(d,steps+1);
    Y(:,1:p) = Z(:,1:p);
    
    for n = p + 1 : steps
        psi_sum = zeros(d,1);
        for i = 1:q+1
            psi_sum = psi_sum + psi(:,:,n-p+i-1)*b(:,i);  
        end
        Y(:,n) = psi_sum - Y(:,n-p:n-1)*a;
        X(:,n+1) = Y(:,n);
    end
    ZZ = reshape(Z,1,d*steps);
    XX = reshape(X(:,2:end),1,d*steps);
    error = ZZ - XX;
end

function C = stab_con(a)
p = size(a,2);
r = floor(p/2);
d1 = [0; 1; -1];
d2 = [1; -1; -1];
A1 = [kron(eye(r),d1)];
A2 = [kron(eye(r),d2)];
A0 = [zeros(3*r,1)];
B = [zeros(2,r),[1;1],zeros(2,r)];
if mod(p,2) == 1
    A = [A1,A0,A2;B];
    c = ones(3*r + 2,1);
else
    A = [A1,A2];
    c = ones(3*r,1);
end
C = max(heaviside(A*a'-c));
if r == 0 
    C = 0; 
end
end

function a = get_a(alpha,beta)
    p = size(alpha,2)+size(beta,2);
    r = floor(p/2);
    syms x
    if p ==1
        a = alpha;
    else
        if mod(p,2) == 0 y = 1; end
        if mod(p,2) == 1 y = x + alpha(r+1); end
        for i = 1 : r
            y = y*poly2sym([1, alpha(i), beta(i)]);
        end
        y = expand(y);
        a = sym2poly(y);
        a = a(2:end)';
    end
end
