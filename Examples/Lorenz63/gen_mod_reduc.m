function [x,a,b] = gen_mod_reduc(Z,PSI)

% Model reduction parameters
p = 4;      % Order of AR term
q = 4;      % Order of X term
r = p/2;    

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
A = get_A(r);                   % Stability condition
C = ones(3*r,1);

ERROR = @(Alpha) Error(Alpha,Z,d,r,p,q,steps,psi,nu)...
    + 10^8*heaviside(A*[Alpha(1:p)]'-C);

Alpha_init = zeros(1, p + nu*(q+1));

options = optimoptions('lsqnonlin',...
                'Display','iter');
                %'UseParallel','true');
alpha_sol = lsqnonlin(ERROR,Alpha_init,[],[],options);

alpha = alpha_sol(1:r);
beta = alpha_sol(r+1:p);
b = alpha_sol(p+1:end);

a = get_a(alpha, beta);
b = reshape(b, nu, q+1);

% Run the reduced model
x = zeros(d,steps);
y = zeros(d,steps);
% We take as the initial conditions the first obj.p observations 
x(:,1:p+1) = Z(:,1:p+1);
y(:,1:p) = Z(:,1:p);

for n = p + 2 : steps
    y_temp = zeros(d,1);
    for i = 1:q+1                   
        y_temp = y_temp + PSI(x(:,n-p+i-2))*b(:,i); 
    end
    % Here we add the autoregressive part. This step could also
    % be accomplished by the 'cascade.' Also, a(1) =
    % a_p-1,...,a(p) = a_0, so I flip it soa s to get
    % a_p-1 y_n-1 + ... + a_0 y_n-p as required
    y(:,n - 1) = y_temp - y(:,n-p-1:n-2)*a;
    x(:,n) = y(:,n - 1); + .001*mvnrnd(zeros(1,d),eye(d))';
end
end


% Form Energy function. This will be our objective function.
function error = Error(Alpha,Z,d,r,p,q,steps,psi,nu)
    alpha = Alpha(1:r);
    beta = Alpha(r+1:2*r);
    b = Alpha(2*r+1:end);
    b = reshape(b,nu,q+1);
    a = get_a(alpha,beta);
    Y = zeros(d,steps);
    Y(:,1:p) = Z(:,1:p);
    
    for n = p + 2 : steps
        Y_temp = zeros(d,1);
        for i = 1:q+1
            Y_temp = Y_temp + psi(:,:,n-p+i-2)*b(:,i);  
        end
        Y(:,n) = Y_temp - Y(:,n-p-1:n-2)*a;
    end
    ZZ = reshape(Z,1,d*steps);
    YY = reshape(Y,1,d*steps);
    error = ZZ - YY;
end

function a = get_a(alpha,beta)
    r = size(alpha,2);
    a = zeros(2*r,1);
    a(4) = alpha(1) + alpha(2);
    a(3) = beta(1) + beta(2) + alpha(1)*alpha(2);
    a(2) = alpha(1)*beta(2) + alpha(2)*beta(1);
    a(1) = beta(1)*beta(2);
end

function A = get_A(r)
d1 = [0; 1; -1];
d2 = [1; -1; -1];
A = [kron(eye(r),d1),kron(eye(r),d2)];
end