close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
dt = 0.1;                      % time step size for time integration
epsilon = 1e-6;                % perturbation amplitude for computing J
tau = 2;                       % time intervals of computing J
N_norm = 5000;                 % number of re-normalizations
N_exp = 10;                    % number of Lyapunov exponents

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% compute Lyapunov exponents
Q = zeros(length(v0),N_exp);
for i = 1:N_exp
  Q(i,i) = 1;
end   % allocate matrix of orthonormal vectors

X = zeros(N_exp,N_norm);       % allocate matrix of Lyapunov exp. history
t = zeros(N_norm,1);           % allocate vector of normalization instances
S=zeros(N_exp,1);

figure
for i = 1:N_norm
    J = Jacobian(v0,tau,epsilon,dt,L,N,symm);
    
    %%% to be completed
        
    % Map the p deviation vectors through the Jacobian
    V = J*Q;                    % n×p

    % Orthonormalize via economy QR
    [Q,R] = qr(V,0);            % Q is n×p, R is p×p upper‐triangular

    % accumulate the local stretch factors
    r = abs(diag(R));           % p×1
    S = S + log(r);             % sum_j=1^i log‖R_jj‖

    % current Lyapunov estimates χ^(i) = S/(i·τ)
    X(:,i) = S/(i*tau);

    
    [v0,~] = KSE_integrate(v0,tau,dt,0,L,N,symm);
    t(i) = i*tau;
    
    if(rem(i,20)==0)           % update figure every 20 re-normalizations
        clf; grid on; hold on
        for q = 1:N_exp
            plot(t(1:i),X(q,1:i),'LineWidth',2)
        end
        xlabel('t'); ylabel('\chi_i')
        labels = arrayfun(@(j) sprintf('\\chi_{%d}',j), 1:N_exp, 'uni',false);
        legend(labels,'Location','best')
        drawnow
    end
end

chi = S/(N_norm*tau);          % p×1 vector of exponents
disp('Leading 10 Lyapunov exponents:')
disp(chi)

T_L = 1/chi(1);                % Lyapunov time
fprintf('Lyapunov time t_L* = %g\n',T_L)
