close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
T_trans = 1000.0;              % transient time period
T_study = 250.0;               % analysis time period
dt = 0.1;                      % time step size for time integration
dt_store = 1.0;                % time intervals of storing a snapshot
epsilon = 1e-2;                % relative ampliture of perturbation

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% transient time integration
[v1000,~] = KSE_integrate(v0,T_trans,dt,0,L,N,symm);

%% Step 2: Reference trajectory u1(x,t)
[v1_snapshots, t] = KSE_integrate(v1000, T_study, dt, dt_store, L, N, symm);
num_snaps = size(v1_snapshots, 2);
u1 = zeros(N, num_snaps);
for j = 1:num_snaps
    u1(:, j) = vector2field(v1_snapshots(:, j), N, symm);
end


%% Step 3: Perturb the final state and run perturbed trajectory u2(x,t)
r = zeros(size(v1000));         % preallocate perturbation vector
for k = 1:length(r)
    eta = rand();               
    r(k) = epsilon * (2*eta - 1) * v1000(k);
end
v0p = v1000 + r;                 % perturbed initial state
[v2_snapshots, ~] = KSE_integrate(v0p, T_study, dt, dt_store, L, N, symm);
% Convert to physical field
u2 = zeros(N, num_snaps);
for j = 1:num_snaps
    u2(:, j) = vector2field(v2_snapshots(:, j), N, symm);
end

%% Step 4: Plot results (Plot A)
u_diff = abs(u1 - u2);

figure;
subplot(3,1,1);
contourf(t, x, u1, 20, 'LineColor', 'none');
title('u_1(x,t) - reference');
xlabel('t'); ylabel('x'); colorbar;

subplot(3,1,2);
contourf(t, x, u2, 20, 'LineColor', 'none');
title('u_2(x,t) - perturbed');
xlabel('t'); ylabel('x'); colorbar;

subplot(3,1,3);
contourf(t, x, u_diff, 20, 'LineColor', 'none');
title('|u_1(x,t) - u_2(x,t)|');
xlabel('t'); ylabel('x'); colorbar;

