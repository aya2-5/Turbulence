close all; clear; clc

%% inputs
addpath('../functions/')       % folder containing functions
L = 38.6;                      % domain length
N = 64;                        % spatial resolution
symm = true;                   % imposed center symmetry
t_study = 2500.0;              % analysis time period
dt = 0.1;                      % time step size for time integration
dt_store = 1.0;                % time intervals of storing a snapshot
T_max = 120.0;                 % maximum T for recurrent flow analysis
T_eqb = 10.0;                  % time interval for computing equilibria
epsilon = 1e-6;                % perturbation amplitude for computing J

%% initial condition
[x,~] = domain(L,N);           % construct the spatial domain
u0 = sin(2.0*pi*x/L);          % initial condition in physical state
v0 = field2vector(u0,N,symm);  % initial state vector

%% time integrate the KSE
[V,t_vec] = KSE_integrate(v0,t_study,dt,dt_store,L,N,symm);

%% computing equilibria
snapshots = [   1,  501, 1001, 1501, 2001, 2501 ];
numSnaps = numel(snapshots);

% Preallocate storage for guesses, solutions, flags, and times
t_snap = zeros(numSnaps,1);
v_guesses = cell(numSnaps,1);
v_sols    = cell(numSnaps,1);
flags      = false(numSnaps,1);

for j = 1:numSnaps
    idx = snapshots(j);
    t_snap(j) = t_vec(idx);              % record time of snapshot
    v_guess   = V(:,idx);                % initial guess from trajectory
    [v_sol, flag] = search4EQ(v_guess, T_eqb, dt, L, N, symm);

    v_guesses{j} = v_guess;
    v_sols{j}    = v_sol;
    flags(j) = (flag==1);
end

% Convert cells to matrices for projection
M_guess = cell2mat(v_guesses');
M_sol   = cell2mat(v_sols');

% Compute E, P, D for guesses and converged equilibria
[E_guess, P_guess, D_guess] = projection(M_guess, N, L, symm);
[E_sol,   P_sol,   D_sol]   = projection(M_sol,   N, L, symm);

% Plot C: initial guesses vs. converged equilibria in physical space
figure('Name','Plot C: Equilibrium Search Results');
for j = 1:numSnaps
    subplot(2,3,j);
    u_guess = vector2field(v_guesses{j}, N, symm);
    u_sol   = vector2field(v_sols{j},    N, symm);
    plot(x, u_guess, '--', 'DisplayName','Initial Guess'); hold on;
    plot(x, u_sol,   '-',  'DisplayName','Converged EQ');
    title(sprintf('t = %.1f, flag = %d', t_snap(j), flags(j)));
    if j == 1, legend('Location','best'); end
    xlabel('x'); ylabel('u(x)');
end

% Display results in a table: snapshot index, time, success, E and P=D of solution
results = table(snapshots', t_snap, flags, E_sol, P_sol, D_sol, ...
    'VariableNames', {'Index','Time','Success','E_eq','P_eq','D_eq'});
disp('Equilibrium Search Results:');
disp(results);


%% equilibria from sinusoidal guesses (Plot D)
Kmax       = 6;
v_sin      = cell(Kmax,1);
v_sin_sol  = cell(Kmax,1);
flags_sin  = false(Kmax,1);

for k = 1:Kmax
    % build the k-th sinusoidal guess in physical space
    u_guess    = sin(k * 2*pi*x/L);
    % convert to state vector
    v_guess    = field2vector(u_guess, N, symm);
    
    % search for equilibrium
    [v_sol, flag] = search4EQ(v_guess, T_eqb, dt, L, N, symm);
    
    % store
    v_sin{k}      = v_guess;
    v_sin_sol{k}  = v_sol;
    flags_sin(k)  = (flag==1);
end

% project the converged solutions
M_sin_sol = cell2mat(v_sin_sol');
[E_sin, P_sin, D_sin] = projection(M_sin_sol, N, L, symm);

% plot Plot D: each guess vs. its equilibrium
figure('Name','Plot D: Sinusoidal Guesses');
for k = 1:Kmax
    subplot(2,3,k);
    plot(x, vector2field(v_sin{k},    N, symm), '--'); hold on;
    plot(x, vector2field(v_sin_sol{k},N, symm), '-');
    title(sprintf('k = %d, flag = %d', k, flags_sin(k)));
    if k==1, legend('guess','equilibrium','Location','best'); end
    xlabel('x'); ylabel('u(x)');
end

% display results in a table
sin_results = table((1:Kmax)', flags_sin, E_sin, P_sin, D_sin, ...
    'VariableNames', {'k','Success','E_eq','P_eq','D_eq'});
disp('Sinusoidal‐guess Equilibrium Search Results:');
disp(sin_results);



% filter only the equilibria that converged to something nonzero
valid_snapshots = flags  & (E_sol>0);
valid_sinusoids = flags_sin & (E_sin>0);

E_eqs    = E_sol(valid_snapshots);
P_eqs    = P_sol(valid_snapshots);
D_eqs    = D_sol(valid_snapshots);

E_sineq  = E_sin(valid_sinusoids);
P_sineq  = P_sin(valid_sinusoids);
D_sineq  = D_sin(valid_sinusoids);


%% projection of the full trajectory and overlay equilibria (Plot E)
[E_traj, P_traj, D_traj] = projection(V, N, L, symm);

figure('Name','Plot E: Attractor in (E,P,D)');
plot3(E_traj, P_traj, D_traj, '.', 'MarkerSize',2);  hold on;

% now overlay only the valid equilibria, in the correct (x=E,y=P,z=D) order:
scatter3(E_eqs,   P_eqs,   D_eqs,   80, 'ro', 'filled');
scatter3(E_sineq, P_sineq, D_sineq, 80, 'ks', 'filled');

xlabel('E');  ylabel('P');  zlabel('D');
legend('Trajectory','Eq (snapshots)','Eq (sinusoids)','Location','best');
view(3);



%% computing UPOs

% Plot F: recurrence indicator r(t,T) = ||v(t+T)-v(t)||/||v(t)||

% assume V is size [nDim × nSnapshots], t_vec is [1 × nSnapshots]
nSnap   = size(V,2);
dtS     = dt_store;                  % your storage interval, e.g. 1.0
maxK    = floor(T_max / dtS);        % number of T‐steps up to T_max
Tvals   = (1:maxK)*dtS;              % row vector of T

% preallocate
r_mat   = nan(nSnap, maxK);

% precompute norms of each snapshot
normV   = sqrt( sum( V.^2, 1 ) );    % 1×nSnap

for k = 1:maxK
    % difference between v(t+T_k) and v(t)
    D = V(:,1+k:end) - V(:,1:end-k);     % size nDim × (nSnap−k)
    ds = sqrt( sum(D.^2,1) );            % 1×(nSnap−k)
    r_mat(1:end-k, k) = ds ./ normV(1:end-k);
end

% take log, avoid log(0)
lnr = log( r_mat + eps );

% now plot as a heatmap
figure('Name','Plot F: ln r(t,T) recurrence map');
imagesc( t_vec, Tvals, lnr.' );    % transpose so rows→T, cols→t
axis xy;                           % put origin at lower left
xlabel('t'); ylabel('T');
title('ln r(t,T) = ln (||v(t+T)-v(t)||/||v(t)||)');
colorbar;

% (optional) if you prefer a mesh/ surface:
% [Tg, tg] = meshgrid(Tvals, t_vec);
% figure;
% surf(tg, Tg, lnr,'EdgeColor','none');
% view(2); colorbar; xlabel('t'); ylabel('T'); title('ln r(t,T)');

%% computing UPOs (Plot G)

% 1) Your five hand-picked minima from the recurrence plot:
t_guesses = [450, 1050, 1230, 1720, 1820];   % example times
T_guesses = [ 59,  39,   78,   31,   58];   % example periods


nG = numel(t_guesses);
po_sols    = cell(nG,1);
flags_po   = false(nG,1);
T_converged= nan(nG,1);

% 2) Loop over each guess, call search4PO
for i = 1:nG
    % find the nearest stored snapshot to t_guesses(i)
    [~, idx] = min(abs(t_vec - t_guesses(i)));
    v0       = V(:,idx);               % initial state guess
    T0       = T_guesses(i);           % period guess

    [v_best, T_best, flag] = search4PO(v0, T0, dt, L, N, symm);

    po_sols{i}     = v_best;
    flags_po(i)    = (flag == 1);
    T_converged(i) = T_best;
end

% 3) Project converged orbits into (E,P,D)
M_po = cell2mat(po_sols(flags_po));    % only the successful ones
[E_po, P_po, D_po] = projection(M_po, N, L, symm);

% 4) Build the results table
Guess    = (1:nG).';
t0       = t_guesses.';
T0       = T_guesses.';
Success  = flags_po;
T_conv   = T_converged;

% For unsuccessful guesses pad with NaNs
E = nan(nG,1); P = nan(nG,1); D = nan(nG,1);
E(flags_po) = E_po;
P(flags_po) = P_po;
D(flags_po) = D_po;

upo_results = table(Guess, t0, T0, Success, T_conv, E, P, D, ...
    'VariableNames',{'#','t_{guess}','T_{guess}','Converged','T_{sol}','E','P','D'});

disp('Unstable Periodic Orbit Search Results:');
disp(upo_results);

%% Plot G: Chaotic attractor + computed UPOs

% 1) Re‐plot the full attractor in (E,P,D)
figure('Name','Plot G: Attractor & UPOs'); clf;
[E_traj, P_traj, D_traj] = projection(V, N, L, symm);
plot3(E_traj, P_traj, D_traj, '.', 'MarkerSize', 2, 'Color', [0.6,0.6,0.6]);
hold on;

% 2) Loop over each successful UPO
colors = lines(sum(flags_po));  % distinct colors for each
ipo = 0;
for i = 1:numel(flags_po)
    if ~flags_po(i), continue; end
    ipo = ipo + 1;

    % time‐march from the converged periodic point for one period
    v_init = po_sols{i};
    Tpo    = T_converged(i);
    [Vpo, ~] = KSE_integrate(v_init, Tpo, dt, dt_store, L, N, symm);

    % project to (E,P,D)
    [Epo, Ppo, Dpo] = projection(Vpo, N, L, symm);

    % plot the UPO orbit
    plot3(Epo, Ppo, Dpo, '-', 'LineWidth', 2, 'Color', colors(ipo,:));
end

% 3) Finish styling
xlabel('E'); ylabel('P'); zlabel('D');
title('Plot G: Chaotic attractor (gray) with UPOs overlaid');
legend_entries = [{'Attractor'}; arrayfun(@(k) sprintf('UPO %d',k), 1:ipo,'uni',0)'];
legend(legend_entries,'Location','best');
view(3);
grid on;


%% ------------------------------------------------------------------------
% Compute leading Floquet multiplier for each converged UPO
nDim    = size(V,1);
nG      = numel(po_sols);
lambda1 = nan(nG,1);

for i = 1:nG
    if flags_po(i)
        % converged periodic point and period
        vp   = po_sols{i};
        Tpo  = T_converged(i);
        
        % 1) compute unperturbed end‐state once
        Vref = KSE_integrate(vp, Tpo, dt, Tpo, L, N, symm);
        vend0 = Vref(:,end);
        
        % 2) power‐method to find leading eigenvalue of DF^T(vp)
        w = randn(nDim,1);
        w = w / norm(w);
        for it = 1:20
            % finite‐difference of the map
            vp_pert = vp + epsilon * w;
            Vp = KSE_integrate(vp_pert, Tpo, dt, Tpo, L, N, symm);
            vend = Vp(:,end);
            
            % action of Jacobian on w
            u = (vend - vend0) / epsilon;
            
            % Rayleigh quotient / norm gives λ
            lambda = norm(u);
            
            % normalize for next iterate
            w = u / lambda;
        end
        
        lambda1(i) = lambda;
    end
end

% 3) Add to results table and display
upo_results.lambda = lambda1;

disp('Unstable Periodic Orbit Search Results (with leading Floquet multiplier):');
disp(upo_results);


%% ------------------------------------------------------------------------
% Compute divergence time‐scale t*_F for each UPO
tFstar = nan(nG,1);
for i = 1:nG
    if flags_po(i) && lambda1(i)~=0
        tFstar(i) = T_converged(i) / log(abs(lambda1(i)));
    end
end

% Add to results table and display
upo_results.tFstar = tFstar;
disp('Unstable Periodic Orbit Search Results (with λ₁ and t*ₙ):');
disp(upo_results);
