clc;
clear;

addpath("../functions/");

%% --- Common Setup ---
params(1).name = 'veldata1.txt';
params(2).name = 'veldata2.txt';
params(3).name = 'veldata3.txt';

%% --- Loading data ---
[u1, sample_freq1, U1] = ...
    load_data(params(1).name, Inf);
[u2, sample_freq2, U2] = ...
    load_data(params(2).name, Inf);
[u3, sample_freq3, U3] = ...
    load_data(params(3).name, Inf);

params(1).u = u1;
params(2).u = u2;
params(3).u = u3;

params(1).U = U1;
params(2).U = U2;
params(3).U = U3;

params(1).f = sample_freq1;
params(2).f = sample_freq2;
params(3).f = sample_freq3;
% You need params(i).L_C from part 2

params(1).dl = params(1).U / params(1).f;
params(2).dl = params(2).U / params(2).f;
params(3).dl = params(3).U / params(3).f;

lmax_autocorr = 10;               % Maximum autocorr. length
C1 = autocorrelation(u1, params(1).dl, lmax_autocorr);
C2 = autocorrelation(u2, params(2).dl, lmax_autocorr);
C3 = autocorrelation(u3, params(3).dl, lmax_autocorr);
C1 = C1 / C1(1);
C2 = C2 / C2(1);
C3 = C3 / C3(1);


% Length vectors
dl1 = params(1).U / params(1).f;
dl2 = params(2).U / params(2).f;
dl3 = params(3).U / params(3).f;


l1 = (0:numel(C1)-1)' * dl1;
l2 = 0:dl2:(params(2).dl*(length(C2) - 1));
l3 = 0:dl3:(params(3).dl*(length(C3) - 1));

warning('Make sure to fill Lc properly')

%Find intersections
params(1).L_C = 1;
params(2).L_C = 1;
params(3).L_C = 1;

for k = 1:3
  Cl = eval(['C' num2str(k)]);
  ll = eval(['l' num2str(k)]);
  idx = find(Cl <= exp(-1), 1);
  if isempty(idx) || idx==1
    Lc = ll(end);
  else
    % interpolate between (idx-1) and idx
    x0 = ll(idx-1);  x1 = ll(idx);
    y0 = Cl(idx-1);  y1 = Cl(idx);
    Lc = x0 + (exp(-1)-y0)*(x1-x0)/(y1-y0);
  end
  params(k).L_C = Lc;
end

%% --- Dissipation and Reynolds numbers ---
nu = 1.5e-5;   % kinematic viscosity [m^2/s], e.g. air at 20 °C
outer_scale = 1;
for i = 1:1:3
    u_i       = params(i).u;               % velocity fluctuation series
    var_u     = mean(u_i.^2);              % ⟨u^2⟩
    LC        = params(i).L_C;             % correlation length [m]
    params(i).epsilon =  0.5 * (var_u)^(3/2) / LC;
    fprintf('Dataset %d: ε = %.3e  ( using ⟨u^2⟩=%.3e, L_C=%.3e )\n', ...
            i, params(i).epsilon, var_u, LC);
   
    lambda = sqrt(15 * nu * var_u / params(i).epsilon);
    % Taylor‐Reynolds number
    Re_lambda = sqrt(var_u) * lambda / nu;
    % store in params and print
    params(i).lambda     = lambda;
    params(i).Re_lambda  = Re_lambda;
    fprintf('Dataset %d: λ = %.3e m,  Re_λ = %.1f\n', ...
            i, lambda, Re_lambda);
    urms  = sqrt(mean(u_i.^2));  
    params(i).Re_outer = urms * LC / nu;   % Re = u′ L_C / ν
    fprintf('Dataset %d: Re_outer = %.1f\n', i, params(i).Re_outer);
end

warning('Make sure that you filled all of the above numbers')