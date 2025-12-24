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


%% --- Frozen flow hypothesis ---
params(1).U = U1;
params(2).U = U2;
params(3).U = U3;
params(1).f = sample_freq1;
params(2).f = sample_freq2;
params(3).f = sample_freq3;

params(1).x = params(1).U * (0:(numel(u1)-1))' / params(1).f;
params(2).x = params(2).U * (0:(numel(u2)-1))' / params(2).f;
params(3).x = params(3).U * (0:(numel(u3)-1))' / params(3).f;

warning('Make sure to properly fill the above parameters')


%% --- 1.2 Report mean velocities from load_data ---
fprintf('\nMean velocities U (m/s):\n');
fprintf('  Dataset 1: %.4f\n', U1);
fprintf('  Dataset 2: %.4f\n', U2);
fprintf('  Dataset 3: %.4f\n\n', U3);

%% --- Plot A ---
figure()
hold on; grid on;
plot(params(1).x, u1 + params(1).U)
plot(params(2).x, u2 + params(2).U)
plot(params(3).x, u3 + params(3).U)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$U(x)$', 'Interpreter', 'latex', 'FontSize', 14)
legend({'$U_1$', '$U_2$', '$U_3$'}, 'Interpreter', 'latex', 'FontSize', 14)

%% --- Turbulence intensity ---
params(1).I = sqrt(mean(u1.^2)) / params(1).U;
params(2).I = sqrt(mean(u2.^2)) / params(2).U;
params(3).I = sqrt(mean(u3.^2)) / params(3).U;

warning('Make sure to fill the turbulence intensity values !');

fprintf('\nTurbulence intensities I (–):\n');
fprintf('  Dataset 1: %.4f\n', params(1).I);
fprintf('  Dataset 2: %.4f\n', params(2).I);
fprintf('  Dataset 3: %.4f\n\n', params(3).I);
