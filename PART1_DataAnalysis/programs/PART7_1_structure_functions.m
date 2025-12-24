clear
clc

addpath("../functions/");

params(1).name = 'veldata1.txt';
params(2).name = 'veldata2.txt';
params(3).name = 'veldata3.txt';

Nsamples = Inf;      % Size of the dataset for processing

%% ====> Parameters to be filled <====

scaling.S2_exp  = 2/3;   % scaling exponents
scaling.S3_exp  = 1;

params(1).lambda =  1.338e-02; 
params(1).L_C    = 0.2807;    

params(2).lambda = 1.003e-02 ;
params(2).L_C    = 0.3402;

params(3).lambda = 9.203e-03;
params(3).L_C    =  0.4002;

warning('Make sure to fill in the relevant values!');

%% --- Common setup ---
Num   = 100;
l_S2  = logspace(log10(0.001), log10(10), Num)';
l_S3  = logspace(log10(0.001), log10(1), Num)';
S2_all = zeros(3, Num);
S3_all = zeros(3, Num);

%% --- Loop through datasets and compute S2 and S3 ---
for i = 1:3
    % --- Read velocity data from file ---
    [u, f, U] = load_data(params(i).name, Nsamples);

    % --- Second Order Structure Function ---
    index = 2;
    S2 = structure_function(u, l_S2, f, U, index);
    SS2 = smooth(S2, 8);
    S2_all(i, :) = SS2;

    % --- Third Order Structure Function ---
    index = 3;
    S3 = structure_function(u, l_S3, f, U, index);
    SS3 = smooth(S3, 3);
    S3_all(i, :) = SS3;
end

%% --- Plot S2 ---
figure;
for i = 1:3
    loglog(l_S2, S2_all(i,:), 'linewidth', 2);
    hold on;
end
loglog(l_S2, 15 * l_S2.^scaling.S2_exp, 'k--', 'linewidth', 1);
xlim([min(l_S2) max(l_S2)]);
ylim([5e-3 100]);
xlabel('$l\;[m]$', 'Interpreter', 'latex', 'fontsize', 14);
ylabel('$S_2(l)\;[m^2/s^2]$', 'Interpreter', 'latex', 'fontsize', 14);
legend({'data1', 'data2', 'data3', 'scaling'}, 'Interpreter', 'latex', 'fontsize', 12);

for i = 1:3
    xline(params(i).lambda, '--', ['$\lambda_{' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    xline(params(i).L_C, '--', ['$L\_C_{' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end


%% --- Plot S3 ---
figure;
for i = 1:3
    loglog(l_S3, S3_all(i,:), 'linewidth', 2);
    hold on;
end
loglog(l_S3, 50 * l_S3.^scaling.S3_exp, 'k--', 'linewidth', 1);
xlim([min(l_S3) max(l_S3)]);
ylim([1e-4 1e2]);
xlabel('$l\;[m]$', 'Interpreter', 'latex', 'fontsize', 14);
ylabel('$S_3(l)\;[m^3/s^3]$', 'Interpreter', 'latex', 'fontsize', 14);
legend({'data1', 'data2', 'data3', 'scaling'}, 'Interpreter', 'latex', 'fontsize', 12);

for i = 1:3
    xline(params(i).lambda, '--', ['$\lambda_{' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    xline(params(i).L_C, '--', ['$L\_C_{' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end
hold on;
for i = 1:3
    % Taylor microscale λ
    xline(params(i).lambda, '--', ...
        ['$\lambda_{' num2str(i) '}$'], ...
        'Interpreter','latex',...
        'FontSize',10,...
        'LabelVerticalAlignment','bottom',...
        'HandleVisibility','off',...
        'LineWidth',1.5);
    % Integral scale L_C
    xline(params(i).L_C, '--', ...
        ['$L_{C,' num2str(i) '}$'], ...
        'Interpreter','latex',...
        'FontSize',10,...
        'LabelVerticalAlignment','bottom',...
        'HandleVisibility','off',...
        'LineWidth',1.5);
end

% Linear fit 

zeta3_data = zeros(1,3);
for i = 1:3
    % define inertial‐range mask (tweak the 1.1/0.9 factors if you like)
    mask = (l_S3 > 1.1*params(i).lambda) & (l_S3 < 0.9*params(i).L_C);
    % fit log|S3| vs. log(l) and extract slope
    p = polyfit( log(l_S3(mask)), log(abs(S3_all(i,mask))) , 1 );
    zeta3_data(i) = p(1);
    fprintf('Dataset %d: measured zeta_3 = %.3f\n', i, zeta3_data(i));
end
% compare to K41 prediction
fprintf('K41 predicts zeta_3 = %.3f\n', scaling.S3_exp);