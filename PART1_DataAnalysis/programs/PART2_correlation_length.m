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

% If not set
params(1).U = U1;
params(2).U = U2;
params(3).U = U3;
params(1).f = sample_freq1;
params(2).f = sample_freq2;
params(3).f = sample_freq3;

%% --- Correlation length ---
params(1).dl = params(1).U / params(1).f;
params(2).dl = params(2).U / params(2).f;
params(3).dl = params(3).U / params(3).f;

warning('Make sure to set dls correctly')

lmax_autocorr = 10;               % Maximum autocorr. length
C1 = autocorrelation(u1, params(1).dl, lmax_autocorr);
C2 = autocorrelation(u2, params(2).dl, lmax_autocorr);
C3 = autocorrelation(u3, params(3).dl, lmax_autocorr);
C1 = C1 / C1(1);
C2 = C2 / C2(1);
C3 = C3 / C3(1);

warning('Make sure to correctly normalize the autocorrelation !');

% Length vectors
dl1 = params(1).U / params(1).f
dl2 = params(2).U / params(2).f
dl3 = params(3).U / params(3).f


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


% 5) print them for your table
fprintf('\nCorrelation lengths L_C (m):\n');
fprintf('  Dataset 1: %.4f\n', params(1).L_C);
fprintf('  Dataset 2: %.4f\n', params(2).L_C);
fprintf('  Dataset 3: %.4f\n\n', params(3).L_C);

%% --- Plot B ---ok 
lmax_plot = 1.1;
lmax_idx1 = find(l1 >= lmax_plot, 1);
lmax_idx2 = find(l2 >= lmax_plot, 1);
lmax_idx3 = find(l3 >= lmax_plot, 1);
if isempty(lmax_idx1)
    lmax_idx1 = length(l1);
end
if isempty(lmax_idx2)
    lmax_idx2 = length(l2);
end
if isempty(lmax_idx3)
    lmax_idx3 = length(l3);
end

figure();
hold on;
plot(l1(1:lmax_idx1), C1(1:lmax_idx1))
plot(l2(1:lmax_idx2), C2(1:lmax_idx2))
plot(l3(1:lmax_idx3), C3(1:lmax_idx3))
plot([0, l1(lmax_idx1)], 1 / exp(1) * [1, 1], ':k')
grid on;

xlabel('$l$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$C_{i}(l)$', 'Interpreter', 'latex', 'FontSize', 14)

% Add LCi to plots
for i = 1:1:3
    xline(params(i).L_C, '--', ['$L_{C' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end

legend({'$C_1(l)$', '$C_2(l)$', '$C_3(l)$', '$e^{-1}$'}, ...
    'Interpreter', 'latex', 'FontSize', 14)


%% --- Integral scale ---
lmax_plot = 5;
lmax_idx1 = find(l1 >= lmax_plot, 1);
lmax_idx2 = find(l2 >= lmax_plot, 1);
lmax_idx3 = find(l3 >= lmax_plot, 1);
if isempty(lmax_idx1)
    lmax_idx1 = length(l1);
end
if isempty(lmax_idx2)
    lmax_idx2 = length(l2);
end
if isempty(lmax_idx3)
    lmax_idx3 = length(l3);
end

% 1.2.3 Cumulative integral
Lint_cum1 = cumtrapz(l1(1:lmax_idx1), C1(1:lmax_idx1));
Lint_cum2 = cumtrapz(l2(1:lmax_idx2), C2(1:lmax_idx2));
Lint_cum3 = cumtrapz(l3(1:lmax_idx3), C3(1:lmax_idx3));

params(1).Lint = Lint_cum1(end);
params(2).Lint = Lint_cum2(end);
params(3).Lint = Lint_cum3(end);

fprintf('\nIntegral scales L_{int} (m):\n');
fprintf('  Dataset 1: %.4f\n', params(1).Lint);
fprintf('  Dataset 2: %.4f\n', params(2).Lint);
fprintf('  Dataset 3: %.4f\n\n', params(3).Lint);

warning('Make sure to fill the cumulative integral properly')

% Plotting
figure();
hold on; grid on;
plot(l1(1:lmax_idx1), Lint_cum1)
plot(l2(1:lmax_idx2), Lint_cum2)
plot(l3(1:lmax_idx3), Lint_cum3)

xlabel('$l$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\int_0^l  \,C_i (s) \,ds $', 'Interpreter', 'latex', ...
        'FontSize', 14)

for i = 1:1:3
    yline(params(i).Lint, '--', ['$L_{int,' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end

legend({'$\int_0^lC_1(s)ds$', '$\int_0^lC_2(s)ds$', '$\int_0^lC_3(s)ds$', ...
    }, ...
    'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 14)
