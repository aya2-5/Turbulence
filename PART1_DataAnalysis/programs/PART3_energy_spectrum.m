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

dl1 = U1 / sample_freq1;
dl2 = U2 / sample_freq2;
dl3 = U3 / sample_freq3;

% Length vectors, from results of questions 1.1 (need to be out back here)
N = length(u1);     % Same for all datasets
l1 = (1:1:N)*dl1;         % params(1).x
l2 = (1:1:N)*dl2;         % params(2).x
l3 = (1:1:N)*dl3;         % params(3).x

L1 = l1(N);
L2 = l2(N);
L3 = l3(N);


%% --- Energy spectrum ---
dx1 = dl1;   dx2 = dl2;   dx3 = dl3;    % grid spacings
L1  = N*dx1; L2  = N*dx2; L3  = N*dx3;  % domain lengths
fprintf('L1 = %.3e m,  L2 = %.3e m,  L3 = %.3e m\n', L1, L2, L3);

norm1 = dx1^2/(2*pi*L1);
norm2 = dx2^2/(2*pi*L2);
norm3 = dx3^2/(2*pi*L3);


fprintf('dx1 = %.3e m, dx2 = %.3e m, dx3 = %.3e m\n', dx1, dx2, dx3);
fprintf('Prefactor1 = %.3e  Prefactor2 = %.3e  Prefactor3 = %.3e\n', ...
        norm1, norm2, norm3);
warning('Make sure to have properly written the normalization')


% Define E(k), k > 0
Ek1 = spectral_energy_density(u1)*norm1;
Ek2 = spectral_energy_density(u2)*norm2;
Ek3 = spectral_energy_density(u3)*norm3;

M  = N/2 + 1;   % number of nonnegative modes
% Nyquistfreq1= sample_freq1/2;
% Nyquistfreq2= sample_freq2/2;
% Nyquistfreq3= sample_freq3/2;

k1 = (0:M-1)*(2*pi/L1);
k2 = (0:M-1)*(2*pi/L2);
k3 = (0:M-1)*(2*pi/L3);

Ek1_pos = Ek1(1:M);
Ek2_pos = Ek2(1:M);
Ek3_pos = Ek3(1:M);

% form the one‐sided E(k)
E1 = Ek1_pos;   E2 = Ek2_pos;   E3 = Ek3_pos;
warning('Define E(k), k > 0 properly, and the k vectors as well')

% Signal filtering
E1(2:M-1) = 2*Ek1_pos(2:M-1);
E2(2:M-1) = 2*Ek2_pos(2:M-1);
E3(2:M-1) = 2*Ek3_pos(2:M-1);



% Define normalized cutoff frequency (0 < Fc < 1).
% Fc = 0.1 means smoothing out fluctuations above 10% of Nyquist.
Fc = 0.1;
Ek1_smooth = movmean(E1, 50);
Ek2_smooth = movmean(E2, 50); 
Ek3_smooth = movmean(E3, 50);
warning('Make sure that you filtered the spectrums')

% time–domain kinetic energy:
E_td1 = 0.5*mean(u1.^2);

% spectral integral JUST TO CHECK:
E_fd1 = trapz(k1, E1);  
rel_err1 = abs(E_td1 - E_fd1)/E_td1;
fprintf("Parseval relative error for data1: %.2f%%\n", rel_err1*100);

%% Plot C
figure()
loglog(k1(2:end), Ek1_smooth(2:end), 'DisplayName','Data 1','LineWidth',1.2);
hold on;
grid on;
loglog(k2(2:end), Ek2_smooth(2:end), 'DisplayName','Data 2','LineWidth',1.2);
loglog(k3(2:end), Ek3_smooth(2:end), 'DisplayName','Data 3','LineWidth',1.2);
yl = ylim;          % get current [ymin ymax]
ylim([yl(1) 1e4]);  % keep the same bottom, set top = 10^4


% --- Add K41 k^(-5/3) reference line ---
k_fit = logspace(log10(k1(2)), log10(k1(end)), 200);
C     = mean(Ek1_smooth(k1>1 & k1<2));   % anchor in the inertial band
hRef  = loglog(k_fit, C*k_fit.^(-5/3), 'k--', ...
               'LineWidth',1, 'DisplayName','K41: $k^{-5/3}$');

% --- Add a text label on the figure ---
x_text = 2; 
y_text = C * x_text^(-5/3) * 1.3;   % a little above the line
text(x_text, y_text, '$k^{-5/3}$', ...
     'Interpreter','latex', 'FontSize',15, 'Color','k', ...
     'HorizontalAlignment','center');


% ————— Re‐compute η_E by 10%‐peak rule —————
peak1 = max(Ek1_smooth); thresh1 = 0.1*peak1;
i_eta1 = find(Ek1_smooth < thresh1,1,'last'); 
k_eta1 = k1(i_eta1);
etaE1 = 2*pi / k_eta1;

peak2 = max(Ek2_smooth); thresh2 = 0.1*peak2;
i_eta2 = find(Ek2_smooth < thresh2,1,'last'); 
k_eta2 = k2(i_eta2);
etaE2 = 2*pi / k_eta2;

peak3 = max(Ek3_smooth); thresh3 = 0.1*peak3;
i_eta3 = find(Ek3_smooth < thresh3,1,'last'); 
k_eta3 = k3(i_eta3);
etaE3 = 2*pi / k_eta3;




%Add estimates for integral and Kolmogorov length scales
%---- Integral scale from the peak of E(k) ----
% drop the k=0 bin to avoid division by zero
k_nz1 = k1(2:end);     E_nz1 = E1(2:end);
k_nz2 = k2(2:end);     E_nz2 = E2(2:end);
k_nz3 = k3(2:end);     E_nz3 = E3(2:end);

% compute L_int by the standard formula
L_intE1 = trapz(k_nz1, E_nz1./k_nz1) / trapz(k_nz1, E_nz1);
L_intE2 = trapz(k_nz2, E_nz2./k_nz2) / trapz(k_nz2, E_nz2);
L_intE3 = trapz(k_nz3, E_nz3./k_nz3) / trapz(k_nz3, E_nz3);

% ---- store in params for marking ----
params(1).Lint_e = L_intE1;   params(2).Lint_e = L_intE2;   params(3).Lint_e = L_intE3;
params(1).eta_e  = etaE1;     params(2).eta_e  = etaE2;     params(3).eta_e  = etaE3;


% ---- print the estimates ----
fprintf('Estimated integral scales [m]:   %.3f, %.3f, %.3f\n', L_intE1, L_intE2, L_intE3);
fprintf('Estimated Kolmogorov scales [m]: %.6f, %.6f, %.6f\n', etaE1, etaE2, etaE3);


warning(['Make sure to add the estimates for the integral and ' ...
    'Kolmogorov length scales'])

% Labels
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$E(k)$', 'Interpreter', 'latex')
legend('Interpreter', 'latex')

% Plot correlation length scales Lci
for i = 1:1:3
    xline(2*pi/params(i).eta_e, '--', ['$\eta_{E' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    xline(2*pi/params(i).Lint_e, '--', ['$L_{int,E' num2str(i) '}$'], ...
        'Interpreter', 'latex', 'fontsize', 10, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end
% xlim([5e-2, 1e4])           % You can ajust this if you need it



%% 1.3.3 – Parseval’s theorem check 
% 1) time‐domain energy
E_td1 = 0.5*mean(u1.^2);
E_td2 = 0.5*mean(u2.^2);
E_td3 = 0.5*mean(u3.^2);

fprintf('dx1 = %.3e m, dx2 = %.3e m, dx3 = %.3e m\n', E_td1, E_td2, E_td3);
% 2) spectral‐domain energy (raw, folded E1)
E_fd1 = trapz(k1, E1);
E_fd2 = trapz(k2, E2);
E_fd3 = trapz(k3, E3);


% spectral‐domain (smoothed) energies
E_fd1_sm = trapz(k1, Ek1_smooth);
E_fd2_sm = trapz(k2, Ek2_smooth);
E_fd3_sm = trapz(k3, Ek3_smooth);


err1_raw  = abs(E_td1  - E_fd1 ) / E_td1  *100;
err1_smoo = abs(E_td1  - E_fd1_sm) / E_td1  *100;

err2_raw  = abs(E_td2  - E_fd2 ) / E_td2  *100;
err2_smoo = abs(E_td2  - E_fd2_sm) / E_td2  *100;

err3_raw  = abs(E_td3  - E_fd3 ) / E_td3  *100;
err3_smoo = abs(E_td3  - E_fd3_sm) / E_td3  *100;

% ————— Print a little results table —————
fprintf('\n     Dataset |   L_int [m]   |   eta [m]    |  Err_raw(%%) | Err_smoo(%%)\n');
fprintf('  ---------------------------------------------------------------\n');
fprintf('      1     |  %8.4e   |  %8.4e   |   %6.3f   |   %.3e\n', ...
        L_intE1, etaE1, err1_raw, err1_smoo);
fprintf('      2     |  %8.4e   |  %8.4e   |   %6.3f   |   %6.3f\n', ...
        L_intE2, etaE2, err2_raw, err2_smoo);
fprintf('      3     |  %8.4e   |  %8.4e   |   %6.3f   |   %6.3f\n\n', ...
        L_intE3, etaE3, err3_raw, err3_smoo);