%% ========================================================================
%  7_2_dissipation_estimate.m
%
%  This script estimates the energy dissipation rate (ε) from the second-
%  and third-order structure functions, S₂(l) and S₃(l), using the three
%  given velocity datasets.

% The estimation formulae are:
% - From S₃(l):  S₃(l) ≈ 4/5 * ε * l      (valid in the inertial range)
% - From S₂(l):  S₂(l) ≈ C₂ * ε^(3/2) * l^(3/2)   with C₂ ≈ 2.1 (Sreenivasan)
%
%  Complete the missing parts below to compute the dissipation rate.
%  Compare your estimates to those based on integral-scale quantities 
%  discussed in Section 1.4.
% ========================================================================

clear
clc

addpath("../functions/");

files = {'veldata1.txt', 'veldata2.txt', 'veldata3.txt'};
Nsamples = 100000;  % Define the number of samples to read

%for i = 1:length(files)
    % --- Load and preprocess data ---
    % Call a function to load and preprocess the velocity data
    % [u, f, U] = load_data( ??? );

    % --- Compute Structure Functions ---
    % Compute the second-order structure function: 
    % S2 = structure_function( ??? ); 

    % Compute the third-order structure function:
    % S3 = structure_function( ??? );

    % --- Estimate ε from S3 ---
    
    % --- Estimate ε from S2 ---
    
    % --- Display the results ---
    % Comparison of the estimates
%end

%% ========================================================================
%  7_2_dissipation_estimate.m
%  Estimates ε from S2 and S3 for three datasets.
% ========================================================================


% --- Structure‐function grids ---
Num   = 100;
l_S2  = logspace(log10(1e-3), log10(10), Num)';   % same as Part 7.1
l_S3  = logspace(log10(1e-3), log10(1),  Num)';

% --- known constants & scales ---
C2      = 2.1;                                  % Sreenivasan prefactor
lambda  = [4.892e-3, 3.885e-3, 3.687e-3];       % Taylor scales
L_C     = [0.2807,   0.3402,   0.4002  ];       % integral scales

% preallocate
epsilon2 = zeros(1,3);
epsilon3 = zeros(1,3);
epsilon_int = zeros(1,3);

fprintf('\n   Dataset    ε_{S3} [m^2/s^3]    ε_{S2} [m^2/s^3]\n');
fprintf('  ------------------------------------------------\n');

for i = 1:3
    %--- load data ---
    [u,f,U] = load_data(files{i}, Nsamples);

  
    %--- compute S2 & S3 ---
    S2 = structure_function(u, l_S2, f, U, 2);
    S3 = structure_function(u, l_S3, f, U, 3);
    S2 = smooth(S2,8);   % optional smoothing
    S3 = smooth(S3,3);

    %--- inertial‐range masks ---
    m2 = (l_S2 > 1.1*lambda(i)) & (l_S2 < 0.9*L_C(i));
    m3 = (l_S3 > 1.1*lambda(i)) & (l_S3 < 0.9*L_C(i));

    %--- ε from four‐fifth law:  S3 ≈ (4/5) ε ℓ  =>  ε = (5/4) S3/ℓ ---
    epsilon3(i) = mean( (5/4) * S3(m3) ./ l_S3(m3) );

    %--- ε from second‐order law:  S2 ≈ C2 ε^(2/3) ℓ^(2/3)
    %    => ε = [ (S2/C2)^(3/2) ] / ℓ
    epsilon2(i) = mean( ( (S2(m2)/C2).^(3/2) ) ./ l_S2(m2) );

    % print it 
    fprintf(' Dataset %d: ε_{S3} = %8.3e  ε_{S2} = %8.3e\n', ...
            i, epsilon3(i), epsilon2(i));
end
