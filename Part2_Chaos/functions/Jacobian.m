%%% DESCRIPTION -----------------------------------------------------------
%   Jacobian of the KSE flow map: J=dv(t)/dv(0)


%%% INPUTS ----------------------------------------------------------------
%   v0          reference point of the Jacobian (column state vector)
%   t           integration time interval
%   epsilon     perturbation magnitude for finite difference derivatrives
%   dt          step size in time integrations
%   L           domain length
%   N           spatial resolution
%   symm        center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   J           Jacobian matrix

function J = Jacobian(v0, t, epsilon, dt, L, N, symm)
% Jacobian of the KSE flow map: J = dv(t)/dv(0)
%
% INPUTS:
%   v0       – reference initial state (n×1)
%   t        – integration time
%   epsilon  – finite‐difference perturbation size
%   dt       – time step for KSE_integrate
%   L, N     – domain length and resolution
%   symm     – center‐symmetry flag
%
% OUTPUT:
%   J        – n×n Jacobian matrix

    n = length(v0);
    J = zeros(n,n);

    % for each component j, do central difference
    for j = 1:n
        e = zeros(n,1);
        e(j) = 1;
        % integrate forward from v0 + ε e_j
        [v_plus, ~]  = KSE_integrate(v0 + epsilon*e, t, dt, 0, L, N, symm);
        % integrate forward from v0 - ε e_j
        [v_minus, ~] = KSE_integrate(v0 - epsilon*e, t, dt, 0, L, N, symm);
        % central‐difference quotient
        J(:,j) = (v_plus - v_minus) / (2*epsilon);
    end
end
