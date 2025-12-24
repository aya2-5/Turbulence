function S = structure_function(u, l, f, U, n)
    % u: vector of perturbative velocity [m/s]
    % l: length for the increment calculation [m]
    % f: frequency of data sampling [Hz]
    % U: mean flow speed [m/s]
    % n : order of the structure function
    % S : vector of the structure function [m/s]
   
 % ensure column vector
    u = u(:);

    % number of separation values
    Nl = numel(l);
    S   = nan(size(l));

    for i = 1:Nl
        % 1) convert spatial lag to time lag via Taylor’s hypothesis:
        %    τ = ℓ / U  [s]
        tau = l(i) / U;

        % 2) convert time lag to sample lag (number of points)
        d = round(tau * f);

        % skip if too small or beyond data length
        if d < 1 || d >= length(u)
            S(i) = NaN;
            continue
        end

        % 3) build the velocity increment series δuₗ = u(t+d) – u(t)
        du = u(1+d:end) - u(1:end-d);

        % 4) compute the n-th order moment (structure function)
        S(i) = mean( du .^ n );
    end
end