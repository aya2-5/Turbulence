function [x, d]=increment(u, l, f, U)
    % u: vector of perturbative velocity [m/s]
    % l: length for the increment calculation [m]
    % f: frequency of data sampling [Hz]
    % U: mean flow speed [m/s]
    % x: vector of lengths [m]
    % d: vector of velocity increments [m/s]

    % spatial sampling interval
    dx = U / f;
    % number of samples corresponding to separation l
    Nlag = round(l / dx);

    % ensure lag is not larger than signal
    N = length(u);
    if Nlag >= N
        error('Specified separation l is too large for input signal length.');
    end

    % compute increments: u(i+Nlag) - u(i)
    M = N - Nlag;
    d = u((1+Nlag):end) - u(1:M);

    % spatial positions corresponding to each increment
    x = (0:(M-1))' * dx;
    
end