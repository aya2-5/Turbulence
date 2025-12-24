%%% DESCRIPTION -----------------------------------------------------------
%   1D spatial domain in physical and Fourier states


%%% INPUTS ----------------------------------------------------------------
%   L   domain length
%   N   spatial resolution


%%% OUTPUTS ---------------------------------------------------------------
%   x   column vector of grid points
%   k   column vector of wave numbers


function [x,k] = domain(L,N)
    dx = L/N;
    x = linspace(0,L-dx,N)' - L/2;
    k = (2*pi/L)*[0:N/2-1 0 -N/2+1:-1]';
end