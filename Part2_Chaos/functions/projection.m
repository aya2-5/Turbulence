%%% DESCRIPTION -----------------------------------------------------------
%   projecting a set of state vectors onto their energy content, energy
%   production rate and energy dissipation rate


%%% INPUTS ----------------------------------------------------------------
%   v       a matrix whose columns are state vectors
%   N       spatial resolution
%   L       domain length
%   symm    center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   E       vec. whose i-th element is energy of the i-th column of v
%   P       vec. whose i-th element is production of the i-th column of v
%   D       vec. whose i-th element is dissipation of the i-th column of v

function [E,P,D] = projection(v,N,L,symm)
    M = size(v,2);
    
    E = zeros(M,1);
    P = zeros(M,1);
    D = zeros(M,1);
    
	[~,k] = domain(L,N);
    
    for i=1:M
        u = vector2field(v(:,i),N,symm);
        U = fft(u);
        
        du_dx = ifft((1j*k).*U,'symmetric');
        d2u_dx2 = ifft((-k.^2).*U,'symmetric');
    
        E(i) = 0.5*L*mean(u.^2);
        P(i) = L*mean(du_dx.^2);
        D(i) = L*mean(d2u_dx2.^2);
    end
end