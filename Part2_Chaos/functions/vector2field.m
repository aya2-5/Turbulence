%%% DESCRIPTION -----------------------------------------------------------
%   1D periodic field corresponding to a minimal vector of state variables


%%% INPUTS ----------------------------------------------------------------
%   v       state vector (column vector of real numbers)
%   N       spatial resolution
%   symm    center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   u       the field in physical state (column vector of real numbers)


%%% REMARKS ---------------------------------------------------------------
%   1-  In forward integration of the KSE, the mean remains constant.
%       We assume the mean value is always zero, hence the mean value of
%       'u' is zero by construction.
%   2-  We do not fill the padded elements in Fourier space, hence no
%       additional dealiasing is needed.


function u = vector2field(v,N,symm)
    Nd = round(N/3);

    U = zeros(N,1);
    U(2:Nd+1) = v(1:Nd)*1j;
    
    if ~symm
        U(2:Nd+1) = U(2:Nd+1) + v(Nd+1:2*Nd);
    end

    U(N-Nd+1:N) = flip(conj(U(2:Nd+1)));
    
    u = ifft(U,'symmetric');
end