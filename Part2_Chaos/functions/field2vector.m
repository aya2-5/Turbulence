%%% DESCRIPTION -----------------------------------------------------------
%   the minimal vector of state variables defining a 1D periodic field


%%% INPUTS ----------------------------------------------------------------
%   u       the field in physical state (column vector of real numbers)
%   N       spatial resolution
%   symm    center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   v       state vector (column vector of real numbers)


%%% REMARKS ---------------------------------------------------------------
%   1-  In forward integration of the KSE, the mean remains constant.
%       We assume the mean value is always zero, hence not stored.


function v = field2vector(u,N,symm)
    Nd = round(N/3);
    
    U = fft(u);
    v = imag(U(2:Nd+1));

    if ~symm
        v = [v;real(U(2:Nd+1))];
    end
end