%%% DESCRIPTION -----------------------------------------------------------
%   shooting-based Newton iterations for computing equilibrium states (EQs)


%%% INPUTS ----------------------------------------------------------------
%   v0      guessed equilibrium solution (column state vector)
%   T       time interval over which EQ is verified not to change
%   dt      reference time step for the required time marchings
%   L       domain length
%   N       spatial resolution
%   symm    center symmetry (true/false boolean)


%%% OUTPUTS ---------------------------------------------------------------
%   v_best  exact equilibrium solution (column state vector)
%   flag    flag == 1: search successful, otherwise: search failed


function [v_best,flag] = search4EQ(v0,T,dt,L,N,symm)
    options = optimoptions('fsolve','Display','iter','MaxIterations',75);

    [v_best,~,flag,~] = fsolve(@(X) recurrence(X),v0,options);
    
    function f = recurrence(v)
        [vT,~] = KSE_integrate(v,T,dt,0,L,N,symm);
        f = vT - v;
    end
end