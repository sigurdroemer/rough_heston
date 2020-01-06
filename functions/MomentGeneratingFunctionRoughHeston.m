function phi = MomentGeneratingFunctionRoughHeston(v_0,alpha,lambda,...
                                                   v_bar,xi,rho,T,u,N)
% Description: Computes the moment generating function of log(S(t)/S(0)) 
% where S(t) is the price of the underlying asset at a time point t and
% where everything is under the rough Heston model as presented 
% in e.g. (Gerhold et al., 2019).
% 
% Using the same notation as in (Gerhold et al., 2019) we assume
% 1/2 < alpha <= 1 and V(0), lambda, v_bar, xi > 0 and -1 <= rho < 0.
%
% The code returns 
%
%   phi(u) := E [ exp(u*log(S(T)/S(0)) ]                              (*)
%
% for an input u where we are using the formulation
%
%   phi(u) = exp ( v_bar*lambda*I_t^1 psi(u,t) 
%                  + v_0*I_t^(1-alpha)psi(u,t))                       (**)
%
% for the actual computation and where the mapping t -> psi(u,t) 
% solves a Volterra integral equation. 
% 
% We refer to (Gerhold et al., 2019) for the details and the notation used.
%
% Parameters:
%   v_0:    [1x1 real] See (Gerhold et al., 2019).
%   alpha:  [1x1 real] See (Gerhold et al., 2019).
%   lambda: [1x1 real] See (Gerhold et al., 2019).
%   v_bar:  [1x1 real] See (Gerhold et al., 2019).
%   xi:     [1x1 real] See (Gerhold et al., 2019).
%   rho:    [1x1 real] See (Gerhold et al., 2019).
%   T:      [1x1 real] Time point to evaluate characteristic function at.
%   u:      [Lx1 real or complex] Point to evaluate (*) at. Note that we 
%           allow complex values s.t. we can also return the characteristic 
%           function.
%   N:      [1x1 integer] Number of steps when solving the Volterra 
%           integral equation.
% 
% Output:
%   phi: [Lx1 complex or real] See (*).
% 
% References: 
%   - Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment
%   explosions in the rough Heston model, Decisions in Economics and
%   Finance (2019) 42:575-608.
% 
%   - Kai Diethelm, Neville J. Ford and Alan D. Freed, Detailed error
%   analysis for a fractional ADAMs method, Numerical Algorithms 36:31-52,
%   2004.
% 

    % Define the Volterra integral equation:
    c1 = 0.5*u.*(u-1);
    c2 =  rho*xi*u - lambda;
    c3 = 0.5*xi.^2;
    R = @(w) (c1 + c2.*w + c3.*w.^2);
    f = @(s,y) (R(y));

    % Solve it:
    [psi,Dalpha_psi] = SolveVIE(f,alpha,T,N);
    
    % Integrate to get the characteristic function, see (**):
    dt = T / N;
    phi = exp( v_bar*lambda*sum(psi(:,1:end-1),2).*dt ...
               + v_0.*sum(Dalpha_psi(:,1:end-1),2).*dt );
    
end

