function u = GetLowerCriticalMomentRoughHeston(lambda,rho,xi,T,T_eval)
% Description: Computes the lower critical moment for the rough Heston 
% model as presented in (Gerhold et al, 2019). Important: See 'Remarks'.
%
% The lower critical moment is defined as 
%
%   u(T) := inf( u real: E[exp(ulog(S(T)/S(0)))] < infinity )
%
% where S(t) is the price of the underlying asset at a time point t.
%
% For the computations we use the moment explosion time function T_star(u) 
% defined as
%
%   T_star(u) := sup(t >= 0: E[S(t)^u] < infinity).
%
% which we estimate using 'Algorithm 1' and 'Algorithm 2' from 
% (Gerhold et al., 2019).
% 
% To find the critical moment for a fixed T > 0 we then use numerical 
% root-finding to solve
%
%   T_star(u) = T
%
% w.r.t. u < 0.
%
% Remarks:
%   - The validity of the computations assume that T_star(u) (where finite) 
%   is strictly increasing for u < 0.
%   - The estimated explosion times may be biased low, see 'Algorithm 2', 
%   and therefore the lower critical moment may be biased high.
%
% Parameters:
%   lambda: [1x1 real] See (Gerhold et al., 2019).
%   rho:    [1x1 real] See (Gerhold et al., 2019). We assume rho < 0.
%   xi:     [1x1 real] See (Gerhold et al., 2019).
%   T:      [function] Explosion time function.
%   T_eval: [Nx1 real] Time point(s) to get lower critical moment(s) for.
%   
% Output:
%   u: [Nx1 real] Lower critical moment(s).
%
% References:
%   (1) Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment
%   explosions in the rough Heston model, Decisions in Economics and
%   Finance (2019) 42:575-608.

% Vectorization by recursive calls:
if numel(T_eval) > 1
    u = NaN(size(T_eval));
    for i=1:numel(T_eval)
        if T(i) > 0
            u(i) = GetLowerCriticalMomentRoughHeston(lambda,rho,xi,...
                                                     T,T_eval(i));
        elseif T(i) == 0
            u(i) = Inf;
        else
            error('GetLowerCriticalMomentRoughHeston: Unsupported input.');
        end
    end
    return;
end

% Determine if critical moment is in (-infinity, lambda/(rho*xi)], i.e.
% case A, or not:
u_sep = lambda / (rho*xi);
case_A = false;
if T(u_sep) >= T_eval
    case_A = true;
end

% Find critical moment:
min_fun = @(u)(abs(T(u) - T_eval));
options = optimset('Display','off','MaxIter',10^4,'MaxFunEvals',10^4);
if case_A
    [u,~,exitflag] = fmincon(min_fun,u_sep,[],[],[],[],[],u_sep,[],...
                             options);
else
    [u,~,exitflag] = fmincon(min_fun,u_sep,[],[],[],[],u_sep+eps,0,[],...
                              options);
end

if exitflag <= 0
    error(['GetLowerCriticalMomentRoughHeston: fmincon did not ', ...
           'converge. Please inspect inputs and/or code.']);
end

end

