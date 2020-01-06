function [T_star,case_A,case_B] = MomentExplosionTimeRoughHeston(alpha,...
                                                         lambda,xi,rho,...
                                                         u,n_max,num_stab)
% Description: Computes the moment explosion time of the rough Heston
% model as presented in e.g. (Gerhold et al., 2019). Letting S(t) denote
% the price of the underlying asset at a time point t what we return is an
% estimate of
%
%   T_star(u) := sup(t >= 0: E[S(t)^u] < infinity).
%
% For the computation we use 'Algorithm 1' and 'Algorithm 2' from (Gerhold
% et al., 2019). Since 'Algorithm 2' technically estimates a lower bound
% note that estimates may be biased low.
%
% Parameters:
%   alpha:    [1x1 real] See notation in (Gerhold et al., 2019).
%   lambda:   [1x1 real] See notation in (Gerhold et al., 2019).
%   xi:       [1x1 real] See notation in (Gerhold et al., 2019).
%   rho:      [1x1 real] See notation in (Gerhold et al., 2019).
%   u:        [Nx1 real] Moment(s) to find explosion time(s) for.
%   n_max:    [1x1 integer] Maximum number of steps in recursive scheme 
%              to estimate explosion time(s).
%   num_stab: [1x1 boolean, optional(default=true)] If true we in some 
%             (rare) cases adjust the algorithm to improve the numerical 
%              stability. See the code.
% 
% Output:
%   T_star: [Nx1 real] Moment explosion time(s) (can be infinite).
%   case_A: [Nx1 real] True if case A, see (Gerhold et al., 2019).
%   case_B: [Nx1 real] True if case B, see (Gerhold et al., 2019).
%   
% References: 
%   - Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment
%   explosions in the rough Heston model, Decisions in Economics and
%   Finance (2019) 42:575-608.
%

if ~exist('num_stab','var') || isempty(num_stab);num_stab = true;end

% Define coefficients:
c1 = 0.5*u.*(u-1);
c2 =  rho*xi*u - lambda;
c3 = 0.5*xi.^2;
e0 = 0.5*c2;
e1 = e0.^2 - c3*c1;

% Determine cases:
case_A = c1 > 0 & e0 >= 0;
case_B = c1 > 0 & e0 < 0 & e1 < 0;

% Set infinite explosion times:
T_star = NaN(size(u,1),1);
T_star(~case_A & ~case_B) = Inf;

% Define more coefficients and relevant functions:
v_fun = @(n)(gamma(alpha*n + 1) ./ gamma(alpha*n - alpha + 1));
d1 = c1*c3;
d2 = c2;

% Recursively compute 'a' coefficients:
a = NaN(size(u,1),n_max);
a(:,1) = d1 / v_fun(1);
large_number = 100;
for n=1:n_max-1
    if alpha*n + 1 < large_number
        v = v_fun(n+1);
    else
        % Use asymptotic approximation of ratio to avoid NaN, 
        % is easily derived from Gautschi's inequality:
        v = ((n+1)*alpha).^(alpha);
    end
    if n > 1
        a(:,n+1) = (1/v)*(d2.*a(:,n)+sum(a(:,1:n-1).*flip(a(:,1:n-1),2),2)); 
    else
        a(:,n+1) = (1/v)*(d2.*a(:,n));
    end
end

% Compute explosion times for case A:
dummy_1 = gamma(alpha)^2 / ( (alpha^(alpha))*gamma(2*alpha) );
dummy_2 = -1./(alpha*(n_max+1));
T_star(case_A) = (a(case_A,end).*(n_max.^(1-alpha)).*dummy_1).^(dummy_2);

% Compute explosion times for case B:
T_star(case_B) = abs(a(case_B,end)).^(-1./(alpha*n_max));

if num_stab
    % For increased numerical stability: When abs(e0) is close to zero we 
    % sometimes experience slow convergence and/or numerical instability, 
    % thus we instead interpolate around this point:
    delta = 0.1;
    u_point = (lambda/(rho*xi));
    idxAdj = find(abs(u - u_point) <= delta); 
    if ~isempty(idxAdj)
        [Tadj, case_A, case_B] = MomentExplosionTimeRoughHeston(alpha,...
                lambda,xi,rho,[u_point-delta;u_point+delta],n_max,false);
        if all(case_A | case_B)
            % Remark: The interpolation is disregarded if 
            % u* = lambda/(rho*xi)+delta is not in case B in which case the 
            % explosion time is infinite and interpolation risks 
            % overestimating the explosion time at u.
            w = ( u(idxAdj) - (u_point  - delta) )./(2.*delta);            
            T_star(idxAdj) = Tadj(1).*(1-w) + w.*Tadj(2);
        end
    end
end

end

