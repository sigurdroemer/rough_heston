function [price, iv] = NumericalIntegrationRoughHeston(s_0,v_0,alpha,lambda,...
                                                v_bar,xi,rho,call,K,T,...
                                                varargin)
% Description: Computes prices of European call and put options under the
% rough Heston model as presented in e.g. (Gerhold et al., 2019) or (El
% Euch & Rosenbaum, 2018) and (El Euch & Rosenbaum 2019).
%
% The model for the asset price S(t) can be stated as
%
%   dS(t) = S(t)*sqrt(V(t))*dW_1(t)
%   V(t) = V(0) + int_0^t K(t-s) lambda*(v_bar - V(s))ds 
%               + int_0^t K(t-s) xi * sqrt(V(s))dW_2(s)
%
% with dW_1(t)dW_2(t) = rho dt and K(t) = (1/gamma(alpha))*t^(alpha-1).
% 
% We assume 1/2 < alpha <= 1 and V(0), lambda, v_bar, xi > 0 and 
% -1 <= rho < 0.
%
% To compute prices we use numerical integration via a Fourier transform
% as presented in e.g. (Gerhold et al., 2019) with an optimal integration 
% contour as proposed in (Lord & Kahl, 2006). To compute the characteristic
% function for the integration we solve a Volterra integral equation which 
% appears in e.g. (Gerhold et al., 2019) using the scheme from 
% (Diethelm, 2004).
%
% Parameters:
%   s_0:        [1x1 real] Initial value of asset price.
%   v_0:        [1x1 real] Initial value of variance process.
%   alpha:      [1x1 real] See equations above.
%   lambda:     [1x1 real] See equations above.
%   v_bar:      [1x1 real] See equations above.
%   xi:         [1x1 real] See equations above.
%   rho:        [1x1 real] See equations above.
%   call:       [1x1 boolean] true (call option), false (put option).
%   K:          [Nx1 real] Strike price(s).
%   T:          [Mx1 real] Expiration(s).
%   r:          [1x1 real] Risk-free interest rate.
%   q:          [1x1 real] Dividend yield.
%   beta:       [1x1 real, optional (default = see text)] Parameter for 
%               shifting the integration contour ideally to lower the 
%               oscillations of the integrand. Must be s.t. the beta'th 
%               moment of S(T) is finite (user responsibility). Default is 
%               an optimal value computed along the lines of 
%               (Lord & Kahl, 2006).
%   N:          [1x1 integer, optional(default=250)] Number of steps used 
%               to solve the Volterra integral equation.
%   ubound:     [1x1 real, optional(default=2000)] Upper truncation of 
%               Fourier integral.
%   midpoint:   [1x1 real, optional(default=10)] Used for multi-domain 
%               integration, see function 'NumericalIntegrationCall' and 
%               (Zhu, 2010).
%   dA:         [1x1 real, optional(default=50)] Same as above.
%   eps:        [1x1 real, optional(default=10^(-6))] Same as above.
%   nmax:       [1x1 integer, optional(default=200)] Number of terms in 
%               recursive scheme to compute moment explosion times which 
%               are needed when searching for an optimal contour 
%               shift, i.e. an optimal 'beta'.
%   disp_iter:  [1x1 boolean, optional(false)] If true and multiple prices 
%               are requested we show the time elapsed of each iteration.
% 
% Output: 
%   price: [NxM real] Prices of call or put options.
%   iv:    [NxM real] Black-scholes implied volatilities.
%
% Example: [price, iv] = NumericalIntegrationRoughHeston(100,0.2^2,0.6,2,...
%                                                  0.2^2,0.2,-0.8,true,...
%                                                  [90;100;110],[0.5;1])
%
% References:
%   - Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment
%   explosions in the rough Heston model, Decisions in Economics and
%   Finance (2019) 42:575-608.
%
%   - Omar El Euch and Mathieu Rosenbaum, Perfect hedging in rough Heston
%   models, Mathematical Finance 28(6), 3813-3856, 2018.
%
%   - Omar El Euch and Mathieu Rosenbaum, The characteristic function of
%   rough Heston models, Mathematical Finance 29(1), 3-38, 2019.
%
%   - Roger lord and Christian Kahl, Optimal Fourier inversion in 
%   semi-analytical option pricing, Tinbergen Institute Discussion Paper 
%   No. 2006-066/2. Available at SRRN:
%   https://papers.ssrn.com/sol3/papers.cfm?abstract_id=921336.
% 
%   - Kai Diethelm, Neville J. Ford and Alan D. Freed, Detailed error
%   analysis for a fractional ADAMs method, Numerical Algorithms 36:31-52,
%   2004.
%

% Parse and check inputs:
p = inputParser;
addOptional(p,'r',0);
addOptional(p,'q',0);
addOptional(p,'beta',[]);
addOptional(p,'N',250);
addOptional(p,'ubound',2000);
addOptional(p,'midpoint',10);
addOptional(p,'dA',50);
addOptional(p,'eps',10^(-6));
addOptional(p,'nmax',200);
addOptional(p,'disp_iter',false);
parse(p,varargin{:});
r = p.Results.r;
q = p.Results.q;
beta = p.Results.beta;
N = p.Results.N;
ubound = p.Results.ubound;
midpoint = p.Results.midpoint;
dA = p.Results.dA;
eps = p.Results.eps;
nmax = p.Results.nmax;
disp_iter = p.Results.disp_iter;

if rho >= 0
    error('NumericalIntegrationRoughHeston: rho >= 0 is current not supported.');
end

if alpha > 1 || alpha <= 1/2
    error('NumericalIntegrationRoughHeston: Alpha must lie in (1/2,1].');
end

if any([v_0,lambda,v_bar,xi]) <= 0
    error(['NumericalIntegrationRoughHeston: v_0, lambda, v_bar and xi ',...
           'must be strictly positive.']);
end


if size(K,2) > 1 || size(T,2) > 1
    error(['NumericalIntegrationRoughHeston: Strike and expiration vectors ', ...
        'must be column vectors.']);
end

% Define characteristic function:
phi = @(u)(MomentGeneratingFunctionRoughHeston(v_0,alpha,lambda,v_bar,...
                                               xi,rho,T,1i*u,N));

% Vectorize by recursive calls:
tot_iter = numel(K)*numel(T);
if tot_iter > 1
    [price,iv] = deal(NaN(size(K,1),size(T,1)));
    ii = 1;
    for i=1:numel(K)
        for j=1:numel(T)
            tic;
            [price(i,j),iv(i,j)] = NumericalIntegrationRoughHeston(s_0,...
                               v_0,alpha,lambda,v_bar,xi,rho,call,K(i),...
                               T(j),'r',r,'q',q,'beta',beta,'N',N,...
                               'ubound',ubound,'nmax',nmax,...
                               'midpoint',midpoint,'eps',eps,'dA',dA);
            tt = toc;
            if disp_iter
                disp([num2str(ii), ' / ', num2str(tot_iter), ...
                   ' prices computed (time: ',num2str(tt),' sec.)']);
            end
            ii = ii + 1;
        end
    end
    return;
end

F = s_0*exp((r-q)*T);
Kadj = K ./ F;
k = log(Kadj);

% Compute optimal beta:
if ~exist('beta','var') || isempty(beta)
    T_star = @(u)(MomentExplosionTimeRoughHeston(alpha,lambda,xi,rho,u,...
                                                 nmax)');
    % Remark: Estimated lower critical moment may be biased high:
    lower_critical = GetLowerCriticalMomentRoughHeston(lambda,rho,xi,...
                                                       T_star,T);
    beta = GetOptimalFourierAlpha(phi,k,lower_critical) + 1;
end

% Compute price of call option assuming S(0) = 1, zero interest rates,
% zero dividends and using the adjusted strike:
price = NumericalIntegrationCall(Kadj,phi,beta,midpoint,ubound,dA,eps);

% Convert to put option if needed:
if ~call
    price = price - 1 + Kadj;
end

if price < 0
    error('NumericalIntegrationRoughHeston: Negative price computed. ');
end

% Scale and adjust back:
price = s_0*exp(-q*T).*price;

% Compute implied volatilities:
iv = blsimpv(s_0,K,r,T,price,'Yield',q,'class',call);

end

