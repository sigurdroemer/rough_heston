function [price, integrand] = NumericalIntegrationCall(K,phi,beta,midpoint,...
                                                 ubound,dA,eps)
% Description: Computes the price of a European call option using 
% numerical integration. We assume zero interest rates and dividends and  
% that S(0) = 1 where S(t) is the time t price of the underlying asset.
% 
% Parameters:
%   call:     [1x1 boolean] true (call option), false (put option).
%   K:        [1x1 real] Strike price.
%   T:        [1x1 real] Expiration.
%   phi:      [function] Characteristic function of log(S(T)). Must be 
%             vectorized.
%   beta:     [1x1 real] Parameter for shifting the contour of the pricing 
%             integral, same notation as in (Gerhold et al., 2019). It is 
%             the user's responsibility to ensure that the beta'th moment 
%             of the underlying asset exists so the computations are valid.
%   midpoint: [1x1 real] Midpoint of multi-domain integration, see function
%             'MultiDomainIntegration' and (Zhu, 2010).
%   ubound:   [1x1 real] Upper truncation of integral.
%   dA:       [1x1 real] Sub-interval step size for multi-domain 
%             integration, see function 'MultiDomainIntegration' and 
%             (Zhu, 2010). Will be rounded down to nearest value s.t. 
%             midpoint - ubound is divisible by dA.
%   eps:      [1x1 real] Tolerance of multi-domain integration, see 
%             function 'MultiDomainIntegration' and (Zhu, 2010).
%
% Output: 
%   price:      [1x1] Price of call option.
%   integrand:  [function] Function to be integrated. Added for inspection.
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
%   - Jianwei Zhu, Applications of Fourier Transform to Smile Modelling,
%   Springer, 2010.
%

    s0 = 1;
    k = log(K);
    alpha = beta - 1;
    psi = @(v) (  phi( v - (alpha + 1)*1i  ) ./ ...
        ( alpha^2 + alpha - v.^2 + 1i*(2*alpha + 1).*v)  );
    integrand = @(v)(real(psi(v).*exp(-1i*v*k')).*(exp(-alpha*k')./pi));
    integrand = @(v)(integrand(v')');
    
    % Remark: We use multi-domain integration to avoid evaluating the 
    % characteristic function at unnecessarily extreme values which  
    % can cause NaNs to appear for some implementations:
    
    N = ceil((ubound - midpoint)/dA);
    price = MultiDomainIntegration(integrand,0,midpoint,ubound,N,eps);
    
    % Correct for residue depending on alpha - see e.g. equation (7) 
    % appearing in (Lord & Kahl, 2006):
    R = 0;
    if alpha <= 0
        R = R + s0;
    end
    if alpha <= -1
        R = R - K;
    end
    if alpha==0
        R = R - 0.5*s0;
    end
    if alpha==-1
        R = R + 0.5*K;
    end
    price = price + R;
    
end

