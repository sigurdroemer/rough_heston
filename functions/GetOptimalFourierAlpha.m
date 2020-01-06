function alpha = GetOptimalFourierAlpha(phi,k,lower_critical_moment)
% Description: Follows (Lord & Kahl, 2006) by finding an optimal contour 
% to compute the pricing integral along when using Fourier inversion 
% methods for options pricing. What is returned is the optimal 'alpha' 
% parameter from that paper.
%
% We assume the initial price of the underlying asset is 1.
%
% Parameters:
%   phi:                   [function] The characteristic function of the 
%                           log asset price.
%   k:                     [1x1 real] Log strike price.
%   lower_critical_moment: [1x1 real] Lower critical moment of the log 
%                           asset price.
%
% Output:
%   alpha: [1x1 real] The optimal alpha value.
%
% References:
%   - Roger lord and Christian Kahl, Optimal Fourier inversion in 
%   semi-analytical option pricing, Tinbergen Institute Discussion Paper 
%   No. 2006-066/2. Available at SRRN:
%   https://papers.ssrn.com/sol3/papers.cfm?abstract_id=921336.
%

psi = @(v,alpha)( real(exp(-1i.*v.*k).*phi( v - 1i*(alpha + 1) ) ...
                  ./ (-(v - 1i.*alpha).*(v - 1i.*(alpha+1)))) );
min_fun = @(alpha)( -alpha*k + 0.5*log( psi(0,alpha).^2 ) );
options = optimset('Display','off','MaxIter',10^4,'MaxFunEvals',10^4);
lb = lower_critical_moment - 1;ub = -1;
[alpha, ~, exitflag] = fmincon(min_fun,(ub-lb)/2,[],[],[],[],lb,ub,[],options);

if exitflag <= 0
    error(['GetOptimalFourierAlpha: fmincon did not converge. Please ', ...
           'inspect inputs and/or code.']);
end

end

