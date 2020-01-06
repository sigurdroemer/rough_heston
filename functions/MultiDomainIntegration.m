function val = MultiDomainIntegration(f,a,b,c,N,eps)
% Description: Implements a method to evaluate the integral
%
%   int_a^c f(x) dx = int_a^b f(x) dx + int_b^c f(x) dx := A_0 + A_1
%
% where a < b <= c and that by first evaluating the integral A_0 using 
% Matlab's own 'integral' function directly and then evaluating the 
% remaining integral A_1 using the multi-domain integration method of 
% (Zhu, 2010) (still with Matlab's 'integral' function).
%
% Parameters:
%   f:   [function] Function to integrate.
%   a:   [1x1 real] Left endpoint of integration interval.
%   b:   [1x1 real] The interval midpoint from where we use the method of
%        (Zhu, 2010).
%   c:   [1x1 real] Right endpoint of integration interval.
%   N:   [1x1 real] Number of sub-intervals to divide [b,c] into.
%   eps: [1x1 real] We stop the integration of [b,c] when the absolute 
%         value of a sub-integral is smaller than eps.
%
% Output: 
%   val: [1x1 real] Value of integral.
%
% References:
%   - Jianwei Zhu, Applications of Fourier Transform to Smile Modelling,
%   Springer, 2010.
%
    
    if b < a || b > c
        error('MultiDomainIntegration: Invalid integration points.');
    end
    
    dz = (c - b) / N;
    val = integral(f,a,b);
    if c > b
        for i=1:N
            valNew = integral(f,b + (i-1)*dz,b + i*dz);
            val = val + valNew;
            if abs(valNew) < eps
                return;
            end
        end
    end
    
end
