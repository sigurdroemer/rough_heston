function [y,Dalpha_y,t] = SolveVIE(f,alpha,T,N,M)
% Description: Solves the Volterra integral equation (VIE)
%              
%               y(t) = (1/gamma(alpha))*int_0^t(t-u)^(alpha - 1)f(u,y(u))du
%
%               for y(t) and that assuming 0 < alpha < 1 with f an input 
%               function. The equation is solved for all t on the interval 
%               0 <= t <= T. It is the user's responsibility to ensure that 
%               the equation actually has a (unique) solution and does not 
%               blow up on the interval [0,T].
%
%               The equation is solved on the equidistant grid
%               0 < T/N < 2*T/N < 3*T/N < ... < N*T/N = T
%
%               The notation is the same as used in (Diethelm, 2004).
% 
% Parameters:
%   f:      [function] Function of two variables from the VIE. We allow it 
%           to return a [Mx1] vector but then it should also accept [Mx1]
%           vectors for its second argument.
%   alpha:  [1x1 real] Number between 0 and 1.
%   T:      [1x1 real] Upper time point to solve VIE on.
%   N:      [1x1 integer] Number steps to discretize [0,T] into.
%   M:      [1x1 integer, optional] Output dimension of f (for speed).
%
% Output: 
%   y:        [Mx(N+1) real] Solution of VIE, 
%             i.e. [y(0),y(T/N),y(2*T/N),...,y(T)]
%   Dalpha_y: [Mx(N+1) real] Fractional (alpha) derivative of solution, 
%             i.e. [D^(alpha)y(0),D^(alpha)y(T/N),...,D^(alpha)y(T)]
%   t:        [1x(N+1) real] Discretization, i.e. [0,T/N,2*T/N,...,T]
%
% References:
%   - Kai Diethelm, Neville J. Ford and Alan D. Freed, Detailed error
%   analysis for a fractional ADAMs method, Numerical Algorithms 36:31-52,
%   2004.
%
    
    if ~exist('M','var') || isempty(M)
        M = size(f(0,0),1);
    end
    
    % Initialization:
    h = T / N;
    t = (0:h:T);
    [y,Dalpha_y] = deal(NaN(M,N+1));
    
    % Define coefficient functions:
    dummy1 = ( (h^alpha) / (alpha*(alpha + 1)) );
    a_0_kp1 = @(k)( dummy1 *( k.^(alpha+1)-(k - alpha).*((k+1).^(alpha))));
    a_j_kp1 = @(j,k)( dummy1*( (k-j+2).^(alpha+1) ...
                      + (k-j).^(alpha+1) - 2*(k-j+1).^(alpha+1) ) );
    a_kp1_kp1 = dummy1;
    b_j_kp1 = @(j,k) (  ((h^alpha)/alpha) * ( (k+1-j).^(alpha) ...
                                            - (k-j).^(alpha)  ) );
    
    % Run scheme:
    y(:,1) = 0;
    Dalpha_y(:,1) = f(0,0);
    for k=0:N-1
        js = (0:1:k); 
        
        % Compute predictor:
        yp = sum(b_j_kp1(js,k).*Dalpha_y(:,1:k+1),2)./gamma(alpha);
        
        % Compute solution:
        if k==0
            y(:,2) = (a_0_kp1(k)*Dalpha_y(:,1) ...
                      + a_kp1_kp1*f(t(k+2),yp))./gamma(alpha);
        else
            y(:,k+2) = (a_0_kp1(k)*Dalpha_y(:,1) ...
                        + sum(Dalpha_y(:,2:k+1).*a_j_kp1(js(2:end),k),2)...
                        + a_kp1_kp1*f(t(k+2),yp))./gamma(alpha);
        end
        
        % Compute fractional derivative:
        Dalpha_y(:,k+2) = f(h,y(:,k+2));
        
    end
    
    if any(any(isnan(y))) || any(any(isnan(Dalpha_y)))
        error('SolveVIE: NaNs produced!');
    end
    
end

