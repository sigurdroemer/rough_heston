%% Clear and add paths:
clear;
currFolder = fileparts(matlab.desktop.editor.getActiveFilename);
idcs   = strfind(currFolder,'\');
projFolder = currFolder(1:idcs(end)-1);
addpath(genpath(projFolder));


%% Finding an optimal alpha
% Here we illustrate the effect of using an optimal alpha (= beta - 1)
% along the lines of (Lord & Kahl, 2006).
%
% References:
%   - Roger lord and Christian Kahl, Optimal Fourier inversion in 
%     semi-analytical option pricing, available at SRRN: 
%     https://papers.ssrn.com/sol3/papers.cfm?abstract_id=921336, 2006.
%

% Set-up:
v_0=0.2^2;v_bar=v_0;alpha_model=0.6;lambda=2;xi=0.2;rho=-0.6;
N = 250;T=0.1;nmax = 100;
phi = @(u)(MomentGeneratingFunctionRoughHeston(v_0,alpha_model,lambda,...
                                                   v_bar,xi,rho,T,1i*u,N));
T_star = @(u)(MomentExplosionTimeRoughHeston(alpha_model,lambda,xi,rho,u,nmax,true));
l_critical = GetLowerCriticalMomentRoughHeston(lambda,rho,xi,T_star,T);
l_critical

% Define integrand:
K = 0.75;
k = log(K);
psi = @(alpha,v) (  phi( v - (alpha + 1)*1i  ) ./ ...
    ( alpha^2 + alpha - v.^2 + 1i*(2*alpha + 1).*v)  );
integrand = @(alpha,v)(real(psi(alpha,v).*exp(-1i*v*k')).*(exp(-alpha*k')./pi));


%% Plot initial value of integrand across different alpha's:
% 'Alpha' here refers to the notation in (Lord & Kahl, 2006) and not the one 
% from the rough Heston model.

delta = 0.01;
alpha_min = l_critical - 1;
alpha_max = -1;delta_min = 1;delta_max = 10;
alpha_test = (alpha_min+delta:1:alpha_max-delta_max)';
val0 = NaN(size(alpha_test));
for i=1:size(alpha_test,1)
    val0(i) = integrand(alpha_test(i),0);
end
figure;
plot(alpha_test,val0,'o-');xlabel('alpha');ylabel('Integrand at 0');
title('Optimal alpha');


%% Solve for optimal one:
alpha_opt = GetOptimalFourierAlpha(phi,log(K),l_critical)

%% Illustrate effect of alpha on oscillations:
v_test = (0:0.1:40)';
alpha_test = [l_critical+0.1,alpha_opt,-10,-3];
figure;
for i=1:size(alpha_test,2)
    plot(v_test,integrand(alpha_test(i),v_test)./integrand(alpha_test(i),0),'-',...
        'DisplayName',['alpha = ', num2str(alpha_test(i))]);
    hold on;
end
legend();
xlabel('x');
ylabel('integrand(x)/integrand(0)');
title('Normalised integrand for different alpha''s');











