%% Add paths:
currFolder = fileparts(matlab.desktop.editor.getActiveFilename);
idcs   = strfind(currFolder,'\');
projFolder = currFolder(1:idcs(end)-1);
addpath(genpath(projFolder));

%% Settings:
v_0=0.2^2;v_bar=v_0;alpha=0.6;lambda=2;xi=0.2;rho=-0.6;

%% How many terms should we use?
u_eval = (-50:1:0)';
n_max = [20,50,100,200,300,500];
num_stab = true;
for i=1:size(n_max,2)
    T_star = MomentExplosionTimeRoughHeston(alpha,lambda,xi,rho,...
                                            u_eval,n_max(i),num_stab);

    plot(u_eval,T_star,'o-','DisplayName',['n_{max} = ', num2str(n_max(i))]);hold on;
end
xlabel('Moment');ylabel('Explosion time');
title('Moment explosion times');
legend();

% Cut-off between case A and case B:
cut_off = lambda / (rho*xi)

% Conclusion: n_max = 200 is sufficient and picking n_max too high results 
% in numerical issues for extreme moments (see n_max = 500 case).


%% Numerical problems with gamma ratio
% In 'MomentExplosionTimeRoughHeston' we need to compute a ratio of two
% gamma functions. For large values this ratio causes numerical problems:
v_fun = @(n)(gamma(alpha*n + 1) ./ gamma(alpha*n - alpha + 1));
n_test = (10:500)';
figure;
plot(n_test,v_fun(n_test));
xlabel('n');ylabel('v(n)');
v_fun(n_test)

% Solution: Use asymptotic limit of ratio for large n - easily derived
% using Gautschi's inequality:

figure;
v_asympt = @(n)((n*alpha).^(alpha));
plot(n_test,v_fun(n_test));hold on;
plot(n_test,v_asympt(n_test));
legend('Direct formula','Asymptotic formula');
xlabel('n');ylabel('v(n)');
title('Gamma ratio');


%% Slow convergence around cut-off point (between case A and B)
% We observe slow convergence around the point that separates cases A and B
% (Gerhold et al., 2019). This could cause problems when trying to find
% the lower critical moment. 
%
% References: 
%   - Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment
%   explosions in the rough Heston model, Decisions in Economics and
%   Finance (2019) 42:575-608.


nmax = [50,100,200,500];
u_eval = (cut_off-0.1:0.001:cut_off+0.1)';
figure;
for i=1:size(nmax,2)
    T_star = MomentExplosionTimeRoughHeston(alpha,lambda,xi,rho,...
                                            u_eval,nmax(i),false);
    plot(u_eval,T_star,'-o','DisplayName',['nmax = ', num2str(nmax(i))]);hold on;
end
xlabel('Moment');ylabel('Explosion time');
title('Moment explosion times around cut-off point');
legend();

% Solution: Linearly interpolate for moments close to lambda/(xi*rho):
nmax = 200;
figure;
u_eval = (cut_off-0.1:0.0001:cut_off+0.1)';
T_star = MomentExplosionTimeRoughHeston(alpha,lambda,xi,rho,...
                                            u_eval,nmax,true);
plot(u_eval,T_star,'-');hold on;
T_star = MomentExplosionTimeRoughHeston(alpha,lambda,xi,rho,...
                                            u_eval,nmax,false);
plot(u_eval,T_star,'-o');
xlabel('Moment');ylabel('Explosion time');
title('Explosion times with and without linear interpolation');
legend('numstab = true','numstab=false');


