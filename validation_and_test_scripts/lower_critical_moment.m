%% Clear and add paths:
clear;
currFolder = fileparts(matlab.desktop.editor.getActiveFilename);
idcs   = strfind(currFolder,'\');
projFolder = currFolder(1:idcs(end)-1);
addpath(genpath(projFolder));


%% Computing the lower critical moment
% Here we illustrate the computation of the lower critical moment for
% different time points. Note that the estimated moments may be biased high
% c.f. the moment explosion times being biased low in case B of (Gerhold et
% al., 2019).
%
% References: 
%   - Stefan Gerhold, Christoph Gerstenecker, Arpad Pinter, Moment
%   explosions in the rough Heston model, Decisions in Economics and
%   Finance (2019) 42:575-608.

alpha=0.6;lambda=2;xi=0.2;rho=-0.8;nmax = 100;

T_star = @(u)(MomentExplosionTimeRoughHeston(alpha,lambda,xi,rho,u,nmax,true));
T_eval = [(0.01:0.01:0.1)';(0.15:0.05:0.5)';(1:0.5:3)'];
l_critical = GetLowerCriticalMomentRoughHeston(lambda,rho,xi,T_star,T_eval);

plot(T_eval,l_critical,'o-');
xlabel('T');ylabel('Moment');title('Lower critical moments');
xlim([-0.1,3]);







