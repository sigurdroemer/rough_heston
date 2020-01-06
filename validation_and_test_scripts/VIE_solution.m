%% Add paths:
currFolder = fileparts(matlab.desktop.editor.getActiveFilename);
idcs   = strfind(currFolder,'\');
projFolder = currFolder(1:idcs(end)-1);
addpath(genpath(projFolder));

%% Solving the Volterra integral equation (VIE):
% Description: We show convergence of the solution of the VIE in the number
% of steps used.

% Settings:
u = -3;rho = -0.6;alpha = 0.6;lambda=2;xi=0.2;
c1 = 0.5*u.*(u-1);
c2 =  rho*xi*u - lambda;
c3 = 0.5*xi.^2;
R = @(w) (c1 + c2.*w + c3.*w.^2);
f = @(s,y) (R(y));
T = 1;
N = [100, 250, 500, 1000, 2000];

[y,Dalpha_y,t] = deal(cell(size(N,2),1));
for i=1:size(N,2)
    [y{i},Dalpha_y{i},t{i}] = SolveVIE(f,alpha,T,N(i));
end

figure;
for i=1:size(y,1)
    plot(t{i},y{i},'DisplayName',['N = ', num2str(N(i))]);hold on;
end
legend();
xlabel('t');
ylabel('y(t)');
title('Solution to VIE');

figure;
for i=1:size(y,1)
    plot(t{i},Dalpha_y{i},'DisplayName',['N = ', num2str(N(i))]);hold on;
end
legend();
xlabel('t');
ylabel('D^{alpha}_t y(t)','interpreter','tex');
title('Fractional derivative');
