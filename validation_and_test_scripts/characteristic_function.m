%% Clear and add paths:
clear all;
currFolder = fileparts(matlab.desktop.editor.getActiveFilename);
idcs   = strfind(currFolder,'\');
projFolder = currFolder(1:idcs(end)-1);
addpath(genpath(projFolder));

%% Computing the characteristic function
v_0=0.2^2;v_bar=v_0;alpha=0.6;lambda=2;xi=0.2;rho=-0.6;T=0.1;

% Check accuracy against the number of steps:
N = [100, 250, 500, 1000, 2000];
phi = @(u,N)(MomentGeneratingFunctionRoughHeston(v_0,alpha,lambda,...
                                                   v_bar,xi,rho,T,1i*u,N));

u_eval = (-10:1:10)';
figure;
for i=1:size(N,2)
    plot(u_eval,imag(phi(u_eval,N(i))),'DisplayName',['N = ', num2str(N(i))]);hold on;
end
legend();
xlabel('u');
ylabel('imaginary part');
title('Characteristic function (imaginary part)');

figure;
for i=1:size(N,2)
    plot(u_eval,real(phi(u_eval,N(i))),'DisplayName',['N = ', num2str(N(i))]);hold on;
end
legend();
xlabel('u');
ylabel('Real part');
title('Characteristic function (real part)');









