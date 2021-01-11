%% Add code to matlab path:
projFolder = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(projFolder));

%% Example 1 (A single smile):
v_0 = 0.15.^2;v_bar=0.15.^2;alpha = 0.6;lambda = 2;xi = .4;rho = -0.6;
s0 = 100;K = (50:5:150)';T = 1;call = true;
[priceC, ivC] = NumericalIntegrationRoughHeston(s0,v_0,alpha,lambda,...
                                        v_bar,xi,rho,call,K,T,...
                                        'disp_iter',true);
                                             
[priceC,ivC]

call = false;
[priceP, ivP] = NumericalIntegrationRoughHeston(s0,v_0,alpha,lambda,...
                                        v_bar,xi,rho,call,K,T,...
                                        'disp_iter',true);
                                    
                                    
% C - P = S - K
priceC - priceP - (s0 - K)
[priceC,priceP]

% Looks ok.
ivC - ivP
[ivC,ivP]