%% Clear and add paths:
clear;
currFolder = fileparts(matlab.desktop.editor.getActiveFilename);
idcs   = strfind(currFolder,'\');
projFolder = currFolder(1:idcs(end)-1);
addpath(genpath(projFolder));

     
%% Define model
s_0=100;v_0=0.15.^2;v_bar=0.15.^2;alpha=0.6;lambda=2;xi=.4;rho=-0.6;


%% Verify against Mathworks code in classical Heston case:
settle = datenum('20190101','yyyymmdd');
maturity = datenum('20200101','yyyymmdd');
ttm = yearfrac(settle, maturity, 0);
K = (50:5:150)';
r = 0.05;q = 0.02;
call = true; T = 1;
price_true = optByHestonNI(r,s_0,settle,maturity,'call',K,v_0,...
                           v_bar,lambda,xi,rho,'DividendYield',q);
iv_true = blsimpv(s_0,K,r,ttm,price_true,'Yield',q);

[price, iv] = NumericalIntegrationRoughHeston(s_0,v_0,1,lambda,...
                                        v_bar,xi,rho,call,K,T,...
                                        'disp_iter',true,'r',r,'q',q);

figure;
yyaxis left
plot(log(K/s_0),iv - iv_true,'o-');
ylabel('Difference in implied volatility');
yyaxis right
plot(log(K/s_0),price - price_true,'o-');
ylabel('Difference in (call option) prices');
xlabel('Log-moneyness');title('Mathworks vs. rHeston code (alpha = 1)');

[price,price_true]
[iv,iv_true]

% Conclusion: Error is around 10^(-4) in implied volatility with default
% settings.


%% Check effect of varying 'eps':
% I.e. the tolerance level used for the multi-domain part of computing the
% pricing integral.

exponent = (0:10);
eps = [10.^(-exponent)];
K = [90,100,110]';T = 0.1;call = true;
[price,iv] = deal(NaN(size(K,1),size(eps,1)));
for i=1:size(eps,2)
    num2str(i);
    try
        [price(:,i), iv(:,i)] = NumericalIntegrationRoughHeston(s_0,v_0,1,lambda,...
                                                v_bar,xi,rho,call,K,T,...
                                                'eps',eps(i));
    catch
        % Large 'eps' values may cause negative prices and thus an error
    end
end

num2cell([eps;price])'
num2cell([eps;iv])'


%% Check effect of varying 'N':
% I.e. the number of steps when solving the Volterra integral equation.
% Remark: This may take a few minutes to run.

N = [100,252,500,1000,2000];
K = [90,100,110]';T = 0.1;call = true;
[price,iv] = deal(NaN(size(K,1),size(N,1)));
for i=1:size(N,2)
    num2str(i)
    [price(:,i), iv(:,i)] = NumericalIntegrationRoughHeston(s_0,v_0,1,...
                                            lambda,v_bar,xi,rho,call,...
                                            K,T,'N',N(i));
end

num2cell([N;price])'
num2cell([N;iv])'

