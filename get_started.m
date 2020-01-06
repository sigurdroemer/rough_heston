%% Add code to matlab path:
projFolder = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(projFolder);

%% Example 1 (A single smiles):
v_0 = 0.15.^2;v_bar=0.15.^2;alpha = 0.6;lambda = 2;xi = .4;rho = -0.6;
s0 = 100;K = (50:5:150)';T = 1;call = true;
[price, iv] = NumericalIntegrationRoughHeston(s0,v_0,alpha,lambda,...
                                        v_bar,xi,rho,call,K,T,...
                                        'disp_iter',true);
                                             
[price,iv]

figure;
yyaxis left
plot(log(K/s0),iv,'o-','DisplayName','Implied volatility');hold on;
ylabel('Implied volatility');
yyaxis right
plot(log(K/s0),price,'o-','DisplayName','Call price');
ylabel('Price');
xlabel('Log-moneyness');title('Rough Heston prices');

% Save figure:
% saveas(gcf,'smile1.jpg');

%% Example 2 (Multiple smiles - explosion of at-the-money skew):
K = (95:1:103)';T = [0.01;0.02;0.1;0.5;1];
[price, iv] = NumericalIntegrationRoughHeston(s0,v_0,alpha,lambda,...
                                        v_bar,xi,rho,call,K,T,'N',252,...
                                        'disp_iter',true);

figure;
c = hsv(size(T,1));
for i=1:size(T,1)
    plot(log(K/s0),iv(:,i),'o-','Color',c(i,:),'linewidth',1.2);hold on;
end
xlabel('Log-moneyness');
ylabel('Implied volatility');
title('Rough Heston smiles');
hleg = legend(cellstr(num2str(T, 'T = %-2.2f')),'location','best');
title(hleg,'Expiry')

% Save figure:
% saveas(gcf,'smile2.jpg');


