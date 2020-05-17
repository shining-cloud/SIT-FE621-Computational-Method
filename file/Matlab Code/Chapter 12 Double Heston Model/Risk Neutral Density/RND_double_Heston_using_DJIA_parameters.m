% Risk Neutral estimation for the double Heston model

clc; clear;

% Market data
S = 129.14;
rf = 0.0010;
q  = 0.0068;
T = [37 72 135 226]./365;
NT = length(T);
trap = 1;

%% Double Heston parameters from DJIA estimation
Dparam = [3.058955819173225, 0.031716846087482, 1.984950961851184, 0.025844087663297,  0.064265856227069;
          1.846734553650094, 0.060513558098155, 0.714885834599290, 0.009236325814930, -0.975011498737402];

% Gauss Laguerre abscissas and weights
[x w] = GenerateGaussLaguerre(32);

% Initiate the cells
Calls  = cell(4,1);
Strike = cell(4,1);
RND = cell(4,1);
K   = cell(4,1);

% Strike increment and ranges
dK = .5;
K{1} = 60:dK:220;
K{2} = 40:dK:300;
K{3} = 18:dK:380;
K{4} = 8:dK:600;

%% Extract the RNDs, areas, and zeros
for t=1:NT
    NK = length(K{t});
    for k=1:NK
        Calls{t}(k) = DoubleHestonPriceGaussLaguerre('C',S,K{t}(k),T(t),rf,q,Dparam,x,w,trap);
    end
    [RND{t} ST{t}] = ExtractRND(K{t},Calls{t});
    Area(t) = trapz(RND{t})*dK;
    Zero(t) = length(find(RND{t}<0));
end

%% Output the results
fprintf('Risk Neutral Density estimation \n');
fprintf('Maturity       Area   NegValues \n');
fprintf('------------------------------- \n');
for t=1:NT
    fprintf(' %5.0f    %10.4f  %6.0f \n',T(t)*365,Area(t),Zero(t));
end
fprintf('------------------------------- \n');

%% Plot the RNDs
plot(ST{1},RND{1},'k-',ST{2},RND{2},'r-',ST{3},RND{3},'b-',ST{4},RND{4},'c')
legend([num2str(T(1)*365) '-day maturity'],...
       [num2str(T(2)*365) '-day maturity'],...
       [num2str(T(3)*365) '-day maturity'],...
       [num2str(T(4)*365) '-day maturity'])
xlabel('Terminal Stock Price');
ylabel('RND');
xlim([0 350]);

