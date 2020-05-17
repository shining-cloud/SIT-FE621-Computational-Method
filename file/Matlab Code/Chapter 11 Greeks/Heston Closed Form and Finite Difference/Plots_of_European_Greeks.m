clc; clear;

%% Settings and parameters
S = 100;
K = 100;
r = 0.05;
q = 0.03;
T = 0.5 ;
kappa = 0;
lambda = 0;
sigma = .3;
v0 = 0.07;
theta = v0;
PutCall = 'C';
trap = 1;


%% Closed form Greeks
[x w] = GenerateGaussLaguerre(32);
T = [1/12:.01:.5];
S = [70:130];

GreekChoice = 'Theta';
rhopos =  0.9;
rhoneg = -0.9;
for s=1:length(S);
    for t=1:length(T);
        GreekPos(s,t) = HestonGreeks(PutCall,S(s),K,T(t),r,q,kappa,theta,sigma,lambda,v0,rhopos,trap,x,w,GreekChoice);
        GreekNeg(s,t) = HestonGreeks(PutCall,S(s),K,T(t),r,q,kappa,theta,sigma,lambda,v0,rhoneg,trap,x,w,GreekChoice);
    end
end
%% Plot the Greeks
mesh(GreekPos,'FaceColor','interp','EdgeColor','k','FaceLighting','phong')
hold on
mesh(GreekNeg,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
alpha(0.5)
hold off
axis tight
view([-300,150,70])

xlabel('Maturity')
ylabel('Strike Price')
set(gca,'XTick',      [ 1   10   20  30  40]);      % Maturity
set(gca,'XTickLabel', [.01  .12  .25  .33  .5 ])        
set(gca,'YTick',      [1   15  30   45   60]);      % Strike
set(gca,'YTickLabel', [70  85  100  115  130]);
set(gca,'YDir','reverse')

