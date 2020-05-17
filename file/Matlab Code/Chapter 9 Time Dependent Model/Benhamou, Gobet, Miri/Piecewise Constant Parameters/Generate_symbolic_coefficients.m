% Benhamou, Gobet, and Miri (2010) symbolic expressions for the
% coefficients of the put price expansion in the Heston model

clc; clear;
syms theta kappa T u s t rho sigma Ti Tip1 v0 v0t v0s Tj Tjm1 
% Ti = T(t), Tip1 = T(t+1), Tj = T(j), Tjm1 = T(j-1)

%% Constant parameters
v0t = (exp(-kappa*t)*(v0-theta)+theta);
f = exp(kappa*t)*v0t*int(exp(-kappa*u),u,t,T);
a1T = rho*sigma*int(f,t,0,T);
a1T = simplify(a1T);
fprintf('Constant coefficient a1T\n');
fprintf('--------------------------------------------\n');
fprintf('a1T = %s\n',char(a1T));
fprintf('\n');

f = exp(kappa*t)*v0t*int(int(exp(-kappa*u),u,s,T),s,t,T);
a2T = (rho*sigma)^2*int(f,t,0,T);
a2T = simplify(a2T);
fprintf('Constant coefficient a2T\n');
fprintf('--------------------------------------------\n');
fprintf('a2T = %s\n',char(a2T));
fprintf('\n');

f = exp(2*kappa*t)*v0t*int(exp(-kappa*s)*int(exp(-kappa*u),u,s,T),s,t,T);
b0T = sigma^2*int(f,t,0,T);
b0T = simplify(b0T);
fprintf('Constant coefficient b0T\n');
fprintf('--------------------------------------------\n');
fprintf('b0T = %s\n',char(b0T));
fprintf('\n');

%% Piecewise constant parameters
% Coefficients for a1T(T+1)
% First element in the sum
f0 = rho*sigma * int(exp(kappa*t)*v0t*int(exp(-kappa*u),u,Ti,Tip1),t,   0,Tj);
f0 = simplify(f0);
% Remaining elements in the sum
f1 = rho*sigma * int(exp(kappa*t)*v0t*int(exp(-kappa*u),u,Ti,Tip1),t,Tjm1,Tj);
f1 = simplify(f1);
% Final term
f2 = rho*sigma * int(exp(kappa*t)*v0t*int(exp(-kappa*u),u,t, Tip1),t,Ti,Tip1);
f2 = simplify(f2);
fprintf('Piecewise constant coefficients for a1T(T+1)\n');
fprintf('--------------------------------------------\n');
fprintf('f0 = %s\n',char(f0));
fprintf('f1 = %s\n',char(f1));
fprintf('f2 = %s\n',char(f2));
fprintf('\n');


% Coefficients for a2T(T+1)
% First element in the first sum
f10 = (rho*sigma)^2 * int(exp(kappa*t)*v0t*int(int(exp(-kappa*u),u,s ,Tip1),s,Ti,Tip1),t,0,   Tj);
f10 = simplify(f10);
% Remaining elements in the first sum
f11 = (rho*sigma)^2 * int(exp(kappa*t)*v0t*int(int(exp(-kappa*u),u,s ,Tip1),s,Ti,Tip1),t,Tjm1,Tj);
f11 = simplify(f11);
% First element in the second sum
f20 = (rho*sigma)^2 * int(exp(kappa*t)*v0t*int(int(exp(-kappa*u),u,Ti,Tip1),s,t ,Ti  ),t,0,   Tj);
f20 = simplify(f20);
% Remaining elements in the second sum
f21 = (rho*sigma)^2 * int(exp(kappa*t)*v0t*int(int(exp(-kappa*u),u,Ti,Tip1),s,t ,Ti  ),t,Tjm1,Tj);
f21 = simplify(f21);
% Final term
f3 = (rho*sigma)^2 * int(exp(kappa*t)*v0t*int(int(exp(-kappa*u),u,s ,Tip1),s,t ,Tip1),t,Ti, Tip1);
f3 = simplify(f3);
fprintf('Piecewise constant coefficients for a2T(T+1)\n');
fprintf('--------------------------------------------\n');
fprintf('f10 = %s\n',char(f10));
fprintf('f11 = %s\n',char(f11));
fprintf('f20 = %s\n',char(f20));
fprintf('f21 = %s\n',char(f21));
fprintf('f3 = %s\n', char(f3));
fprintf('\n');
 
% b0T(T+1)
% First element in the first sum
f10 = sigma^2 * int(exp(2*kappa*t)*v0t*int(exp(-kappa*s)*int(exp(-kappa*u),u,s ,Tip1),s,Ti,Tip1),t,0,   Tj);
f10 = simplify(f10);
% Remaining elements in the first sum
f11 = sigma^2 * int(exp(2*kappa*t)*v0t*int(exp(-kappa*s)*int(exp(-kappa*u),u,s ,Tip1),s,Ti,Tip1),t,Tjm1,Tj);
f11 = simplify(f11);
% First element in the second sum
f20 = sigma^2 * int(exp(2*kappa*t)*v0t*int(exp(-kappa*s)*int(exp(-kappa*u),u,Ti,Tip1),s,t ,Ti ) ,t,0,   Tj);
f20 = simplify(f20);
% Remaining elements in the second sum
f21 = sigma^2 * int(exp(2*kappa*t)*v0t*int(exp(-kappa*s)*int(exp(-kappa*u),u,Ti,Tip1),s,t ,Ti ) ,t,Tjm1,Tj);
f21 = simplify(f21);
% Final term
f3 = sigma^2 * int(exp(2*kappa*t)*v0t*int(exp(-kappa*s)*int(exp(-kappa*u),u,s ,Tip1),s,t ,Tip1),t,Ti, Tip1);
f3 = simplify(f3);
fprintf('Piecewise constant coefficients for b0T(T+1)\n');
fprintf('--------------------------------------------\n');
fprintf('f10 = %s\n',char(f10));
fprintf('f11 = %s\n',char(f11));
fprintf('f20 = %s\n',char(f20));
fprintf('f21 = %s\n',char(f21));
fprintf('f3 = %s\n', char(f3));


