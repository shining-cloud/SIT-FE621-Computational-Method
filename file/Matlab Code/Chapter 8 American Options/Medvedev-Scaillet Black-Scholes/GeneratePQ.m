clc; clear;

% Medvedev-Scaillet coefficients for p0n and q0n
% Finds coefficients up to 7 terms
% Can be easily modified to find more terms

p0n = inline('(n-2*i)*(n-2*i-1)/(2*i+2)*pi','n','i','pi');
q0n = inline('(qi*(n-1-2*i)*(n-2-2*i) + 2*pi*(n-2*i-2))/(2*n-2*i-2)','n','i','qi','pi');

Q = {'q01','q02','q03','q04','q05','q06','q07'};
P = {'p01','p02','p03','p04','p05','p06','p07'};

N = 7;

%% Coefficients for p0n
p0 = zeros(N,N+1);
for n=1:N
    p0(n,N-n+1) = 1;
    i = 0;
    for j=N-n+3:2:N+1
        p0(n,j) = p0n(n,i,p0(n,j-2));
        i = i+1;
    end
end

fprintf('Coefficients for p0n(theta)\n')
fprintf('p0n  theta^7     theta^6   theta^5    theta^4    theta^3    theta^2    theta^1    theta^0\n')
fprintf('-----------------------------------------------------------------------------------------\n')
for n=1:N
    fprintf('%s %5.0f %10.0f %10.0f %10.0f %10.0f %10.0f %10.0f %10.0f\n', cell2mat(P(n)), p0(n,:))
end
fprintf('-----------------------------------------------------------------------------------------\n')
fprintf('\n');

%% Coefficients for q0n
q0 = zeros(N,N+1);
for n=1:N
    q0(n,N-n+2) = 1;
    i = 0;
     for j=N-n+4:2:N+1
         q0(n,j) = q0n(n,i,q0(n,j-2),p0(n,j-1));
         i = i+1;
     end
end

fprintf('Coefficients for q0n(theta)\n')
fprintf('q0n  theta^7     theta^6   theta^5    theta^4    theta^3    theta^2    theta^1    theta^0\n')
fprintf('-----------------------------------------------------------------------------------------\n')
for n=1:N
    fprintf('%s %5.0f %10.0f %10.0f %10.0f %10.0f %10.0f %10.0f %10.0f\n', cell2mat(Q(n)), q0(n,:))
end
fprintf('-----------------------------------------------------------------------------------------\n')
clear p0n q0n


%% First set of polynomials
n = 1;
fprintf('Set 1 polynomials\n')

% Polynomials and their derivatives
[p0n q0n p1n q1n p1n_1 q1n_1 p1n_2 q1n_2 p0n_1 q0n_1 p0n_2 q0n_2 ...
 dp1n d2p1n dp0n_1 dp1n_1 dq1n d2q1n dq0n_1 dq1n_1] = PQ(n,N,p0,q0);
p01 = p0n;
p11 = p1n;
q01 = q0n;
q11 = q1n;

syms theta
fprintf('-----------------\n')
fprintf('p01 = %s \n',char(p01))
fprintf('p11 = %s \n',char(p11))
fprintf('q01 = %s \n',char(q01))
fprintf('q11 = %s \n',char(q11))
fprintf('\n')


%% Second set of polynomials 
n = 2;
fprintf('Set 2 polynomials\n')

% Polynomials and their derivatives
[p0n q0n p1n q1n p1n_1 q1n_1 p1n_2 q1n_2 p0n_1 q0n_1 p0n_2 q0n_2 ...
 dp1n d2p1n dp0n_1 dp1n_1 dq1n d2q1n dq0n_1 dq1n_1] = PQ(n,N,p0,q0);
p02 = p0n;
q02 = q0n;

% Solve the first PDE for p12
syms sigma mu C1 C0 r pi21 sigma_
sigma_ = (sigma^2-2*mu)/sigma;
PDE1 = d2p1n + theta*dp1n - n*p1n + sigma_*C1*dp0n_1 + sigma_*dp1n_1 - 2*r*C0*p0n_2 - 2*r*p1n_2;
pi21 = solve(PDE1,pi21);
p12 = pi21;
fprintf('p12 is of the form p12 = %s \n', char(p1n))

% Solve the second PDE for q12
q12 = 0;

fprintf('-----------------\n')
fprintf('p02 = %s \n', char(p02))
fprintf('p12 = %s \n', char(p12))
fprintf('q02 = %s \n', char(q02))
fprintf('q12 = %1.0f \n',q12)
fprintf('\n')
fprintf('\n')

%% Third set of polynomials
n = 3;
fprintf('Set 3 polynomials\n')

% Polynomials and their derivatives
[p0n q0n p1n q1n p1n_1 q1n_1 p1n_2 q1n_2 p0n_1 q0n_1 p0n_2 q0n_2 ...
 dp1n d2p1n dp0n_1 dp1n_1 dq1n d2q1n dq0n_1 dq1n_1] = PQ(n,N,p0,q0);
p03 = p0n;
q03 = q0n;

% Solve the first PDE for p13
syms pi31 C2 C1
PDE1 = d2p1n + theta*dp1n - n*p1n + sigma_*C2*dp0n_1 + sigma_*dp1n_1 - 2*r*C1*p0n_2 - 2*r*p1n_2;
pi31 = solve(PDE1,pi31);
p13 = pi31*theta;
fprintf('p13 is of the form p13 = %s \n', char(p1n))

% Solve the second PDE for q13
syms qi31 
PDE2 = -(n+1)*q1n - theta*dq1n + d2q1n + 2*dp1n + sigma_*C2*p0n_1 ...
     + sigma_*C2*dq0n_1 - sigma_*C2*theta*q0n_1 + sigma_*p1n_1 + sigma_*dq1n_1 - sigma_*theta*q1n_1 ...
     - 2*r*(C1*q0n_2 + q1n_2);
qi31 = solve(PDE2,qi31);
q13 = simplify(subs(qi31));
fprintf('q13 is of the form q13 = %s \n', char(q1n))

fprintf('-----------------\n')
fprintf('p03 = %s \n',char(p03))
fprintf('p13 = %s \n',char(p13))
fprintf('q03 = %s \n',char(q03))
fprintf('q13 = %s \n',char(q13))
fprintf('\n')
fprintf('\n')
 
%% Fourth set of polynomials
n = 4;
fprintf('Set 4 polynomials\n')

% Polynomials and their derivatives
[p0n q0n p1n q1n p1n_1 q1n_1 p1n_2 q1n_2 p0n_1 q0n_1 p0n_2 q0n_2 ...
 dp1n d2p1n dp0n_1 dp1n_1 dq1n d2q1n dq0n_1 dq1n_1] = PQ(n,N,p0,q0);
p04 = p0n;
q04 = q0n;

% Solve the first PDE for p14
syms pi41 C3 C2
PDE1 = d2p1n + theta*dp1n - n*p1n + sigma_*C3*dp0n_1 + sigma_*dp1n_1 - 2*r*C2*p0n_2 - 2*r*p1n_2;
PDE1 = subs(PDE1);
collect(PDE1,theta);
pi41 = solve('(-2*pi41-2*r*C2+3*(sigma^2-2*mu)/sigma*C3)*theta^2','pi41');
pi42 = solve('+2*pi41-(sigma^2-2*mu)/sigma^2*(-sigma^2*C2+2*C2*mu+r*C1*sigma)-4*pi42+3*(sigma^2-2*mu)/sigma*C3+r*(-sigma^2+2*mu)*C1/sigma-2*r*C2','pi42');
pi42 = simplify(subs(pi42));
p14 = pi41*theta^2 + pi42;
fprintf('p14 is of the form p14 = %s \n', char(p1n))

% Solve the second PDE for q14
syms qi41
PDE2 = -(n+1)*q1n - theta*dq1n + d2q1n + 2*dp1n + sigma_*C3*p0n_1 ...
     + sigma_*C3*dq0n_1 - sigma_*C3*theta*q0n_1 + sigma_*p1n_1 + sigma_*dq1n_1 - sigma_*theta*q1n_1 ...
     - 2*r*(C2*q0n_2 + q1n_2);
qi41 = solve(PDE2,qi41);
qi41 = subs(qi41);
qi41 = subs(qi41);
qi41 = simplify(subs(qi41));
q14 = qi41*theta;
fprintf('q14 is of the form q14 = %s \n', char(q1n))

fprintf('-----------------\n')
fprintf('p04 = %s \n',char(p04))
fprintf('p14 = %s \n',char(p14))
fprintf('q04 = %s \n',char(q04))
fprintf('q14 = %s \n',char(q14))
fprintf('\n')
fprintf('\n')

 
%% Fifth set of polynomials
n = 5;
fprintf('Set 5 polynomials\n')

% Polynomials and their derivatives
[p0n q0n p1n q1n p1n_1 q1n_1 p1n_2 q1n_2 p0n_1 q0n_1 p0n_2 q0n_2 ...
 dp1n d2p1n dp0n_1 dp1n_1 dq1n d2q1n dq0n_1 dq1n_1] = PQ(n,N,p0,q0);
p05 = p0n;
q05 = q0n;

% Solve the first PDE for P15
syms pi51 C4 C3
PDE1 = d2p1n + theta*dp1n - n*p1n + sigma_*C4*dp0n_1 + sigma_*dp1n_1 - 2*r*C3*p0n_2 - 2*r*p1n_2;
PDE1 = subs(PDE1);
collect(PDE1,theta);
pi51 = solve('-2*pi51-2*r*C3+4*(sigma^2-2*mu)/sigma*C4','pi51');
pi52 = solve('6*pi51-4*pi52+2*r*(-sigma^2*C2+2*C2*mu+r*C1*sigma)/sigma-1/2*(2*sigma^2-4*mu)/sigma^2*(2*r*C2*sigma-3*C3*sigma^2+6*C3*mu)+12*(sigma^2-2*mu)/sigma*C4-6*r*C3','pi52');
pi52 = subs(pi52);
pi52 = simplify(pi52);
p15 = pi51*theta^3 + pi52*theta;
fprintf('P15 is of the form P15 = %s \n', char(p1n))

% Solve the second PDE for Q15
syms qi51
PDE2 = -(n+1)*q1n - theta*dq1n + d2q1n + 2*dp1n + sigma_*C4*p0n_1 ...
      + sigma_*C4*dq0n_1 - sigma_*C4*theta*q0n_1 + sigma_*p1n_1 + sigma_*dq1n_1 - sigma_*theta*q1n_1 ...
      - 2*r*(C3*q0n_2 + q1n_2);
PDE2 = subs(PDE2);
collect(PDE2,theta);
qi51 = solve('1/48*(sigma^2-2*mu)/sigma^4*(-72*C3*sigma^4+144*C3*mu*sigma^2+48*r*C2*sigma^3+sigma^6*C1-6*sigma^4*C1*mu+12*sigma^2*C1*mu^2-8*C1*mu^3)+4*(sigma^2-2*mu)/sigma*C4-8*qi51-1/2*(sigma^2-2*mu)/sigma^2*(2*r*C2*sigma-3*C3*sigma^2+6*C3*mu)-6*(r*C3*sigma-2*C4*sigma^2+4*C4*mu)/sigma-2*r*C3','qi51');
qi52 = solve('1/2*(-12*r*C3*sigma^2+24*C4*sigma^3-48*C4*sigma*mu-4*r*C2*sigma^3+8*r*sigma*C2*mu+2*r^2*sigma^2*C1+3*C3*sigma^4-12*C3*mu*sigma^2+12*C3*mu^2)/sigma^2-6*qi52+2*qi51-2*r*(2*C3-1/4*(-2*pi31*sigma-2*sigma^2*C2+4*C2*mu-pi21*sigma^2+2*pi21*mu+2*r*C1*sigma)/sigma)+8*(sigma^2-2*mu)/sigma*C4-1/48*(sigma^2-2*mu)/sigma^4*(-72*C3*sigma^4+144*C3*mu*sigma^2+48*r*C2*sigma^3+sigma^6*C1-6*sigma^4*C1*mu+12*sigma^2*C1*mu^2-8*C1*mu^3)+1/4*(sigma^2-2*mu)/sigma^3*(-4*r*C2*sigma^2+6*C3*sigma^3-12*C3*sigma*mu+sigma^4*C2-4*sigma^2*C2*mu-2*sigma^3*r*C1+4*C2*mu^2+4*mu*r*C1*sigma)','qi52');
qi52 = subs(qi52);
qi51 = simplify(qi51);
qi52 = simplify(qi52);
q15 = qi51*theta^2 + qi52;
fprintf('q15 is of the form q15 = %s \n', char(q1n))

fprintf('-----------------\n')
fprintf('p05 = %s \n',char(p05))
fprintf('p15 = %s \n',char(p15))
fprintf('q05 = %s \n',char(q05))
fprintf('q15 = %s \n',char(q15))
fprintf('\n')
fprintf('\n')


%% Sixth set of polynomials
n = 6;
fprintf('Set 6 polynomials\n')

% Polynomials and their derivatives
[p0n q0n p1n q1n p1n_1 q1n_1 p1n_2 q1n_2 p0n_1 q0n_1 p0n_2 q0n_2 ...
 dp1n d2p1n dp0n_1 dp1n_1 dq1n d2q1n dq0n_1 dq1n_1] = PQ(n,N,p0,q0);
p06 = p0n;
q06 = q0n;

% Solve the first PDE for P16
syms pi61 pi62 C5 C4
PDE1 = d2p1n + theta*dp1n - n*p1n + sigma_*C5*dp0n_1 + sigma_*dp1n_1 - 2*r*C4*p0n_2 - 2*r*p1n_2;
PDE1 = subs(PDE1);
PDE1 = subs(PDE1);
AA = collect(PDE1,theta);
pi61 = solve('(-2*pi61+5*(sigma^2-2*mu)/sigma*C5-2*r*C4)','pi61');
pi62 = solve('(12*pi61-4*pi62-2*r*(-r*C2*sigma+3/2*sigma^2*C3-3*C3*mu)/sigma+30*(sigma^2-2*mu)/sigma*C5-12*r*C4+(sigma^2-2*mu)/sigma^2*(-3*r*C3*sigma+6*C4*sigma^2-12*C4*mu))','pi62');
pi63 = solve('-6*r*C4+2*pi62-2*r*(-r*C2*sigma^2+3/2*C3*sigma^3-3*C3*sigma*mu+1/4*sigma^4*C2-sigma^2*C2*mu-1/2*sigma^3*r*C1+C2*mu^2+mu*r*C1*sigma)/sigma^2+15*(sigma^2-2*mu)/sigma*C5-6*pi63+(sigma^2-2*mu)/sigma^3*(-3*r*C3*sigma^2+6*C4*sigma^3-12*C4*sigma*mu-r*C2*sigma^3+2*r*sigma*C2*mu+1/2*r^2*sigma^2*C1+3/4*sigma^4*C3-3*C3*mu*sigma^2+3*C3*mu^2)','pi63');
p16 = pi61*theta^4 + pi62*theta^2 + pi63;
fprintf('p16 is of the form p16 = %s \n', char(p1n))

% Solve the second PDE for q15
syms qi61
PDE2 = -(n+1)*q1n - theta*dq1n + d2q1n + 2*dp1n + sigma_*C5*p0n_1 ...
      + sigma_*C5*dq0n_1 - sigma_*C5*theta*q0n_1 + sigma_*p1n_1 + sigma_*dq1n_1 - sigma_*theta*q1n_1 ...
      - 2*r*(C4*q0n_2 + q1n_2);
PDE2 = subs(PDE2);
PDE2 = subs(PDE2);
collect(PDE2,theta);
qi61 = solve('(-10*qi61-(-20*C5*sigma^2+40*C5*mu+8*r*C4*sigma)/sigma+5*(sigma^2-2*mu)/sigma*C5-(sigma^2-2*mu)/sigma^5*(1/384*sigma^8*C1-1/48*sigma^6*C1*mu+1/16*sigma^4*C1*mu^2-1/12*sigma^2*C1*mu^3+1/24*C1*mu^4+2*C4*sigma^5-4*C4*sigma^3*mu-r*C3*sigma^4)+(sigma^2-2*mu)/sigma^2*(-r*C3*sigma+2*C4*sigma^2-4*C4*mu)-2*r*C4)','qi61');
qi62 = solve('(6*qi61-8*qi62+25*(sigma^2-2*mu)/sigma*C5-(sigma^2-2*mu)/sigma^5*(-3*C3*sigma^4*mu-1/3*mu^2*r*C1*sigma^2+4/3*mu*r*C2*sigma^3-1/384*sigma^8*C1-1/24*C1*mu^4+1/6*sigma^4*r^2*C1+1/2*sigma^3*C2*mu^2+3/4*sigma^6*C3-1/3*sigma*C2*mu^3-2/3*sigma^5*r*C2+2/3*r*C3*sigma^3*C2*mu+4*C4*sigma^5+1/3*sigma^4*mu*r*C1+1/6*r*C3*sigma^4*C1*mu-8*C4*sigma^3*mu-1/3*r*C3*sigma^5*C2-2*r*C3*sigma^4+3*C3*mu^2*sigma^2-1/6*r*C3*sigma^2*C1*mu^2-1/24*r*C3*sigma^6*C1-1/12*sigma^6*r*C1+1/24*sigma^7*C2+1/48*sigma^6*C1*mu-1/16*sigma^4*C1*mu^2-1/4*sigma^5*C2*mu+1/12*sigma^2*C1*mu^3+1/3*r^2*C3*sigma^4*C1)+(4*(15/2*C5*sigma^2-15*C5*mu-3*r*C4*sigma)*sigma-6*r*C3*sigma^3+6*sigma^4*C4-24*sigma^2*C4*mu+12*mu*r*C3*sigma+24*C4*mu^2+2*r^2*sigma^2*C2+30*C5*sigma^3-60*C5*sigma*mu-12*r*C4*sigma^2)/sigma^2+(sigma^2-2*mu)/sigma^3*(-3*r*C3*sigma^2+6*C4*sigma^3-12*C4*sigma*mu-r*C2*sigma^3+2*r*sigma*C2*mu+1/2*r^2*sigma^2*C1+3/4*C3*sigma^4-3*C3*mu*sigma^2+3*C3*mu^2)-10*r*C4+1/24*r*C4*(48*r*C2*sigma^3-72*C3*sigma^4+144*C3*mu*sigma^2+sigma^6*C1-6*sigma^4*C1*mu+12*sigma^2*C1*mu^2-8*C1*mu^3)/sigma^3-(1/192*sigma^2-1/96*mu)/sigma^5*(-sigma^8*C1+8*sigma^6*C1*mu-24*sigma^4*C1*mu^2+32*sigma^2*C1*mu^3-16*C1*mu^4-768*C4*sigma^5+1536*C4*sigma^3*mu+384*r*C3*sigma^4))','qi62');
q16 = qi61*theta^3 + qi62*theta;
fprintf('q16 is of the form q16 = %s \n', char(q1n))

fprintf('-----------------\n')
fprintf('p06 = %s \n',char(p06))
fprintf('p16 = %s \n',char(p16))
fprintf('q06 = %s \n',char(q06))
fprintf('q16 = %s \n',char(q16))
fprintf('\n')
fprintf('\n')


%% Find the "C" coefficients
fprintf('"C" Coefficients\n')
fprintf('----------------\n')
syms pdf cdf y K
Exp1 = sigma*y*K;
C1 = (Exp1 - p11*cdf - q11*pdf)/(p01*cdf + q01*pdf);
C1 = subs(C1,theta,y);
fprintf('C1 = %s \n', char(C1))

Exp2 = -sigma^2*y^2*K/2;
C2 = (Exp2 - p12*cdf - q12*pdf)/(p02*cdf + q02*pdf);
C2 = subs(C2,theta,y);
C2 = simplify(C2);
fprintf('C2 = %s \n', char(C2))

Exp3 = sigma^3*y^3*K/6;
C3 = (Exp3 - p13*cdf - q13*pdf)/(p03*cdf + q03*pdf);
C3 = subs(C3,theta,y);
C3 = simplify(C3);
fprintf('C3 = %s \n', char(C3))

Exp4 = -sigma^4*y^4*K/24;
C4 = (Exp4 - p14*cdf - q14*pdf)/(p04*cdf + q04*pdf);
C4 = subs(C4,theta,y);
C4 = simplify(C4);
fprintf('C4 = %s \n', char(C4))

Exp5 = sigma^5*y^5*K/120;
C5 = (Exp5 - p15*cdf - q15*pdf)/(p05*cdf + q05*pdf);
C5 = subs(C5,theta,y);
C5 = simplify(C5);
fprintf('C5 = %s \n', char(C5))
fprintf('\n')
