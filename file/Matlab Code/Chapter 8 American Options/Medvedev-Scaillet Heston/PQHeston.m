function [p0n       q0n ...        p0(n)                 q0(n)
          p1n       q1n ...        p1(n)                 q1(n)
          p1n_1     q1n_1 ...      p1(n-1)               q1(n-1)
          p1n_2     q1n_2 ...      p1(n-2)               q1(n-2)
          p0n_1     q0n_1 ...      p0(n-1)               q0(n-2)
          p0n_2     q0n_2 ...      p0(n-2)               q0(n-2)
          dp1n      d2p1n ...      dp1(n)/dtheta         d^2p1(n)/dtheta^2
          dp0n_1    dp1n_1 ...     dp1(n-1)/dtheta       dp1(n-1)/dtheta
          dq1n      d2q1n ...      dq1(n)/dtheta         d^2q1(n)/dtheta^2
          dq0n_1    dq1n_1 ...     dq0(n-1)/dtheta       dq1(n-1)/dtheta
          dp0n_2    dp1n_2 ...     dp0(n-2)/dtheta       dp1(n-2)/dtheta
          d2p0n_1   d2p1n_1 ...    d^2p0(n-1)/dtheta^2   d^2p1(n-1)/dtheta^2
          d2p0n_2   d2p1n_2 ...    d^2p0(n-2)/dtheta^2   d^2p1(n-2)/dtheta^2
          dq0n_2    dq1n_2 ...     dq0(n-2)/dtheta       dq1(n-2)/dtheta
          d2q0n_2   d2q1n_2...     d^2q0(n-2)/dtheta^2   d^2q1(n-2)/dtheta^2
          d2q0n_1   d2q1n_1 ...    d^2q0(n-1)/dtheta^2   d^2q1(n-1)/dtheta^2
          dp0n      dq0n ...       dp0(n)/dtheta         dq0(n)/dtheta
          d2p0n     d2q0n ...      d^2p0(n)/dtheta^2     d^2q0(n)/dtheta^2
		  ] = PQHeston(n,N,p0,q0)
      
% Use symbolic toolbox to obtain the polynomials 
% Expands up to N=7 terms

syms theta
thet = [theta^7 theta^6 theta^5 theta^4 theta^3 theta^2 theta^1 theta^0];

% Generate the coefficients p0n and q0n
for i=1:N
    pp(i,:) = sym(p0(i,:));
    qq(i,:) = sym(q0(i,:));
end

% Generate the polynomials p0n and q0n
p0n = dot(pp(n,:),thet);
q0n = dot(qq(n,:),thet);


%% Generate coefficients for p1n
syms pi11 pi12 pi13 pi14 pi15 pi16 pi17  % add more as needed
syms pi21 pi22 pi23 pi24 pi25 pi26 pi27  % add more as needed
syms pi31 pi32 pi33 pi34 pi35 pi36 pi37  % add more as needed
syms pi41 pi42 pi43 pi44 pi45 pi46 pi47  % add more as needed
syms pi51 pi52 pi53 pi54 pi55 pi56 pi57  % add more as needed
syms pi61 pi62 pi63 pi64 pi65 pi66 pi67  % add more as needed
syms pi71 pi72 pi73 pi74 pi75 pi76 pi77  % add more as needed

pis = [pi11 pi12 pi13 pi14 pi15 pi16 pi17; ...
       pi21 pi22 pi23 pi24 pi25 pi26 pi27; ...
       pi31 pi32 pi33 pi34 pi35 pi36 pi37; ...
       pi41 pi42 pi43 pi44 pi45 pi46 pi47; ...
       pi51 pi52 pi53 pi54 pi55 pi56 pi57; ...
       pi61 pi62 pi63 pi64 pi65 pi66 pi67; ...
       pi71 pi72 pi73 pi74 pi75 pi76 pi77];
syms pi

for i=1:N
    for k=1:N-i+2
        pi(i,k) = 0;
    end
    j = N-i+3;
    k = 1;
    while j<=N+1
        pi(i,j) = pis(i,k);
        k = k+1;
        j = j+2;
    end
end

% Generate the polynomial p1n
p1n = dot(pi(n,:)',thet);

%% Generate coefficients q1n

syms thetalong
syms qi11 qi12 qi13 qi14 qi15 qi16 qi17 qi18 qi19 qi110 qi111 qi112 qi113 qi114 qi115 qi116 qi117 % add more as needed
syms qi21 qi22 qi23 qi24 qi25 qi26 qi27 qi28 qi29 qi210 qi211 qi212 qi213 qi214 qi215 qi216 qi217 % add more as needed
syms qi31 qi32 qi33 qi34 qi35 qi36 qi37 qi38 qi39 qi310 qi311 qi312 qi313 qi314 qi315 qi316 qi317 % add more as needed
syms qi41 qi42 qi43 qi44 qi45 qi46 qi47 qi48 qi49 qi410 qi411 qi412 qi413 qi414 qi415 qi416 qi417 % add more as needed
syms qi51 qi52 qi53 qi54 qi55 qi56 qi57 qi58 qi59 qi510 qi511 qi512 qi513 qi514 qi515 qi516 qi517 % add more as needed
syms qi61 qi62 qi63 qi64 qi65 qi66 qi67 qi68 qi69 qi610 qi611 qi612 qi613 qi614 qi615 qi616 qi617 % add more as needed
syms qi71 qi72 qi73 qi74 qi75 qi76 qi77 qi78 qi79 qi710 qi711 qi712 qi713 qi714 qi715 qi716 qi717 % add more as needed

qis = [qi11 qi12 qi13 qi14 qi15 qi16 qi17 qi18 qi19 qi110 qi111 qi112 qi113 qi114 qi115 qi116 qi117; ...
       qi21 qi22 qi23 qi24 qi25 qi26 qi27 qi28 qi29 qi210 qi211 qi212 qi213 qi214 qi215 qi216 qi217; ...
       qi31 qi32 qi33 qi34 qi35 qi36 qi37 qi38 qi39 qi310 qi311 qi312 qi313 qi314 qi315 qi316 qi317; ...
       qi41 qi42 qi43 qi44 qi45 qi46 qi47 qi48 qi49 qi410 qi411 qi412 qi413 qi414 qi415 qi416 qi417; ...
       qi51 qi52 qi53 qi54 qi55 qi56 qi57 qi58 qi59 qi510 qi511 qi512 qi513 qi514 qi515 qi516 qi517; ...
       qi61 qi62 qi63 qi64 qi65 qi66 qi67 qi68 qi69 qi610 qi611 qi612 qi613 qi614 qi615 qi616 qi617; ...
       qi71 qi72 qi73 qi74 qi75 qi76 qi77 qi78 qi79 qi710 qi711 qi712 qi713 qi714 qi715 qi716 qi717];
syms qi

thetlong = [theta^16 theta^15 theta^14 theta^13 theta^12 theta^11 theta^10 theta^9 ...
            theta^8  theta^7  theta^6  theta^5  theta^4  theta^3  theta^2  theta^1 theta^0];

qi = qis;
for i=1:N
    for j=1:17;
        qi(i,j) = 0;
    end
end

qi(2,16) = qis(2,1);

qi(3,13) = qis(3,1);
qi(3,15) = qis(3,2);
qi(3,17) = qis(3,3);

qi(4,10) = qis(4,1);
qi(4,12) = qis(4,2);
qi(4,14) = qis(4,3);
qi(4,16) = qis(4,4);

qi(5,7) = qis(5,1);
qi(5,9) = qis(5,2);
qi(5,11) = qis(5,3);
qi(5,13) = qis(5,4);
qi(5,15) = qis(5,5);
qi(5,17) = qis(5,6);

qi(6,4) = qis(6,1);
qi(6,6) = qis(6,2);
qi(6,8) = qis(6,3);
qi(6,10) = qis(6,4);
qi(6,12) = qis(6,5);
qi(6,14) = qis(6,6);
qi(6,16) = qis(6,7);

qi(7,1) = qis(7,1);
qi(7,3) = qis(7,2);
qi(7,5) = qis(7,3);
qi(7,7) = qis(7,4);
qi(7,9) = qis(7,5);
qi(7,11) = qis(7,6);
qi(7,13) = qis(7,7);
qi(7,15) = qis(7,8);
qi(7,17) = qis(7,9);


% Generate the polynomial q1n
q1n = dot(qi(n,:)',thetlong);


%% Generate the past polynomials p0n-1 and p0n-2, q0n-1 and q0n-2, etc.
if n>=2
    p0n_1 = dot(pp(n-1,:) ,thet);
    q0n_1 = dot(qq(n-1,:) ,thet);
    p1n_1 = dot(pi(n-1,:)',thet);
    q1n_1 = dot(qi(n-1,:)',thetlong);
else
    p0n_1 = 0;
    q0n_1 = 0;
    p1n_1 = 0;
    q1n_1 = 0;
end
if n>=3
    p0n_2 = dot(pp(n-2,:) ,thet);
    q0n_2 = dot(qq(n-2,:) ,thet);
    p1n_2 = dot(pi(n-2,:)',thet);
    q1n_2 = dot(qi(n-2,:)',thetlong);
else
    p0n_2 = 0;
    q0n_2 = 0;
    p1n_2 = 0;
    q1n_2 = 0;
end

%% Required derivatives
% n
dp0n   = diff(p0n,theta);         % First order derivative of p0n
dq0n   = diff(q0n,theta);         % First order derivative of q0n
dp1n   = diff(p1n,theta);         % First order derivative of p1n
d2p1n  = diff(p1n,theta,2);       % Second order derivative of p1n
dq1n   = diff(q1n,theta);         % First derivative of q1n
d2q1n  = diff(q1n,theta,2);       % Second derivative of q1n
d2p0n  = diff(p0n,theta,2);       % Second derivative of p0n
d2q0n  = diff(q0n,theta,2);       % Second derivative of q0n

% n-1
dp0n_1  = diff(p0n_1,theta);      % First derivative of p0n-1
dp1n_1  = diff(p1n_1,theta);      % First derivative of p1n-1
dq0n_1  = diff(q0n_1,theta);      % First derivative of q0n-1
dq1n_1  = diff(q1n_1,theta);      % First derivative of q1n-1
d2p0n_1 = diff(p0n_1,theta,2);    % Second derivative of p0n-1
d2p1n_1 = diff(p1n_1,theta,2);    % Second derivative of p1n-1
d2q0n_1 = diff(q0n_1,theta,2);    % Second derivative of q0n-1
d2q1n_1 = diff(q1n_1,theta,2);    % Second derivative of q1n-1

% n-2
dp0n_2  = diff(p0n_2,theta);      % First derivative of p0n-2
dp1n_2  = diff(p1n_2,theta);      % First derivative of p1n-2
d2p0n_2 = diff(p0n_2,theta,2);    % Second derivative of p0n-2
d2p1n_2 = diff(p1n_2,theta,2);    % Second derivative of p1n-2

dq0n_2  = diff(q0n_2,theta);      % First derivative of q0n-2
dq1n_2  = diff(q1n_2,theta);      % First derivative of q1n-2
d2q0n_2 = diff(q0n_2,theta,2);    % Second derivative of q0n-2
d2q1n_2 = diff(q1n_2,theta,2);    % Second derivative of q1n-2

