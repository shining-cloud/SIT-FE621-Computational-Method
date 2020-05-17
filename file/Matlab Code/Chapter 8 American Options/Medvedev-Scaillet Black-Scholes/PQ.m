function [p0n   q0n ...
          p1n   q1n ...
          p1n_1 q1n_1 ...
          p1n_2 q1n_2 ...
          p0n_1 q0n_1 ...
          p0n_2 q0n_2 ...
          dp1n  d2p1n ...
          dp0n_1 dp1n_1 ...
          dq1n    d2q1n ...
          dq0n_1  dq1n_1] = PQ(n,N,p0,q0)

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
syms qi11 qi12 qi13 qi14 qi15 qi16 qi17  % add more as needed
syms qi21 qi22 qi23 qi24 qi25 qi26 qi27  % add more as needed
syms qi31 qi32 qi33 qi34 qi35 qi36 qi37  % add more as needed
syms qi41 qi42 qi43 qi44 qi45 qi46 qi47  % add more as needed
syms qi51 qi52 qi53 qi54 qi55 qi56 qi57  % add more as needed
syms qi61 qi62 qi63 qi64 qi65 qi66 qi67  % add more as needed
syms qi71 qi72 qi73 qi74 qi75 qi76 qi77  % add more as needed

qis = [qi11 qi12 qi13 qi14 qi15 qi16 qi17; ...
       qi21 qi22 qi23 qi24 qi25 qi26 qi27; ...
       qi31 qi32 qi33 qi34 qi35 qi36 qi37; ...
       qi41 qi42 qi43 qi44 qi45 qi46 qi47; ...
       qi51 qi52 qi53 qi54 qi55 qi56 qi57; ...
       qi61 qi62 qi63 qi64 qi65 qi66 qi67; ...
       qi71 qi72 qi73 qi74 qi75 qi76 qi77];
syms qi

for i=1:N;
    for k=1:N-i+3;
        qi(i,k) = 0;
    end
    j = N-i+4;
    k = 1;
    while j<=N+1
        qi(i,j) = qis(i,k);
        k = k+1;
        j = j+2;
    end
    qi = qi(:,1:N+1);
end

% Generate the polynomial q1n
q1n = dot(qi(n,:)',thet);

%% Generate the past polynomials p0n-1 and p0n-2, q0n-1 and q0n-2, etc.
if n>=2
    p0n_1 = dot(pp(n-1,:) ,thet);
    q0n_1 = dot(qq(n-1,:) ,thet);
    p1n_1 = dot(pi(n-1,:)',thet);
    q1n_1 = dot(qi(n-1,:)',thet);
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
    q1n_2 = dot(qi(n-2,:)',thet);
else
    p0n_2 = 0;
    q0n_2 = 0;
    p1n_2 = 0;
    q1n_2 = 0;
end

%% Required derivatives
dp1n   = diff(p1n,theta);         % First order derivative of p1n
d2p1n  = diff(dp1n,theta);        % Second order derivative of p1n
dp0n_1 = diff(p0n_1,theta);       % First derivative of p0n-1
dp1n_1 = diff(p1n_1,theta);       % First derivative of p1n-1

dq1n   = diff(q1n,theta);          % First derivative of q1n
d2q1n  = diff(dq1n,theta);         % Second derivative of q1n
dq0n_1 = diff(q0n_1,theta);        % First derivative of q0n-1
dq1n_1 = diff(q1n_1,theta);        % First derivative of q1n-1
