 function y = HestonDE(NG,NP,CR,F,Hi,Lo,S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule)

% Differential Algorithm for Heston parameter estimation
% NG = Number of generations (iterations)
% NP = Number of population members
% CR = Crossover ratio (=0.5)
% F  = Threshold (=0.8)
% Hi = Vector of upper bounds for the parameters
% Lo = Vector of lower bounds for the parameters
% S  = Spot Price
% K1 = First strike point for the FRFT
% rf = Risk free rate
% q  = Dividend Yield
% MktPrice = Market quotes for prices
% K = Vector of Strikes
% T = Vector of Maturities
% PutCall = Matrix of 'P'ut or 'C'all
% MktIV = Market quotes for implied volatilities
% ObjFun = Type of Objective Function
%    1 = MSE
%    2 = RMSE
%    3 = IVMSE
%    4 = Christoffersen, Jacobs, Heston (2009)
% Bisection method settings
%    a = Lower limit
%    b = upper limit
%    Tol = Tolerance
%    MaxIter = Max number of iterations
% trap = 1 is "Little Trap" c.f., 0 is Heston c.f.
% FRFT Settings
%    N = Number of points
%    uplimit = Upper integration point
%    eta = integration grid size
%    alpha = dampening factor
%    rule = 'T'rapezoidal or 'S'impsons


kappaU = Hi(1);  kappaL = Lo(1);
thetaU = Hi(2);  thetaL = Lo(2);
sigmaU = Hi(3);  sigmaL = Lo(3);
v0U    = Hi(4);  v0L    = Lo(4);
rhoU   = Hi(5);  rhoL   = Lo(5);

% Step1.  Generate the population matrix of random parameters
P = [kappaL + (kappaU-kappaL)*rand(1,NP); ...
      thetaL + (thetaU-thetaL)*rand(1,NP); ...
	  sigmaL + (sigmaU-sigmaL)*rand(1,NP); ...
	     v0L + (   v0U-   v0L)*rand(1,NP); ...
	    rhoL + (  rhoU-  rhoL)*rand(1,NP)];

% Generate the random numbers outside the loop
U = rand(5,NP,NG);

% Loop through the generations
for k=1:NG
	% Loop through the population
	for i=1:NP
		% Select the i-th member of the population
		P0 = P(:,i);
        % Set the condition so that the "while" statement 
        % is executed at least once
        Pnew = -ones(1,5);
        Condition = sum(Lo < Pnew & Pnew < Hi);
        while Condition < 5
            % Select random indices for three other distinct members
            I = randperm(NP);
            I = I(find(I~=i));
            r = I(1:3);
            % The three distinct members of the population
            Pr1 = P(:,r(1));
            Pr2 = P(:,r(2));
            Pr3 = P(:,r(3));
            R = randperm(5);
            Pnew = zeros(1,5);
            % Steps 2 and 3.  Mutation and recombination
            for j=1:5
                Ri = R(1);
                u = U(j,i,k);
                if u<=CR || j==Ri
                    Pnew(j) = Pr1(j) + F*(Pr2(j) - Pr3(j));
                else
                    Pnew(j) = P0(j);
                end
            end
            Condition = sum(Lo < Pnew & Pnew < Hi);
        end
        % Step 4.  Selection
        % Calculate the objective function for the i-th member and for the candidate
        f0   = HestonObjFunFRFT(P0,  S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule);
        fnew = HestonObjFunFRFT(Pnew,S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule);
  		% Verify whether the candidate should replace the i-th member
		% in the population and replace if conditions are satisfied
		if fnew < f0
			P(:,i) = Pnew;
		end
	end
  	disp(['Completed loop for generation ' num2str(k)])
end

% Calculate the objective function for each member in the updated population
for i=1:NP
    f(i) = HestonObjFunFRFT(P(:,i),S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule);
end

% Find the member with the lowest objective function
J = find(f==min(f));

% Return the selected member/parameter
y = P(:,J)';

