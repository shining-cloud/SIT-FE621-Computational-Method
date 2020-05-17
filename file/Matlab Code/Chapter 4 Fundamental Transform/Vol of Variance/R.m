function y = R(p,q,X,Z,T)

% "R" function for Lewis (2000) vol of vol expansion

if p==2 & q==0
	y = T*(0.5*(X/Z)^2 - 0.5/Z - 1/8);
elseif p==1 & q==1
	y = -X/Z + 0.5;
elseif p==1 & q==2
	y = (X/Z)^2 - X/Z - 0.25/Z*(4-Z);
elseif p==2 & q==2
	y = T*(0.5*(X/Z)^4 - 0.5*(X/Z)^3 - 3*(X^2/Z^3) + 1/8*(X/Z^2)*(12+Z) + 1/32/Z^2*(48-Z^2));
end



