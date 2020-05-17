function y = BisecMSIV(S,K,rf,q,T,a,b,MktPrice,Tol,MaxIter)

lowCdif  = MktPrice - MSPriceBS(S,K,T,a,rf,q);
highCdif = MktPrice - MSPriceBS(S,K,T,b,rf,q);

if lowCdif*highCdif > 0
	y = -1;
else
	for x=1:MaxIter
		midP = (a+b)/2;
        midCdif = MktPrice - MSPriceBS(S,K,T,midP,rf,q);
		if abs(midCdif)<Tol
			break
		else
			if midCdif>0
				a = midP;
			else
				b = midP;
			end
		end
	end
	y = midP;
end


