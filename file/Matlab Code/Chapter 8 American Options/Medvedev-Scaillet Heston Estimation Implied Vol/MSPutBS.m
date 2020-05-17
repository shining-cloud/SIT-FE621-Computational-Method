function Price = MSPutBS(y,theta,K,sigma,r,q,T)

% Black-Scholes American put price using the Medvedev-Scaillet
% approximation with 5 terms
% To obtain the coefficients in this function, run the Matlab file
% Generating_PQ_coefficients.m, which uses the symbolic calculator in
% Matlab to obtain the coefficients in Appendix A and B of their paper.

mu = r-q;

% The "C" coefficients evaluated at theta = y
cdf = normcdf(y);
pdf = normpdf(y);
C1 = sigma*y*K/(y*cdf+pdf) ;
C2 = -1/2*(cdf*sigma^2*C1-2*mu*C1*cdf+sigma^3*y^2*K)/sigma/(cdf*y^2+cdf+y*pdf) ;
C3 = 1/24*(-24*y*cdf*sigma^3*C2+48*y*cdf*sigma*C2*mu+24*y*cdf*sigma^2*r*C1-24*pdf*sigma^3*C2+48*pdf*C2*mu*sigma+24*pdf*r*C1*sigma^2-3*pdf*sigma^4*C1+12*pdf*C1*sigma^2*mu-12*pdf*C1*mu^2+4*sigma^5*y^3*K)/sigma^2/(cdf*y^3+3*y*cdf+pdf*y^2+2*pdf) ;
C4 = -1/48*(-48*cdf*sigma^3*y^2*r*C2+72*cdf*sigma^4*y^2*C3-144*cdf*sigma^2*y^2*C3*mu-48*cdf*sigma^3*r*C2+72*cdf*sigma^4*C3-144*cdf*sigma^2*C3*mu+12*cdf*sigma^5*C2-48*cdf*sigma^3*C2*mu-24*cdf*sigma^4*r*C1+48*cdf*sigma*C2*mu^2+48*cdf*sigma^2*mu*r*C1+72*y*pdf*sigma^4*C3-144*y*pdf*C3*mu*sigma^2-48*y*pdf*r*C2*sigma^3-y*pdf*sigma^6*C1+6*y*pdf*sigma^4*C1*mu-12*y*pdf*sigma^2*C1*mu^2+8*y*pdf*C1*mu^3+2*sigma^7*y^4*K)/sigma^3/(cdf*y^4+6*cdf*y^2+3*cdf+pdf*y^3+5*y*pdf) ;
C5 = 1/1920*(16*sigma^9*y^5*K+5*pdf*sigma^8*C1-80*pdf*sigma^7*C2-7680*pdf*sigma^5*C4+40*pdf*y^2*sigma^6*C1*mu-120*pdf*y^2*sigma^4*C1*mu^2+160*pdf*y^2*sigma^2*C1*mu^3+7680*pdf*y^2*C4*sigma^3*mu+1920*pdf*y^2*r*sigma^4*C3+1920*pdf*r*sigma^5*C2-960*pdf*sigma^3*C2*mu^2+480*pdf*sigma^5*C2*mu+240*pdf*sigma^6*r*C1-5760*pdf*sigma^2*C3*mu^2+15360*pdf*C4*sigma^3*mu+640*pdf*sigma*C2*mu^3+5760*pdf*sigma^4*C3*mu+3840*pdf*r*sigma^4*C3-960*pdf*sigma^4*mu*r*C1+960*pdf*mu^2*r*C1*sigma^2-40*pdf*sigma^6*C1*mu+120*pdf*sigma^4*C1*mu^2-160*pdf*sigma^2*C1*mu^3-960*pdf*r^2*sigma^4*C1+80*pdf*C1*mu^4-1440*pdf*sigma^6*C3-3840*pdf*r*sigma^3*C2*mu+1920*y^3*cdf*sigma^4*r*C3+5760*y*cdf*sigma^4*C3*mu-5760*y*cdf*sigma^2*C3*mu^2-3840*y^3*cdf*sigma^5*C4-11520*y*cdf*sigma^5*C4-1440*y*cdf*sigma^6*C3-3840*pdf*y^2*sigma^5*C4-5*pdf*y^2*sigma^8*C1+7680*y^3*cdf*sigma^3*C4*mu+5760*y*cdf*sigma^4*r*C3+23040*y*cdf*sigma^3*C4*mu+1920*y*cdf*sigma^5*r*C2-3840*y*cdf*sigma^3*r*C2*mu-960*y*cdf*sigma^4*r^2*C1-80*pdf*y^2*C1*mu^4)/sigma^4/(cdf*y^5+10*cdf*y^3+15*y*cdf+pdf*y^4+9*pdf*y^2+8*pdf) ;

% Set 1 polynomials
p01 = theta ;
p11 = 0 ;
q01 = 1 ;
q11 = 0 ;

% Set 2 polynomials
p02 = theta^2+1 ;
p12 = -1/2*(-sigma^2+2*mu)*C1/sigma ;
q02 = theta ;
q12 = 0 ;

% Set 3 polynomials
p03 = theta^3+3*theta ;
p13 = -theta*(-sigma^2*C2+2*C2*mu+r*C1*sigma)/sigma ;
q03 = theta^2+2 ;
q13 = -1/8*(-8*sigma^3*C2+16*C2*mu*sigma+8*r*C1*sigma^2-sigma^4*C1+4*sigma^2*mu*C1-4*C1*mu^2)/sigma^2 ;

% Set 4 polynomials
p04 = theta^4+6*theta^2+3 ;
p14 = -1/2*(2*r*C2*sigma-3*C3*sigma^2+6*C3*mu)*theta^2/sigma+1/4*(-4*r*C2*sigma^2+6*C3*sigma^3-12*C3*sigma*mu+sigma^4*C2-4*sigma^2*C2*mu-2*sigma^3*r*C1+4*C2*mu^2+4*mu*r*C1*sigma)/sigma^2 ;
q04 = theta^3+5*theta ;
q14 = -1/48*(-72*sigma^4*C3+144*C3*mu*sigma^2+48*r*C2*sigma^3+sigma^6*C1-6*sigma^4*C1*mu+12*sigma^2*C1*mu^2-8*C1*mu^3)/sigma^3*theta ;

% Set 5 polynomials
p05 = theta^5+10*theta^3+15*theta ;
p15 = -(r*C3*sigma-2*C4*sigma^2+4*C4*mu)/sigma*theta^3+1/4*(-12*r*C3*sigma^2+24*C4*sigma^3-48*C4*sigma*mu-4*r*C2*sigma^3+8*r*sigma*C2*mu+2*r^2*sigma^2*C1+3*sigma^4*C3-12*C3*mu*sigma^2+12*C3*mu^2)/sigma^2*theta ;
q05 = theta^4+9*theta^2+8 ;
q15 = -1/384*(-sigma^8*C1+8*sigma^6*C1*mu-24*sigma^4*C1*mu^2+32*sigma^2*C1*mu^3-16*C1*mu^4-768*sigma^5*C4+1536*sigma^3*C4*mu+384*r*sigma^4*C3)/sigma^4*theta^2+1/384*(192*sigma^4*mu*r*C1-192*mu^2*r*C1*sigma^2+16*sigma^7*C2+192*sigma^3*C2*mu^2-128*sigma*C2*mu^3+192*r^2*sigma^4*C1-96*sigma^5*C2*mu-48*sigma^6*r*C1+768*r*sigma^3*C2*mu+288*sigma^6*C3-sigma^8*C1-16*C1*mu^4+1536*sigma^5*C4-768*r*sigma^4*C3-1152*sigma^4*C3*mu-384*r*sigma^5*C2+8*sigma^6*C1*mu-24*sigma^4*C1*mu^2+32*sigma^2*C1*mu^3+1152*sigma^2*C3*mu^2-3072*sigma^3*C4*mu)/sigma^4 ;

% The Black-Scholes American put approximation
cdf = normcdf(theta);
pdf = normpdf(theta);

Price = (C1*(p01*cdf + q01*pdf) + p11*cdf + q11*pdf)*T^(1/2) ...
      + (C2*(p02*cdf + q02*pdf) + p12*cdf + q12*pdf)*T^(2/2) ...
      + (C3*(p03*cdf + q03*pdf) + p13*cdf + q13*pdf)*T^(3/2) ...
      + (C4*(p04*cdf + q04*pdf) + p14*cdf + q14*pdf)*T^(4/2) ...
      + (C5*(p05*cdf + q05*pdf) + p15*cdf + q15*pdf)*T^(5/2);
