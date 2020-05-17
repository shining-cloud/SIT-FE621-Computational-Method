function xhat = FRFT(x,beta)

% Fractional fast Fourier transform
% Uses the Matlab functions fft() and ifft() for the FFT and inverse FFT
% Input:
%   x = a row vector of values
%   beta = increment parameter

N = length(x);

% Construct the y and z vectors
y = [exp(-i.*pi.*(0:N-1).^2.*beta).*x, zeros(1,N)];
z = [exp( i.*pi.*(0:N-1).^2.*beta)   , exp(i.*pi.*(N:-1:1).^2.*beta)];

% FFT on y and z
Dy = fft(y);
Dz = fft(z);

% h vectors
h = Dy.*Dz;
ih = ifft(h);

% e vector
e = [exp(-i.*pi.*(0:N-1).^2.*beta), zeros(1,N)];

% The FRFT vector
xhat = e.*ih;
xhat = xhat(1:N);
