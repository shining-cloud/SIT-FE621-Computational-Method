% Illustration of the put payoff

clc; clear;

% Real and Imaginary axes
R = [-4:.1:4];
I = [-.9:.01:-.1];

% Strike price
K = 50;

%% Generate the real part of the put payoff
for x=1:length(R);
    for y=1:length(I);
        kr = R(x);
        ki = I(y);
        k = kr + i*ki;
        w(x,y) = real(-K^(i*k+1)/(k^2 - i*k));
    end
end

%% Generate the surface
surf(w,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
view([-300,130,70])
axis tight
camlight left
set(gca,'XTickLabel', [-0.9:.5:-0.1]);  % Imaginary
set(gca,'YTickLabel', [-4:1:4]);        % Real

numYticks = 5;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),numYticks))
numXticks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),numXticks))

xlabel('Imag(k)')
ylabel('Real(k)')



