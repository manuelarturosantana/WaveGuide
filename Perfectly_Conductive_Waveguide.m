%%%%% File to compute solution to wave equation inside a wave guide.

% Create x,y,omega grid points
xVals = linspace(0,pi,50);
yVals = linspace(0,0.1,50);
omegaVals = linspace(-5,5,100);

% Comute sine series of g(x)
% For simplicity lets assume that g(x) is a sum of sin functions.
g = @(x,n) sin(n * x);
%Function to take the inverse transform of Assume c = 1,sigma=1
f = @(omega,y,n,t) exp(1i .* sqrt(omega.^2 - n^2) .* y) .* exp(pi^2 .* omega.^2) .* exp(-1i .* omega .* t);

% Numerically Compute inverse Fourier Transform at all points
v = zeros(50,50);
for x = 1:50
    for y = 1:50
        for n = 1:20
            fVals = f(omegaVals,y,n,0.01);
            %Assume sigma = 1, and 1/ 2 * pi factor comes in the inverse
            %transform
            v(x,y) = v(x,y) + (1 / (2 * sqrt(pi))) * g(xVals(x),n) * trapz(omegaVals,fVals);
       end       
    end
end
imagesc(real(v))