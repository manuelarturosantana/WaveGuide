%%%%% File to compute solution to wave equation inside a wave guide.

% Create x,y,omega grid points
xVals = linspace(0,pi,100);
yVals = linspace(0,1,100);
omegaVals = linspace(-5,5,100);
times = linspace(0,5,200);
solutions = zeros(length(xVals),length(yVals),length(times));

% Comute sine series of g(x)
% For simplicity lets assume that g(x) is a sum of sin functions.
%g = @(x,n) sin(n * x);
%Function to take the inverse transform of Assume c = 1,sigma=1
f = @(omega,y,n,t) exp(1i .* sqrt(omega.^2 - n^2) .* y) .* exp(-pi^2 .* omega.^2) .* exp(-1i .* omega .* t);

% Numerically Compute inverse Fourier Transform at all points
for t = 1:length(times)
    v = zeros(length(xVals),length(yVals));
    for y = 1:length(yVals)
        for n = 1:1
            fVals = f(omegaVals,yVals(y),n,times(t));
            %Assume sigma = 1, and 1/ 2 * pi factor comes in the inverse
            %transform
            v(:,y) = v(:,y) + (1 / (2 * sqrt(pi))) * trapz(omegaVals,fVals);
       end       
    end
%     figure(t)
    solutions(:,:,t) = real(v);
%     imagesc(real(v),[0,max(solutions(:,:,1),[],'all')])
%     colorbar
end

vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
limits = [0,max(solutions(:,:,1),[],'all')];
for ind = 1:length(times)
   z = solutions(:,:,ind);
   imagesc(z,limits),colormap(hot) 
   drawnow
   F(ind) = getframe(gcf); 
   writeVideo(vidfile,F(ind));
end
close(vidfile)


