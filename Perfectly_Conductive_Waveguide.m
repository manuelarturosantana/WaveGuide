%%%%% File to compute solution to wave equation inside a wave guide.

% Create x,y,omega grid points
xVals = linspace(0,pi,50);
yVals = linspace(0,20,500);
%Naming conventions follow Anderson paper
a = -10; b = 10; A = b-a / 2; delta = b + a / 2; P = 2 * A; alpha = P / (2 * pi); 
M = 100;  
omegaVals = linspace(-delta,delta,m);
times = linspace(0,40,1000);
solutions = zeros(length(xVals),length(yVals),length(times));

%Function to take the inverse transform of Assume c = 1,sigma=1
f = @(omega,y,n,t) exp(1i .* sqrt(omega.^2 - n^2) .* y) .* exp(-(1/4) .* (omega -2).^2) .* exp(-1i .* omega .* t);

% Numerically Compute inverse Fourier Transform at all points
%for t=1:length(times)
parfor (t = 1:length(times),4)
    v = zeros(length(xVals),length(yVals));
    for x = 1:length(xVals)
        for y = 1:length(yVals)
            for n = 1:1
                fVals = f(omegaVals + delta,yVals(y),n,times(t));
                c = fft(fVals) / M;
                %Assume sigma = 1, and 1/ 2 * pi factor comes in the inverse
                %transform
                %Calcuate the InverseFT via the trapezoidal method
                %v(x,y) = v(x,y) + sin(n *x) * (1 / (2 * sqrt(pi))) * trapz(omegaVals,fVals);
                
                %Calculate using the sinc method
                M_vals = -(M/2):(M/2)-1;
                for m = 1:length(c)
                    v(x,y) = v(x,y) + c(m) * 2 * A * sinc((2 * P / A) * alpha * times(t) - M_vals(m));
                end
                %Multiply by the exponential and the sin val
                v(x,y) = v(x,y) *sin(n * x) * exp(-1i * delta * times(t));
            end 
        end
    end
%     figure(t)
    solutions(:,:,t) = real(v);
%     imagesc(real(v),[0,max(solutions(:,:,1),[],'all')])
%     colorbar
end

vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
limits = [min(solutions(:,:,1),[],'all'),max(solutions(:,:,1),[],'all')];
%limits = [0,max(solutions(:,:,1),[],'all')];
for ind = 1:length(times)
   z = solutions(:,:,ind);
   z = [zeros(1000,length(yVals));z;zeros(1000,length(yVals))];
   mycolormap = customcolormap([linspace(0,1,7)], {'#042a8a','#0f0f0f','#008080','#0f0f0f','#ff8000','#0f0f0f','#eb0202'});
   imagesc(z,limits),colormap(mycolormap)
   colorbar
   drawnow
   F(ind) = getframe(gcf); 
   writeVideo(vidfile,F(ind));
end
close(vidfile)


