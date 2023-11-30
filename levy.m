%================ Levy flight vs Brownian motion in 2d ==================%
%                                                                        %
%   The particle is initially located at [0,0].                          %
%   Assumptions:                                                         %
%       * diffusion is homogeneous                                       %
%       * step size ~ L(s,alpha)                                         %
%         where L(s,alpha) = |s|^(-1 - alpha) for s large, alpha > 0     %
%                                                                        %
%   Written by: Xige Yang                                                %
%========================================================================%

clear all;
close all;
clc;

% Parameters
T = 80;           % total time
T0 = 0.;            % initial time
dt = .1;            % time step
alpha = 1;        
x1 = 0; y1 = 0;
x2 = 0; y2 = 0;
sigma = 1;       % Standard deviation of short-range diffusion
jump = 5;

% Initialize
path_b = [];
path_l = [];

% Add plots
figure;
s=get(gcf, 'Position');
s(3)=500;
s(4)=200;
set(gcf, 'Position', s);
% For Browian Motion
subplot(1,2,1);
plot(x1,y1);
xlabel('x -->');
ylabel('y -->');
title(sprintf('Brownian motion at time = %f', T0));
% For Random Walk
subplot(1,2,2)
plot(x2,y2);
xlabel('x -->');
ylabel('y -->');
title(sprintf('Random walk at time = %f', T0));

% Iteration to store positions of particles
for t = dt:dt:T
    %% Brownian motion part
    x_prev = x1; y_prev = y1;
    
    % updating positions
    step = sigma*randn(1,2);
    x1 = x_prev + step(1);
    y1 = y_prev + step(2);
    path_b = [path_b; norm(step)];
    
    % updating plots
    subplot(1,2,1);
    line([x_prev,x1],[y_prev,y1]);
    title(sprintf('Brownian motion at time = %f', t));
    
    %% Levy flight part
    x_prev = x2; y_prev = y2;
    theta = 2*pi*rand();    % moving direction, homogeneity assumption
    
    % creating levy distribution (changerable)
    top = gamma(1+alpha)*sin(pi*alpha/2);
    bottom = gamma((1+alpha)/2)*alpha*2^((alpha-1)/2);
    sigma_u = (top/bottom)^(1/alpha);
    u = sigma_u*randn(); 
    v = sigma*randn();
    step = u./(abs(v).^(1/alpha));
    path_l = [path_l; norm(step)];
    
    % updating position
    x2 = x_prev + step * cos(theta);
    y2 = y_prev + step * sin(theta);
    subplot(1,2,2);
    line([x_prev,x2],[y_prev,y2]);
    title(sprintf('Levy flight at time = %f', t));
    drawnow;
    
    saveInd = floor(t/dt);
    if mod(saveInd,10) == 0
        saveas(gcf,sprintf('Imgs_vs/Compare_%03d.png',floor(saveInd/10)));
    end
    
end

figure;
subplot(1,2,1);
hb = histogram(path_b);
title('Brownian motion steps')
subplot(1,2,2);
hl = histogram(path_l);
title('Levy flight steps')