% 1D deblurring example. A given input signal xtrue is convolved with a
% Gaussian kerne and white noise is added. The problem is to estimate the
% input from the noisy blurred sample.

% Discretization of the interval [0,1] by N points

% input: 'boxcar', 'triangle', oscillation'
inputfunc = 'boxcar';

% prior: 'ord1fixed', 'ord1interval', 'ord1jumps','ord1falsejumps;
prior = 'ord1falsejumps';

N = 160;
t = linspace(1/N,1,N);

% Defining the Gaussian blurring kernel

width = 0.05;
a = 1/sqrt(2*pi*width^2)*exp(-0.5*(1/width^2)*t.^2);
A = 1/N*toeplitz(a);

% Plotting the kernel at two different locations

figure(1)
plot(t,A(10,:),'k-','LineWidth',2)
text(0.2,0.7*max(A(10,:)),'n=10')
hold on
plot(t,A(100,:),'r-','LineWidth',2)
text(0.8,0.7*max(A(100,:)),'n=100')
hold off

pause
% Defining different input signals to claculate the data

if strcmp(inputfunc,'boxcar')
    % boxcar input
    xtrue = zeros(N,1);
    n = 40;
    xtrue(50:50+n) = ones(n+1,1);
    figure(2)
    plot(t,xtrue,'k-','LineWidth',2)
    axis([0,1,-0.1,1.1*max(xtrue)])
elseif strcmp(inputfunc,'oscillation')
   % smooth oscillation input
   xtrue = 10*(t'-0.5).*exp(-0.5*1e3*(t'-0.5).^2);
   figure(2)
   plot(t,xtrue,'k-','LineWidth',2)
   axis([0,1,-0.2,0.2])
elseif strcmp(inputfunc,'triangle')
    % triangular pulse
    t0 = 0.7;
    alpha = 0.6;
    beta = alpha*t0/(1-t0);
    xtrue = alpha*t.*(t<=t0) + beta*(1-t).*(t>t0);
    xtrue = xtrue';
    figure(2)
    plot(t,xtrue,'k-','LineWidth',2)
    axis([0,1,-0.02,0.8])
end
hold on    

% Claculating the blurred noisy singnal. The standard deviation of the
% additive noise is defined percentagewise of the maximum of the noiseless
% signal.

noise = 4;       % Noise level in percentages of the max. of noiseless signal
y0 = A*xtrue;    % Noiseless signal
sigma = max(y0)*noise/100;                % STD of the additive noise
y = y0 + sigma*randn(N,1);

figure(2)
plot(t,y,'r-','LineWidth',2)
hold off

pause

% MAP solution with different Gaussian priors

% White noise prior

%STD = 0.1;    % Reasonable values: Boxcar: STD = 0.5, Smooth: 0.001 ... 0.01
%L = 1/STD*eye(N);

% First order smoothness prior, information about the variance included

if strcmp(prior,'ord1fixed')
    % First order prior, variance fixed
    L = eye(N) - diag(ones(N-1,1),-1);
    tau = 100;
    theta = tau*1/(4*N^2);
    D12 = 1/sqrt(theta)*eye(N);
elseif strcmp(prior,'ord1interval')
    % First order, varince larger over an interval
    L = eye(N) - diag(ones(N-1,1),-1);
    theta = 1e-6*ones(N,1);
    Iint = find((t>0.4)&(t<0.6));
    theta(Iint) = 25/(N^2)*ones(length(Iint),1);
    D12 = diag(1./sqrt(theta));
elseif strcmp(prior,'ord1jumps')
    % First order, information about the whereabouts of the jumps
    L = eye(N) - diag(ones(N-1,1),-1);
    theta = 1e-4*ones(N,1);
    theta(48:52) = 0.5*ones(5,1);
    theta(88:92) = 0.5*ones(5,1);
    D12 = diag(1./sqrt(theta));
    elseif strcmp(prior,'ord1extrajumps')
    % First order, information about the whereabouts of the jumps. An extra
    % false jump assumed
    L = eye(N) - diag(ones(N-1,1),-1);
    theta = 1e-4*ones(N,1);
    theta(48:52) = 0.5*ones(5,1);
    theta(88:92) = 0.5*ones(5,1);
    theta(108:112) = 0.5*ones(5,1);
    D12 = diag(1./sqrt(theta));
elseif strcmp(prior,'ord1falsejumps')
    % First order, incorrect information about the jumps
    L = eye(N) - diag(ones(N-1,1),-1);
    theta = 1e-4*ones(N,1);
    theta(38:42) = 0.5*ones(5,1);
    theta(98:102) = 0.5*ones(5,1);
    theta(68:72) = 0.5*ones(5,1);
    D12 = diag(1./sqrt(theta));
end


% Five random draws from the prior distribution

w = randn(N,5);
xdraw = L\(D12\w);
figure(3)
plot(t,xdraw, 'k-','LineWidth',2)

pause

% Calculating the MAP estimate and posterior variances

DL = D12*L;
xmean = [(1/sigma)*A;DL]\[(1/sigma)*y;zeros(N,1)];
Gamma = inv((1/sigma^2)*A'*A + DL'*DL); 

% Plotting the MAP estimate and the 2*STD envelope

% Defining different shades of blue for plotting
 shades = [176 224 230;
             135 206 235;
             135 206 255;
             126 192 238;
             108 166 205];
   shades = 1/255*shades;

STD = sqrt(diag(Gamma));
xhigh = xmean + 2*STD;
xlow = xmean - 2*STD;

% Adding the known zero initial value

t = [0;t'];
xlow = [0;xlow];
xhigh = [0;xhigh];
xmean = [0;xmean];
xtrue = [0;xtrue];

figure(4)
fill([t;t(N+1:-1:1)],[xlow;xhigh(N+1:-1:1)],shades(1,:))
hold on
plot(t,xmean,'r-','LineWidth',2)
plot(t,xtrue,'k-','LineWidth',1)
hold off

