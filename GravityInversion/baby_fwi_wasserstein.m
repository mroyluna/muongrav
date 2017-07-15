function void = baby_fwi_wasserstein(n_amp)
clf
t = linspace(0,10,301)';
m = 0;
g = fwd_model(m,t);
g = g + n_amp*2*(rand(size(t))-0.5);

for method = [1 2]
    m = 1.6;
    gam = 1e-5;
    grad = 10;
    delta = 1e-6;
    i = 1;
    while ((i < 50) && gam > 0)
        f = fwd_model(m,t);
        fp = fwd_model(m+delta,t);
        if(method == 1)
            d = compute_Wasserstein_metric(f,g,t);
            dp = compute_Wasserstein_metric(fp,g,t);    
        else
            % L-2
            d = sum((f-g).^2); 
            dp = sum((fp-g).^2); 
        end
        grad = (dp-d)/delta;
        if (i == 1)
        else
            % Slightly damped (and regularized) step size update
            gamold = gam;
            gam = 0.75*abs((m - mold)*(grad-gradold)/(abs(grad-gradold)^2+1e-4));
        end
        mold = m;
        gradold = grad;
        m = m - gam*grad;
        % Constrain the parameter in the model
        if(abs(m) > 5)
            m = mold - 0.1*gam*grad;
        end
        i = i+1;
        figure(1)
        subplot(2,2,method)
        plot(t,g,'r',t,f,'--','linewidth',2)
        set(gca,'fontsize',18)
        xlabel('t')
        legend('g(t)','f(t)')
        figure(1)
        subplot(2,2,2+method)
        plot(i,m,'r*')
        hold on
        drawnow
        pause
    end
end

function rec = fwd_model(m,t);
% This is the forward model, it returns the data at a reciever at
% the times specified in the vector t. 

% Currently, the model parameter m is the phase shift of the train of Gaussians

centers = [3 5 7];
widths = [0.3 0.4 0.2];
centers  = centers + m;
rec = zeros(size(t));
for i = 1:3
    rec = rec + exp(-((t-centers(i))/widths(i)).^2);
end
rec = rec + 1;


function dw1 = compute_Wasserstein_metric(f,g,t);
h = t(2)-t(1);
n = length(t);
% When multiplied onto a vector this matrix will 
% return an approximation to the CDF of that vector at 
% the times t+h.
IM = h*tril(ones(n,n));
G = IM*g;
F = IM*f;
% We normalize the signals to have unit mass.
g = g/G(end);
f = f/F(end);
G = G/G(end);
F = F/F(end);
% We perform inverse iteration (i.e. we swap the roles of (x,y)) to
% find the inverse CDF.

FI = interp1(F,t,F,'linear','extrap');
GI = interp1(G,t,F,'linear','extrap');

% Finally we sum up the Wasserstein metric 
dw1 = sum((GI-FI).^2);
