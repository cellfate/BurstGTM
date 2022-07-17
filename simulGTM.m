function [x,t] = simulGTM(param)
%% This code simulates scRNA-seq data under the GTM.
% Inputs:
%    param: A structure contains model parameters 
%    param.kon and param.ron: OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon-1) * e^(-ron * t) / gamma(kon)
%    param.koff and param.roff: ON dwell time distribution f_on(t) = roff^(koff) * t^(koff-1) * e^(-roff * t) / gamma(koff)
%    param.mu: Transcriptional rate
%    param.delta: Degradation rate
%    param.x0; The initial point of the simulation
%    param.tottime: The total time of the simulation
%
% Outputs:
%    x - Time series of [OFF ON mRNA]
%    t - Corresponds to time of x

% parameters setting
koff = param.koff;
roff = param.roff;
kon = param.kon;
ron = param.ron;
mu = param.mu;
delta = param.delta;
x = param.x0;
tottime = param.tottime;

% reaction matrix
r_mu  = [-1  1  0;
    1 -1  0;
    0  0  1;
    0  0 -1];

% time
t = 0;

% simulation the reaction
while (t(end) < tottime)
    if x(end,1) == 1 % OFF state
        tau_off = gamrnd(kon,1/ron);
        t_temp = t(end) + tau_off; % OFF dwell time
        while t(end) < t_temp % degration reaction only
            a_0 = delta * x(end,3);
            r1 = rand;
            tau_1 = (1/a_0) * log(1/r1);
            t = [t;t(end) + tau_1];
            x = [x;x(end,:) + r_mu(4,:)];
        end
        t = [t(1:end-1);t_temp];
        x = [x(1:end-1,:);x(end-1,:) + r_mu(1,:)];% OFF -> ON

    elseif x(end,1) == 0 % ON state
        tau_on = gamrnd(koff,1/roff);
        t_temp = t(end) + tau_on; % OFF dwell time
        while t(end) < t_temp % generation and degration reaction
            a_mu(1) = mu * x(end,2);
            a_mu(2) = delta * x(end,3);
            a_0  = sum(a_mu);
            r1  = rand;
            tau_2 = (1/a_0) * log(1/r1);

            r2 = rand;
            for iter = 1:2
                if (sum(a_mu(1:iter)) >= r2*a_0)
                    next_mu = iter;
                    break;
                end
            end
            t = [t;t(end) + tau_2];
            x = [x;x(end,:) + r_mu(next_mu + 2,:)];
        end
        t = [t(1:end-1);t_temp];
        x = [x(1:end-1,:);x(end-1,:) + r_mu(2,:)];
    end
end

