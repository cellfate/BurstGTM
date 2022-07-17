%% Estimate the parameters of synthesized data for a demo
clear;clc;
addpath(genpath(pwd))

% Parameters setting
% OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon-1) * e^(-ron * t) / gamma(kon)
% ON dwell time distribution f_on(t) = roff^(koff) * t^(koff-1) * e^(-roff * t) / gamma(koff)
param_true.kon = 3;
param_true.ron = 0.5;
param_true.koff = 2;
param_true.roff = 0.5;
param_true.mu = 20;
param_true.delta = 1;
param_true.x0 = [1,0,0];
param_true.tottime = 2000;

% Simulation algorithm for GTM.
[x,t] = simulGTM(param_true);
tq = 1500:0.1:param_true.tottime;
xq = interp1(t,x(:,3),tq,'previous');
data = xq;

figure
histogram(data,'normalization','pdf')

% Inference algorithm
statis_data = statisData(data);
statis_ther = statisGTM(param_true,4);
eps0 = sqrt(sum(log(statis_data./statis_ther).^2));
rho = @(s) sqrt(sum(log(statis_data./s).^2));
f = @(k) statisGTM(k,4);
N = 1000;
T = 6;
epsilon = 1;
prior = @() [5*rand(),logunif(-1,1),5*rand(),logunif(-1,1),50*rand,1];
proposal_sigma = 0.2;
proposal = @(x) lognrnd(x,proposal_sigma);
proposal_pdf = @(kon_post,kon_prior,ron_post,ron_prior,koff_post,koff_prior,roff_post,roff_prior,mu_post,mu_prior)...
    lognpdf(mu_post,log(mu_prior),proposal_sigma) * lognpdf(kon_post,log(kon_prior),proposal_sigma) *...
    lognpdf(ron_post,log(ron_prior),proposal_sigma) * lognpdf(koff_post,log(koff_prior),proposal_sigma) *...
    lognpdf(roff_post,log(roff_prior),proposal_sigma);
[result,~] = ABCSMCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,1);
figureResult(data,result(:,end),param_true)

save(sprintf('results/example/%d_%.1f_%d_%.1f_%d_%d.mat',param_true.kon,param_true.ron,param_true.koff,param_true.roff,param_true.mu))