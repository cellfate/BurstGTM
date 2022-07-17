function [result,flag] = ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene)
%% Rejection Sampler for approximate Bayesian computaion
% Inputs:
%    N - the number of ABC posterior samples
%    prior - function that generates iid samples from the parameter joint
%        prior distribution
%    f - function that computes statictis given a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - the discrepancy acceptance threshold
%    T - number of rounds
%    gene - sequence number of the gene.
%
% Outputs:
%    result - a matrix of ABC prior samples

result = cell(N,T);
result0 = [];
total_time = 0;
flag = true;
while total_time < N && size(result0,2) < 10*N
    tic;
    % generate trial from the prior
    theta_trial = prior();
    param.kon = theta_trial(1);
    param.ron = theta_trial(2);
    param.koff = theta_trial(3);
    param.roff = theta_trial(4);
    param.mu = theta_trial(5);
    param.delta = theta_trial(6);
    % compute theorical statictis of parameters
    static_theo = f(param);
    dist = rho(static_theo);
    param.dist = dist;
    % accept or reject
    if dist <= epsilon
        ind = size(result0,2)+1;
        result0(ind).kon = param.kon;
        result0(ind).ron = param.ron;
        result0(ind).koff = param.koff;
        result0(ind).roff = param.roff;
        result0(ind).mu = param.mu;
        result0(ind).delta = param.delta;
        result0(ind).dist = param.dist;
%         disp(ind)
    end
    elapsedTime = toc;
    total_time = total_time + elapsedTime;
end

if total_time > 300
    flag = false;
    fprintf('Gene %d :wrong!!\n',gene);
    pause(5)
else
    [~,index0] = sort([result0.dist]);
    result0 = result0(index0(1:N));
    for index2 = 1:N
        result(index2,1) = {result0(index2)};
    end
end
end
