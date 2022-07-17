function [result,flag] = ABCSMCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,gene)
%% Sequential Monte Carlo Sampler for approximate Bayesian computaion
% Inputs:
%    N - the number of particles
%    prior - prior distribution sampler
%    f - function that generates simulated data give a parameters set
%    rho - discrepancy metric, treated as a function of simulated data only
%    epsilon - a sequence of discrepancy acceptance thresholds
%    T - number of rounds
%    proposal - proposal kernel process, generates new trials
%    proposal_pdf - the PDF of the proposal kernel
%    gene - sequence number of the gene.
%
% Outputs:
%    result - the sequence of particles

% initialise
[result,flag] = ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene);
W = (1/N)*ones(N,T);

% sequential sampling
if flag == true
    for t = 2:T
        epsilon = prctile(cellfun(@(c)c.dist,result(:,t-1)),50);
%         fprintf('The eps of %d rounds is %f\n',t,epsilon);
        % generate new generation of particles
        for i = 1:N
            % rejections steps
            r = inf;
            total_time = 0;
            while total_time < 20 && r > epsilon
                tic;
                % sample from weighted particles
                j = randsample(N,1,true,W(:,t-1));
                param_temp = result{j,t-1};
                
                % generate new particle based on proposal kernel
                param_proposal = proposal(log([param_temp.kon,param_temp.ron,param_temp.koff,param_temp.roff,param_temp.mu]));
                param.kon = param_proposal(1);
                param.ron = param_proposal(2);
                param.koff = param_proposal(3);
                param.roff = param_proposal(4);
                param.mu = param_proposal(5);
                param.delta = 1;
                result{i,t} = param;
                static_temp = f(param);
                r = rho(static_temp);
                result{i,t}.dist = r;
                elapsedTime = toc;
                total_time = total_time + elapsedTime;
            end
            
            if total_time > 20
                flag = false;
                break;
            end
            
            % recompute particle weight using optimal backward kernel
            back_K = 0;
            for j = 1:N
                back_K = back_K + W(j,t-1)*proposal_pdf(result{i,t}.kon,result{j,t-1}.kon,...
                    result{i,t}.ron,result{j,t-1}.ron,...
                    result{i,t}.koff,result{j,t-1}.koff,...
                    result{i,t}.roff,result{j,t-1}.roff,...
                    result{i,t}.mu,result{j,t-1}.mu);
            end
            W(i,t) = (1/5*1/5*1/30*1/(4*result{i,t}.ron * result{i,t}.roff))/back_K;
        end
        
        if flag == false
            break;
        end
        
        % resample
        if t < T
            result_rs  = result(:,t);
            W(:,t) = W(:,t)./sum(W(:,t));
            J = randsample(N,N,true,W(:,t));
            result(:,t) = result_rs(J);
        end
        
        %   re-set weights
        W(:,t) = 1/N;
    end
end
