function theta = MLInference(vals,theta0)
vals = vals(:);
% Parameters range
lb = [1e-3;1e-3;1];
ub = [1e3;1e3;1e4];
A = [];
b = [];
Aeq = [];
beq = [];

theta = fmincon(@(theta) LogLikelihood(theta,vals),theta0,A,b,Aeq,beq,lb,ub);

% nonlcon = [];
% options = optimoptions(@fmincon,'HessianApproximation','lbfgs');
% theta = fmincon(@(theta) LogLikelihood(theta,vals),theta0,A,b,Aeq,beq,lb,ub,nonlcon,options);

end

function LogLik = LogLikelihood(theta,vals)
kon = theta(1);
koff = theta(2);
ksyn = theta(3);
% logprob = zeros(1,length(vals));
[bp,wf] = gaussJacob(50,kon-1,koff-1);
A = 1/beta(kon,koff) * 2^(1-kon-koff);
p = exp( repmat(vals,1,50).*log(repmat(ksyn*(1+bp')/2,length(vals),1)+eps)...
    - repmat(gammaln(vals+1),1,50)...
    - repmat(ksyn*(1+bp')/2,length(vals),1) );
LogLik = -sum(log(A*p*wf + eps));

% for k = 1:length(vals)
%     logprob(k) = A * sum(poisspdf(vals(k),ksyn*(1+bp)/2).*w) + 1e-10;
% end
% LogLik = -sum(log(logprob));

end