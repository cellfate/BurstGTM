%% master program for MEFs scRNA-seq data
clear;clc;
addpath(genpath(pwd))

tic;
if isempty(gcp('nocreate'))
    parpool(4);
end
toc;
%% load preprocessed data
data_all = csvread('data/MEF_QC_all.csv');
save_folder  = fullfile(pwd,'results/results_MEF');

% check the uninferred gene
gene_number_wait = [];
for gene_number = 1:size(data_all,1)
    filename = sprintf('//result_gene_%d.mat',gene_number);
    if exist([save_folder,filename])==0
        gene_number_wait = [gene_number_wait,gene_number];
    end
end
fprintf('开始运行了\n');

%% main program
parfor infer_index = 1:length(gene_number_wait)
    % Delete 5% of the tail data
    gene = gene_number_wait(infer_index);
    data = data_all(gene,:);
    data = data(data>=0);
    [~,inter] = mink(data,floor(0.95*length(data)));
    data = data(inter);
    data_mean = mean(data);
    data_var = var(data);
    data_noise = data_var/data_mean^2;
    
    % inference
    statis_data = statisData(data);
    rho = @(s) sqrt(sum(log(statis_data./s).^2));
    f = @(k) statisGTM(k,4);
    N = 1000;
    T = 5;
    epsilon = 1;
    % [kon ron koff roff mu delta]
    % OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon-1) * e^(-ron * t) / gamma(kon)
    % ON dwell time distribution f_on(t) = roff^(koff) * t^(koff-1) * e^(-roff * t) / gamma(koff)
    prior = @() [5*rand(),logunif(-1,1),5*rand(),logunif(-1,1),100*rand,1];
    proposal_sigma = 0.2;
    proposal = @(x) lognrnd(x,proposal_sigma);
    proposal_pdf = @(kon_post,kon_prior,ron_post,ron_prior,koff_post,koff_prior,roff_post,roff_prior,mu_post,mu_prior)...
        lognpdf(mu_post,log(mu_prior),proposal_sigma) * lognpdf(kon_post,log(kon_prior),proposal_sigma) *...
        lognpdf(ron_post,log(ron_prior),proposal_sigma) * lognpdf(koff_post,log(koff_prior),proposal_sigma) *...
        lognpdf(roff_post,log(roff_prior),proposal_sigma);
    [result,~] = ABCSMCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,gene);
    
    % save result
    filename = sprintf('//result_gene_%d',gene);
    parsave([save_folder,filename],gene,data,result);
    fprintf('模拟已接受基因%d的结果\n',gene);
end