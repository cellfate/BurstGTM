%% collect result
clear;clc;
genefiles = dir('results\results_MEF\*.mat');

results = zeros(length(genefiles),2);
for i = 1:length(genefiles)
    s = genefiles(i).name;
    s = str2num(s(isstrprop(s,'digit')));
    filename = sprintf('results\\results_MEF\\%s',genefiles(i).name);
    load(filename);
    
    kon = cellfun(@(c) c.kon,result(:,end));
    ron = cellfun(@(c) c.ron,result(:,end));
    koff = cellfun(@(c) c.koff,result(:,end));
    roff = cellfun(@(c) c.roff,result(:,end));
    mu = cellfun(@(c) c.mu,result(:,end));
    dist = cellfun(@(c) c.dist,result(:,end));
    
    tau_off = kon./ron;
    tau_on = koff./roff;
    bf = 1./(tau_off + tau_on);
%     bf = 1./(tau_off);
    bs = mu.*tau_on;
    
    [f, burst] = ksdensity([bf, bs]);
    bf_center = burst(find(f == max(f)),1);
    bs_center = burst(find(f == max(f)),2);
    
    results(s,:) = [bf_center,bs_center];
    disp(i)
end
csvwrite('results\results_MEF.csv',results);