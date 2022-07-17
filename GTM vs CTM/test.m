clear;clc;
% Simulation algorithm for GTM.
% example1: 3, 0.5, 2, 0.5, 20
% example2: 5,   3, 1, 0.8, 30
param_true.kon = 3;
param_true.ron = 0.5;
param_true.koff = 2;
param_true.roff = 0.5;
param_true.mu = 20;
param_true.delta = 1;
param_true.x0 = [1,0,0];
param_true.tottime = 2000;
burst_true = [1/(param_true.kon/param_true.ron + param_true.koff/param_true.roff), param_true.mu*param_true.koff/param_true.roff];

num_sim = 1000;
results_all = zeros(num_sim,2);
for i = 1:num_sim
    [x,t] = queuingGeneOnOffModel(param_true);
    tq = 1000:0.1:param_true.tottime;
    xq = interp1(t,x(:,3),tq,'previous');
    
    % figure
    % histogram(xq,'normalization','pdf')
    
    data = xq;
    [theta_est,theta0] = burstInference(data);
    
    koff_est = theta_est(1);
    kon_est = theta_est(2);
    ksyn_est = theta_est(3);
    burst_est = [1/(1/koff_est + 1/kon_est), ksyn_est/kon_est];
    results_all(i,:) = burst_est;
end

log_burst_true = log10(burst_true);

[f, logbrust] = ksdensity([log10(results_all(:,1)),log10(results_all(:,2))]);
bf_center = logbrust(find(f == max(f)),1);
bs_center = logbrust(find(f == max(f)),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the inferred results
max_expr = max(data+1);
[bp,w] = gaussJacob(50,koff_est-1,kon_est-1);
prob_est = zeros(1,max_expr+1);
for m = 0:max_expr
    prob_est(m+1) = 1/beta(koff_est,kon_est) * 2^(1-koff_est-kon_est) * sum(poisspdf(m,ksyn_est*(1+bp)/2).*w);
end

Color = {'#EE2201'};
FaceColor = {'#00837E'};

f1 = figure;
histogram(data,0-0.5:max_expr+0.5,'Normalization','pdf','EdgeColor','none','FaceColor',FaceColor{1})
xlim([-0.5 max_expr])
hold on
plot(0:max_expr,prob_est,'r','LineWidth',1,'Color',Color{1},'LineWidth',1.5)
hold off
legend('GTM','CTM')
xlabel('mRNA')
ylabel('Prob.')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);

set(f1,'position',[400 400 200 120]);

f2 = figure;
tau_off_true = gamrnd(param_true.kon,1/param_true.ron,1,2000);
tau_off_est = exprnd(1./koff_est,1,2000);
y_tau_off_est = exppdf(linspace(0,max(tau_off_est),100),1./koff_est);

subplot(1,2,1)
histogram(tau_off_true,'Normalization','pdf','EdgeColor','none','FaceColor',FaceColor{1})
hold on
plot(linspace(0,max(tau_off_est),100),y_tau_off_est,'Color',Color{1},'LineWidth',1.5)
hold off
% xlim([0,max(max(tau_off_true),max(tau_off_est))])
xlim([0,5])
xlabel('OFF state dwell time')
ylabel('Prob.')
axis square
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);

tau_on_true = gamrnd(param_true.koff,1/param_true.roff,1,2000);
tau_on_est = exprnd(1./kon_est,1,2000);
y_tau_on_est = exppdf(linspace(0,max(tau_on_est),100),1./kon_est);

subplot(1,2,2)
histogram(tau_on_true,'Normalization','pdf','EdgeColor','none','FaceColor',FaceColor{1})
hold on
plot(linspace(0,max(tau_on_est),100),y_tau_on_est,'Color',Color{1},'LineWidth',1.5)
hold off
% xlim([0,max(max(tau_on_true),max(tau_on_est))])
xlim([0,5])
xlabel('ON state dwell time')
ylabel('Prob.')
axis square
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);

set(f2,'position',[400 400 264 132]);

f3 = figure;
scatplot(log10(results_all(:,1)),log10(results_all(:,2)),'circles')
hold on
plot(log10(burst_true(1)),log10(burst_true(2)),'*','MarkerEdgeColor','red','Markersize',10)
plot(bf_center,bs_center,'x','MarkerEdgeColor','red','Markersize',10)
xlabel('log(bf)')
ylabel('log(bs)')
% xlim([-0.55,-0.05])
% ylim([1.15,1.65])
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
axis square
axis on
box on
set(f3,'position',[400 400 178.5 178.5]);

save GTMvsCTM_1.mat