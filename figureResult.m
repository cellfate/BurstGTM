function  [] = figureResult(data,result,param_true)
kon = cellfun(@(c) c.kon,result);
ron = cellfun(@(c) c.ron,result);
koff = cellfun(@(c) c.koff,result);
roff = cellfun(@(c) c.roff,result);
mu = cellfun(@(c) c.mu,result);
dist = cellfun(@(c) c.dist,result);

tau_off = kon./ron;
tau_on = koff./roff;
bf = 1./(tau_on + tau_off);
bs = mu .* tau_on;

param_true.tau_off = param_true.kon./param_true.ron;
param_true.tau_on = param_true.koff./param_true.roff;
param_true.bf = 1./(param_true.tau_on + param_true.tau_off);
param_true.bs = param_true.mu .* param_true.tau_on;


[f, tau] = ksdensity([log10(tau_off), log10(tau_on)]);

tau_offCenter = tau(find(f == max(f)), 1);
tau_onCenter = tau(find(f == max(f)), 2);
Idx = knnsearch(log10([tau_off,tau_on]),[tau_offCenter,tau_onCenter]);

[~,min_ind] = min(dist);
param_est = result{min_ind,1};
param_est.x0 = [1,0,0];
param_est.tottime = 2000;
[x_est,t_est] = simulGTM(param_est);
tq_est = 1000:0.1:param_est.tottime;
xq_est = interp1(t_est,x_est(:,3),tq_est,'previous');
data_est = xq_est;
pts = linspace(1,max(data)+10,100);
[f_true,xi_true] = ksdensity(data + 1,pts,'Support','positive','Bandwidth',0.3); 

%%%%%%
Color = {'#EE2201'};
FaceColor = {'#00837E','#4DBBD4'};

f1 = figure;
histogram(data_est,'normalization','pdf','BinEdges',0:max(data_est),'EdgeColor','none','FaceColor',FaceColor{1})
hold on
plot(xi_true-0.5,f_true,'r','LineWidth',1,'Color',Color{1},'LineWidth',1.5)
xlim([0 max(data_est)+5])
xlabel('mRNA')
ylabel('PDF')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
set(f1,'position',[400 400 200 120]);

f2 = figure;
h = histogram(mu,'normalization','pdf','EdgeColor','none','FaceColor',FaceColor{2});
[~, max_ind] = max(h.Values);
muCenter = mean(h.BinEdges(max_ind:max_ind + 1));
hold on;
plot([param_true.mu param_true.mu], get(gca, 'YLim'), '-r', 'LineWidth', 1) 
plot([muCenter muCenter], get(gca, 'YLim'), '--k', 'LineWidth', 1) 
xlabel('\mu')
ylabel('PDF')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
set(f2,'position',[400 400 120 120]);

f3 = figure;
subplot(1,2,1)
scatplot(log10(tau_off),log10(tau_on),'circles')
hold on
plot([log10(param_true.tau_off) log10(param_true.tau_off)], get(gca, 'YLim'), '-r', 'LineWidth', 1) 
plot([tau_offCenter tau_offCenter], get(gca, 'YLim'), '--k', 'LineWidth', 1) 
plot(get(gca, 'XLim'), [log10(param_true.tau_on) log10(param_true.tau_on)], '-r', 'LineWidth', 1) 
plot(get(gca, 'XLim'), [tau_onCenter tau_onCenter], '--k', 'LineWidth', 1)
xlim(get(gca, 'XLim'))
ylim(get(gca, 'YLim'))
xlabel('log(<\tau_{off}>)')
ylabel('log(<\tau_{on}>)')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
axis square
axis on
box on

subplot(1,2,2)
scatplot(log10(bf),log10(bs),'circles')
hold on
plot([log10(param_true.bf) log10(param_true.bf)], get(gca, 'YLim'), '-r', 'LineWidth', 1) 
[f, log_burst] = ksdensity([log10(bf), log10(bs)]);
log_bfCenter = log_burst(find(f == max(f)), 1);
log_bsCenter = log_burst(find(f == max(f)), 2);
plot([log_bfCenter log_bfCenter], get(gca, 'YLim'), '--k', 'LineWidth', 1) 
plot(get(gca, 'XLim'), [log10(param_true.bs) log10(param_true.bs)], '-r', 'LineWidth', 1) 
plot(get(gca, 'XLim'), [log_bsCenter log_bsCenter], '--k', 'LineWidth', 1)

xlim(get(gca, 'XLim'))
ylim(get(gca, 'YLim'))
xlabel('log(bf)')
ylabel('log(bs)')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
axis square
box on
set(f3,'position',[400 400 400 200]);
end
