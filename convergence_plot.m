t = load('table_fig_convergence.mat').t;

markers = {'-o', '-^', '-d'};
leg = {};
fig = figure; 
for i = 1:size(t,2)
    plot(0:30, t(i,:).fval{1}, markers{i}, 'DisplayName', [num2str(t.sumL(i)/10^6), 'Mbits']); hold on
end 

legend show
xlabel('Iteration')
ylabel('Weighted sum energy consumption (J)')
grid on

saveas(fig, 'figure_convergence.eps', 'epsc')
