t = load('table_fig_M_eps.mat').t;

fig = figure();
fig.Units = 'pixels';
fig.Position = [1 1 500 500];

M_list = sort(unique(t.M), 'descend');
epsilon = uniquetol(t.epsilon, 10^-10);

markers = {'o-', '^-', 'd-'};
for m = 1:length(M_list)
    M = M_list(m);
    epsilon = t(t.M == M,:).epsilon;
    fval = t(t.M == M,:).fval;
    plot(epsilon, fval, markers{m}, 'DisplayName', ['M = ', num2str(M)]); hold on
end

ax = gca;
ax.XScale = 'log';
legend show
xlabel('Decoding error')
ylabel('Weighted sum energy consumption (J)')
grid on

saveas(fig, 'figure_diffM_diffEps.eps', 'epsc')
