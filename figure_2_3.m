t_fig2 = load('table_fig_sumL.mat').t_fig2;
t_fig3 = load('table_fig_T.mat').t_fig3;

% FIG 2 - sumL

fig2 = figure();
fig2.Units = 'pixels';
fig2.Position = [1 1 500 500];


% PROPOSED
x = t_fig2(ismember(t_fig2.method, 'Proposed algorithm'),:).L;
fval_proposed = t_fig2(ismember(t_fig2.method, 'Proposed algorithm'),:).fval;
plot(x/10^6, fval_proposed, 'o-', 'DisplayName', 'Proposed algorithm'); hold on

% FIXED B
x = t_fig2(ismember(t_fig2.method, 'Fixed bandwidth'),:).L;
fval_B = t_fig2(ismember(t_fig2.method, 'Fixed bandwidth'),:).fval;
plot(x/10^6, fval_B, '^-', 'DisplayName', 'Fixed bandwidth'); hold on

% FIXED q
x = t_fig2(ismember(t_fig2.method, 'Fixed trajectory'),:).L;
fval_q = t_fig2(ismember(t_fig2.method, 'Fixed trajectory'),:).fval;
plot(x/10^6, fval_q, 'd-', 'DisplayName', 'Fixed trajectory'); hold on

% FIXED m
x = t_fig2(ismember(t_fig2.method, 'Fixed blocklength'),:).L;
fval_m = t_fig2(ismember(t_fig2.method, 'Fixed blocklength'),:).fval;
plot(x/10^6, fval_m, 'x-', 'DisplayName', 'Fixed blocklength'); hold on


legend show
xlabel('The total amount of task input-bits (Mbits)')
ylabel('Weighted sum energy consumption (J)')
grid on

saveas(fig2, 'figure_2.eps', 'epsc')

disp('Figure 2')
disp(['Fixed B uses on average ', num2str(mean(fval_B ./ fval_proposed * 100) - 100), '% more energy than proposed'])
disp(['Fixed q uses on average ', num2str(mean(fval_q ./ fval_proposed * 100) - 100), '% more energy than proposed'])
disp(['Fixed m uses on average ', num2str(mean(fval_m ./ fval_proposed * 100) - 100), '% more energy than proposed'])


% FIG 3 - T

fig3 = figure();
fig3.Units = 'pixels';
fig3.Position = [1 1 500 500];

% PROPOSED
x = t_fig3(ismember(t_fig3.method, 'Proposed algorithm'),:).T;
fval_proposed = t_fig3(ismember(t_fig3.method, 'Proposed algorithm'),:).fval;
plot(x, fval_proposed, 'o-', 'DisplayName', 'Proposed algorithm'); hold on

% FIXED B
x = t_fig3(ismember(t_fig3.method, 'Fixed bandwidth'),:).T;
fval_B = t_fig3(ismember(t_fig3.method, 'Fixed bandwidth'),:).fval;
plot(x, fval_B, '^-', 'DisplayName', 'Fixed bandwidth'); hold on

% FIXED q
x = t_fig3(ismember(t_fig3.method, 'Fixed trajectory'),:).T;
fval_q = t_fig3(ismember(t_fig3.method, 'Fixed trajectory'),:).fval;
plot(x, fval_q, 'd-', 'DisplayName', 'Fixed trajectory'); hold on

% FIXED m
x = t_fig3(ismember(t_fig3.method, 'Fixed blocklength'),:).T;
fval_m = t_fig3(ismember(t_fig3.method, 'Fixed blocklength'),:).fval;
plot(x, fval_m, 'x-', 'DisplayName', 'Fixed blocklength'); hold on

legend show
xlabel('Mission period T (sec)')
ylabel('Weighted sum energy consumption (J)')
ylim([25 60])
grid on

saveas(fig3, 'figure_3.eps', 'epsc')

disp('Figure 3')
disp(['Fixed B uses on average ', num2str(mean(fval_B ./ fval_proposed * 100) - 100), '% more energy than proposed'])
disp(['Fixed q uses on average ', num2str(mean(fval_q ./ fval_proposed * 100) - 100), '% more energy than proposed'])
disp(['Fixed m uses on average ', num2str(mean(fval_m ./ fval_proposed * 100) - 100), '% more energy than proposed'])

