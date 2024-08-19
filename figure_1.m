
% Load the table
table_proposed = load('table_fig_trajectory.mat').t;

% Get all unique row names
selected_rows = table_proposed.name;

% Initialize the figure
figure;

% Define colors and markers for each row (adjust as needed)
col_uavs = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0 0.4470 0.7410], [0.4660 0.6740 0.1880]};
marker_styles = {'d-', 'x-', 'o-', 's-'};  % Different markers for different rows


% legend 
L(1) = plot(table_proposed.w{1}(:,1), table_proposed.w{1}(:,2), 'k*'); hold on
L(2) = plot(nan, nan, 'd-k'); hold on
L(3) = plot(nan, nan, 'x-k'); hold on


% Loop through all selected rows and plot their trajectories
for i = 1:length(selected_rows)
    % Find the row corresponding to the current selected_row
    row_idx = strcmp(table_proposed.name, selected_rows{i});
    
    % Check if the row exists
    if sum(row_idx) == 1
        % Extract the 'w' values and 'sol' struct for the selected row
        w_values = table_proposed.w{row_idx};
        sol_values = table_proposed.sol(row_idx).q;

        % Set up parameters
        K = size(sol_values, 1); 
        N = size(sol_values, 3); 

        % Plot the initial points
        plot(w_values(:,1), w_values(:,2), 'k*'); hold on
        a = [1:size(w_values, 1)]'; b = num2str(a); c = cellstr(b);
        dx = 0.1; dy = 0; % displacement so the text does not overlay the data points
        text(w_values(:,1) + dx, w_values(:,2) + dy, c);

        % Plot the trajectory
        for k = 1:K
            coords = zeros(N, 2); % Preallocate coords array
            for n = 1:N
                coords(n,:) = [sol_values(k,1,n), sol_values(k,2,n)];
            end
            plot(coords(:,1), coords(:,2), marker_styles{i}, 'Color', col_uavs{k}); hold on
        end
    else
        disp(['The selected row "', selected_rows{i}, '" does not exist in the table.']);
    end
end

% Set axis limits
xlim([-1 11])
ylim([-1 11])

% Add grid and labels
grid on
xlabel('x(m)')
ylabel('y(m)')

% Add title indicating the plot includes all rows
title('Trajectory Plot: All Configurations');

% Add legend for all rows
leg = {'IoT devices','Non-uniform', 'Uniform'};
legend(leg, 'Location', 'Best');

% Save the figure
saveas(gcf,'figure_all_rows.eps', 'epsc')