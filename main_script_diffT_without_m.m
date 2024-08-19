clear variables

%% PARAMETERS

%Number of UAVs
K = 4;
%Number of users
I = 16;

%Location of users
x_user = linspace(0,10,4);
y_user = linspace(0,10,4);
params.w = table2array(combinations(x_user, y_user));

%UAVs height: H=10m
params.H = 10;
%Total bandwidth 30 MHz 
params.B_total = 30*10^6;
% Vmax 10 m/s
params.V_max = 10;
%Noise power -50 dBm
sigma2_dbm = -70;
params.sigma2 = 10^(sigma2_dbm/10) * 10^-3; % * params.B_total;
% dmin 1 m
params.d_min = 1;
%Number of time slot 20 
N = 20;
%θ1, θ2 0.00614, 15.976 
params.theta1 = 0.00614; 
params.theta2 = 15.976;
% α, β 0.2, 0.4
params.alpha = 0.2;
params.beta = 0.4;
%ηLoS, ηNLoS 20, 1 
eta_los = 20; 
eta_nlos = 1;
% gk 8GHz
params.g_k = 8*10^9;
%a, b 9.61, 0.16 
params.a = 9.61;
params.b = 0.16;
% fi 1GHz
params.f_i = 10^9;

% task size in Mbits
params.L_i = 4 * 10^6 * ones(I,1);
params.C_i = 40 * ones(I,1);

% capacitance coefficients
params.kappa_user = 10^-26;
params.kappa_uav = 10^-26;

% % time in seconds
% T = linspace(2, 2.5, 6);
% params.delta = T / N;

% path loss
c = 3*10^8;
f_c = 2.5*10^9; 
params.A = eta_los - eta_nlos;
params.C = 20 * log10(4*pi*f_c/c) + eta_nlos;
params.beta_0 = 10^(-params.C/10);

% total blocklength 
params.M_total = 1000;

% error
epsilon = 10^-7;
params.Q_inv = sqrt(2) * erfcinv(2*epsilon);

% coverage area
phi_max = pi/4;
params.R_c = params.H*tan(phi_max);


%% Different tasks

% time in seconds
T = [1.5,1.9];

for t = 1:length(T)

    params.delta = T(t) / N;
    
    %% Initial values
    
    %Initial and final location of UAVs: !!TODO
    params.q_I = [0,5; 5,0; 10,5; 5,10];
    params.q_F = [5,10; 0,5; 5,0; 10,5];
    
    %Initial location
    q_op = zeros(K,2,N);
    d_op_square = zeros(K,I,N);
    for k = 1:K
        q_op(k,:,:) = [linspace(params.q_I(k,1),params.q_F(k,1),N)', linspace(params.q_I(k,2),params.q_F(k,2),N)']';
        for i = 1:I
            for n = 1:N
                d_op_square(k,i,n) = sum((q_op(k,:,n)-params.w(i,:)).^2);
            end
        end
    end
    
    %Initial bandwidth distribution
    B_op = 5*10^5 .* ones(K,I,N);

    % Initial power
    p_op = 1 * ones(K,I,N);
    
    % Initial blocklength
    m_op = params.M_total / I * ones(1,I,N);
    
    %Inital CPU frequencies for users
    f_op = 10^6 * ones(I,N);

    %Initial amount of offloading task bit
    l_op = repmat((params.L_i - params.delta .* sum(f_op,2) ./ params.C_i)' ./ (K*(N-1)), [K, 1, N-1]);
    
    %Inital CPU frequencies for UAVs
    g_op = l_op .* repmat(params.C_i', [K, 1, N-1]) / params.delta; %params.g_k/I * ones(K,I,N);
    g_op(:,:,N) = g_op(:,:,N-1);


    % find feasible solution for the linear problem
    x0 = struct();
    x0.f = f_op;
    x0.g = g_op;
    x0.l = l_op;
    x0.p = p_op;
    % x0.m = m_op;
    lin = true;
    obj = false;
    [f_op, g_op, l_op, p_op, ~, ~] = prob2_function_without_m(params, B_op, d_op_square, lin, obj, x0);

    % find feasible solution for the nonlinear problem
    x0 = struct();
    x0.f = f_op;
    x0.g = g_op;
    x0.l = l_op;
    x0.p = p_op;
    % x0.m = m_op;
    lin = false;
    obj = false;
    [f_op, g_op, l_op, p_op, ~, ~] = prob2_function_without_m(params, B_op, d_op_square, lin, obj, x0);


    %% Problem 

    %tolerance for main algorithm
    max_tol = 10^-5;
    max_iter = 15;
    max_iter_alg1 = 10;

    tol = 1;
    iteration = 1; 
    prob2_output_list = {};
    prob4_output_list = {};
    while tol > max_tol && iteration <= max_iter

        fprintf('Iteration # %i \n',iteration )
        fprintf('############################################## \n')

        %% P2

        x0 = struct();
        x0.f = f_op;
        x0.g = g_op;
        x0.p = p_op;
        % x0.m = m_op;
        x0.l = l_op;

        lin = false;
        obj = true;
        [f_op, g_op, l_op, p_op, prob2_output, fval] = prob2_function_without_m(params, B_op, d_op_square, lin, obj, x0);
        fval2 = fval;
        prob2_output_list{iteration} = prob2_output;


        %% P4

        % initialize variable t <= R_bar
        % c4 = repmat(sqrt(1./m_op), [K,1,1]) .* params.Q_inv ./ log(2);
        R_bar_op = log2(1 + p_op .* params.beta_0 ./ ((params.H^2 + d_op_square) * params.sigma2)); % - c4);
        t_op = R_bar_op; 
        
        iteration_alg1 = 1;
        tol_alg1 = 1;
        while tol_alg1 > max_tol && iteration_alg1 <= max_iter_alg1

            x0 = struct();
            x0.B = B_op;
            x0.q = q_op; %(:,:,2:end-1);
            x0.t = t_op;
            % x0.v_bar = v_bar_op;

            [B_op, q_op, t_op, d_op_square, prob4_output, fval] = prob4_function_without_m(params, l_op, p_op, x0);     
            fval4 = fval;
            prob4_output_list{iteration, iteration_alg1} = prob4_output;

            obj_iter_alg1(iteration_alg1) = fval4;
            if iteration_alg1 ~= 1
                tol_alg1 = abs(obj_iter_alg1(iteration_alg1-1) - obj_iter_alg1(iteration_alg1));
            end

            iteration_alg1 = iteration_alg1 + 1;


        end  

        obj_final = fval2 + fval4; 
        obj_iter(iteration) = obj_final;
        if iteration ~= 1
            tol = abs(obj_iter(iteration-1) - obj_iter(iteration));
        end

        iteration = iteration + 1;

    end

    new_date = replace(char(datetime), ':', '-');
    save(['without_m_prob2_output_', new_date, '.mat'], 'prob2_output_list')
    save(['without_m_prob4_output_', new_date, '.mat'], 'prob4_output_list')
    save(['without_m_params_', new_date, '.mat'], 'params')

    % figure
    % plot(params.w(:,1), params.w(:,2), '*'); hold on
    % % xlim([0,10])
    % % ylim([0,10])
    % for k = 1:K
    %     for n = 1:N
    %         coords(n,:) = [q_op(k,1,n),q_op(k,2,n)];
    %     end
    %     plot(coords(:,1), coords(:,2),'.-'); hold on
    % end

end