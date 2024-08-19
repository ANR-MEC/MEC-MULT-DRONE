function [f_op, g_op, l_op, p_op, prob2_output, true_fval] = prob2_function_without_m(params, B_op, d_op_square, lin, obj, x0)    
    %aux indicates if it is auxiliary linear problem to find x0
    
    [K, I, N] = size(B_op);

    options = optimoptions(@fmincon, 'Algorithm','interior-point', ...
        'Display', 'iter-detailed', 'PlotFcn', {@optimplotfirstorderopt, @optimplotfvalconstr}, ...
        'MaxIterations', 200, ...
        'OptimalityTolerance', 10^-6, ...
        'ConstraintTolerance', 10^-6, ...
        ...'OutputFcn',@outfun, ...
        ...'StepTolerance', 10^-4, ...
        ...'HessianApproximation', 'finite-difference', ...
        'ScaleProblem', true, ...
        ...'SubproblemAlgorithm','cg', ...
        ...'EnableFeasibilityMode',true, ...
        'SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true ...
    );

    options_lin = optimoptions(@linprog, 'MaxIterations', 10000, ...
        'Display', 'iter', 'ConstraintTolerance', 10^-6);

    prob2 = optimproblem('ObjectiveSense','minimize');

    %variables
    l = optimvar('l',[K,I,N-1],'Type','continuous','LowerBound',0,'UpperBound',inf);
    % m = optimvar('m',[1,I,N],'Type','continuous','LowerBound',realmin,'UpperBound',inf);
    p = optimvar('p',[K,I,N],'Type','continuous','LowerBound',0,'UpperBound',inf);
    f = optimvar('f',[I,N],'Type','continuous','LowerBound',0,'UpperBound',inf);
    g = optimvar('g',[K,I,N],'Type','continuous','LowerBound',0,'UpperBound',inf);
    
    % distance
    d_op = sqrt(d_op_square);

    % path loss
    p_los = 1 ./ (1 + params.a * exp(-params.b * (atan(params.H./d_op) - params.a)));
    L = 20*log10(sqrt(params.H^2+d_op_square)) + params.A*p_los + params.C;
    
    % channel power gain
    h = 10.^(-L/10);
    
    % short blocklength capacity
    gamma = (p .* h)/params.sigma2;
    % blocklength = repmat(sqrt(1 ./ m_op), [K,1,1]) * params.Q_inv / log(2);
    % blocklength = 1 ./ sqrt(x0.m) * params.Q_inv / log(2);
    % taylor_blocklength = blocklength - 1 ./ (2*x0.m.^(3/2)) * params.Q_inv / log(2) .* (m - x0.m);
    R = B_op .* ( log(1 + gamma) / log(2) );
    

    % binary indicating if user i is inside coverage area of uav k
    x_bin = d_op <= params.R_c;
    M_k = sum(x_bin,2);
    Delta_t = params.delta ./ M_k;
    Delta_t_mat = repmat(Delta_t, [1,I,1]);

    % energies
    E_off = sum(p .* x_bin .* Delta_t_mat, 'all');
    E_local = sum(params.delta * params.kappa_user * f.^3, 'all');
    % E_local_taylor = E_local + sum(params.delta * params.kappa_user * 3*x0.f.^2 .* (f - x0.f), 'all');
    E_exe = sum(params.delta * params.kappa_uav * g.^3, 'all');
    % E_exe_taylor = E_exe + sum(params.delta * params.kappa_uav * 3*x0.g.^2 .* (g - x0.g), 'all');
    % E_exe = sum(Delta_t_mat .* params.kappa_uav .* g.^3, 'all');
    E_user = E_local + E_off;
    

    % constraints

    if ~lin
        prob2.Constraints.c3 = l - R(:,:,1:N-1) .* Delta_t_mat(:,:,1:N-1) <= 0;
        prob2.Constraints.c3_2 = 0 <= R(:,:,N) .* Delta_t_mat(:,:,N);

        % prob2.Constraints.c6 = sum(m, 2) <= params.M_total;

    else

        % m_op_aux = repmat(m_op, [K,1,1]);
        % m_op = x0.m;
        % c4 = sqrt(1./m_op_aux) .* params.Q_inv ./ log(2);

        prob2.Constraints.c_aux1 = l ./ (B_op(:,:,1:N-1) .* Delta_t_mat(:,:,1:N-1)) <= 10;

    %     prob2.Constraints.c_aux1 = diff(f) == 0;
    %     prob2.Constraints.c_aux2 = diff(f,1,2) == 0;
    % 
    %     % prob2.Constraints.c_aux3 = diff(l) == 0;
    %     % prob2.Constraints.c_aux4 = diff(l,1,2) == 0;
    %     % prob2.Constraints.c_aux5 = diff(l,1,3) == 0;
    % 
    %     prob2.Constraints.c_aux6 = l >= 1000;
    %     prob2.Constraints.c_aux7 = l <= 2*10^5;
    %     % prob2.Constraints.c_aux8 = g >= 1000;
    % 
        prob2.Constraints.c_aux_Eoff = E_off <= 500;

    end
    
    prob2.Constraints.c7 = params.delta ./ params.C_i .* sum(f, 2) + sum(sum(l, 1), 3)' == params.L_i;

    prob2.Constraints.c8 = optimconstr(K,N);
    for k = 1:K
        for n = 2:N
            prob2.Constraints.c8(k,n) = params.delta * sum(sum(g(k,:,1:n) ./ repmat(params.C_i', [1,1,n]), 2),3) <= sum(sum(l(k,:,1:n-1), 2),3);
        end
    end

    prob2.Constraints.c9 = params.delta * sum(sum(g ./ repmat(params.C_i', [K,1,N]), 2),3) == sum(sum(l(:,:,1:N-1), 2),3);

    prob2.Constraints.c18b1 = f <= params.f_i;
    prob2.Constraints.c18b2 = sum(g,2) <= params.g_k;

    % prob2.Constraints.c18c = l(:,:,N) == 0;

    % objective
    if lin
        prob2.Objective = 0; %-sum(f, 'all') - sum(g, 'all') - sum(m, 'all') + sum(l, 'all');
        problem = prob2struct(prob2, x0, ...
            'ConstraintFunctionName', 'generatedConstraintsProb2', ...
            'ObjectiveFunctionName', 'generatedObjectiveProb2');
        problem.options = options_lin;
        [sol,fval,exitflag,output] = linprog(problem);
        % [sol,fval,exitflag,output,lambda] = solve(prob2, x0, Options=options_lin);
    else
        if obj
            obj_func = E_user + params.alpha .* E_exe;
        else
            obj_func = 0;
        end
        
        % obj = E_off;
        prob2.Objective = obj_func;
        problem = prob2struct(prob2, x0, ...
            'ConstraintFunctionName', 'generatedConstraintsProb2', ...
            'ObjectiveFunctionName', 'generatedObjectiveProb2');
        problem.options = options;

        % H = @(x,lam)(hessianProb2(x, params, lam, Delta_t_mat, B_op, h));
        % problem.options.HessianFcn = H;

        [sol,fval,exitflag,output] = fmincon(problem);
        % [sol,fval,exitflag,output,lambda] = solve(prob2, x0, Options=options);        
    end

    f_op = reshape(sol(varindex(prob2).f), [I, N]);
    g_op = reshape(sol(varindex(prob2).g), [K, I, N]);
    l_op = reshape(sol(varindex(prob2).l), [K, I, N-1]);
    if lin
        % m_op_aux = repmat(m_op, [K,1,1]);
        % m_op = x0.m;
        % B_op = reshape(sol(varindex(prob2).B), [K, I, N]);
        % p_op = x0.p;
        % c4 = sqrt(1./m_op_aux) .* params.Q_inv ./ log(2);
        p_op(:,:,1:N-1) = 2*params.sigma2 ./ h(:,:,1:N-1) .* (2.^(l_op ./ (B_op(:,:,1:N-1) .* Delta_t_mat(:,:,1:N-1))) - 1);
        p_op(:,:,N) = 0; %2*params.sigma2 ./ h(:,:,N) .* (2.^c4(:,:,N) - 1);
    else
        % m_op = reshape(sol(varindex(prob2).m), [1, I, N]);
        p_op = reshape(sol(varindex(prob2).p), [K, I, N]);
    end

    sol_struct = struct(f = f_op, g = g_op, l = l_op, p = p_op);

    prob2_output = struct('prob', prob2, 'sol', sol_struct, 'fval', fval, ...
        'exitflag', exitflag, 'output', output);

    % veloc = zeros(K,N-1);
    % for n = 1:N-1
    %     veloc(:,n) = sqrt(sum((q_op(:,:,n+1) - q_op(:,:,n)).^2, 2)) / params.delta;   
    % end
    true_fval = fval; % + sum(params.beta * params.delta * (params.theta1*veloc.^3 + params.theta2 ./ veloc), 'all');