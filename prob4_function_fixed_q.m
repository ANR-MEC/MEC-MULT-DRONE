function [B_op, prob4_output] = prob4_function_fixed_q(params, l_op, p_op, m_op, d_op_square, fval4, x0)    
    %aux indicates if it is auxiliary linear problem to find x0

    [K, I, N] = size(p_op);

    options = optimoptions(@linprog, ...
        'Display', 'iter', ...
        'OptimalityTolerance', 10^-6, ...
        'ConstraintTolerance', 10^-6 ...
    );

    prob4 = optimproblem('ObjectiveSense','minimize');

    %variables
    B = optimvar('B',[K,I,N],'Type','continuous','LowerBound',0,'UpperBound',Inf);

    l_op(:,:,N) = 0;
    d_op = sqrt(d_op_square);

    x_bin = d_op <= params.R_c;
    M_k = sum(x_bin,2);
    Delta_t = params.delta ./ M_k;
    Delta_t_mat = repmat(Delta_t, [1,I,1]);

    % constraints

    prob4.Constraints.c4 = sum(x_bin .* B,1) <= params.B_total;


    % path loss
    p_los = 1 ./ (1 + params.a * exp(-params.b * (atan(params.H./d_op) - params.a)));
    L = 20*log10(sqrt(params.H^2+d_op_square)) + params.A*p_los + params.C;
    
    % channel power gain
    h = 10.^(-L/10);
    
    % short blocklength capacity
    gamma = (p_op .* h)/params.sigma2;
    blocklength = repmat(sqrt(1 ./ m_op), [K,1,1]) * params.Q_inv / log(2);
    % blocklength = 1 ./ sqrt(x0.m) * params.Q_inv / log(2);
    % taylor_blocklength = blocklength - 1 ./ (2*x0.m.^(3/2)) * params.Q_inv / log(2) .* (m - x0.m);
    R = B .* ( log(1 + gamma) / log(2) - blocklength );

    prob4.Constraints.c3 = l_op - R .* Delta_t_mat <= 0;


    % objective
    obj = 0; %sum(sum(params.beta * params.delta * (params.theta1*sqrt(veloc_square).^3 + params.theta2 ./ v_bar)));
    prob4.Objective = obj; 

    problem = prob2struct(prob4, x0, ...
        'ConstraintFunctionName', 'generatedConstraintsProb4', ...
        'ObjectiveFunctionName', 'generatedObjectiveProb4');
    problem.options = options;
    [sol,~,exitflag,output] = linprog(problem);

    B_op = reshape(sol(varindex(prob4).B), [K, I, N]);


    sol_struct = struct(B = B_op);
    prob4_output = struct('sol', sol_struct, 'fval', fval4, ...
        'exitflag', exitflag, 'output', output);

    % true_veloc = sqrt(evaluate(veloc_square, struct(q = q_op)));
    % true_fval = sum(sum(params.beta * params.delta * (params.theta1*true_veloc.^3 + params.theta2 ./ true_veloc)));

end