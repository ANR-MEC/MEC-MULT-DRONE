function [B_op, q_op, t_op, d_op_square, prob4_output, true_fval] = prob4_function(params, l_op, p_op, m_op, x0)    
    %aux indicates if it is auxiliary linear problem to find x0

    [K, I, N] = size(p_op);

    options = optimoptions(@fmincon, 'Algorithm','interior-point', ...
        'Display', 'iter-detailed', 'PlotFcn', {@optimplotfirstorderopt, @optimplotfvalconstr}, ...
        'MaxIterations', 200, ...
        'OptimalityTolerance', 10^-6, ...
        'ConstraintTolerance', 10^-6, ...
        'StepTolerance', 10^-6, ...
        ...'OutputFcn',@outfun, ...
        ...'HessianApproximation', 'finite-difference', ...
        ...'ScaleProblem', true, ...
        ...'SubproblemAlgorithm','cg', ...
        ...'EnableFeasibilityMode',true, ...
        'SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true ...
    );

    prob4 = optimproblem('ObjectiveSense','minimize');

    %variables
    B = optimvar('B',[K,I,N],'Type','continuous','LowerBound',0,'UpperBound',Inf);
    q = optimvar('q',[K,2,N],'Type','continuous','LowerBound',-1,'UpperBound',Inf);
    t = optimvar('t',[K,I,N],'Type','continuous','LowerBound',0,'UpperBound',Inf);
    v_bar = optimvar('v_bar',[K,N-1],'Type','continuous','LowerBound',0,'UpperBound',Inf);

    l_op(:,:,N) = 0;

    % initial v_bar
    x0.v_bar = zeros(K, N-1);
    for n = 1:N-1
        x0.v_bar(:,n) = sqrt(sum((x0.q(:,:,n+1) - x0.q(:,:,n)).^2, 2)) / params.delta;   
    end

    d_square = optimexpr(K,I,N);
    d_op_square = zeros(K,I,N);
    for k = 1:K
        for i = 1:I
            for n = 1:N
                d_square(k,i,n) = sum((q(k,:,n)-params.w(i,:)).^2);
                d_op_square(k,i,n) = sum((x0.q(k,:,n)-params.w(i,:)).^2);
            end
        end
    end
    d_op = sqrt(d_op_square);

    x_bin = d_op <= params.R_c;
    M_k = sum(x_bin,2);
    Delta_t = params.delta ./ M_k;
    Delta_t_mat = repmat(Delta_t, [1,I,1]);

    % constraints

    prob4.Constraints.c4 = sum(x_bin .* B,1) <= params.B_total;

    prob4.Constraints.c12_1 = q(:,:,1) == params.q_I;
    prob4.Constraints.c12_2 = q(:,:,N) == params.q_F;

    veloc_square = optimexpr(K,N-1);
    % for k = 1:K
    for n = 1:N-1
        veloc_square(:,n) = sum((q(:,:,n+1) - q(:,:,n)).^2, 2) / params.delta;   
    end
    % end

    prob4.Constraints.c13 = veloc_square <= params.V_max^2;

    % constraint 21
    % grad_quadratic_B = -2 * (x0.B + x0.t) ./ (x0.B.^3 .* x0.t);
    % grad_quadratic_t = -2 * (x0.B + x0.t) ./ (x0.B .* x0.t.^3);
    % taylor_quadratic = (1./x0.B + 1./x0.t).^2 + grad_quadratic_B .* (B - x0.B) + grad_quadratic_t .* (t - x0.t);
    G = (1./B + 1./t).^2 - 1./x0.B.^2 - 1./x0.t.^2 + 2./x0.B.^3 .* (B - x0.B) + 2./x0.t.^3 .* (t - x0.t);
    prob4.Constraints.c21 = l_op .* G <= 2*Delta_t_mat;

    % constraint 24
    c4 = repmat(sqrt(1./m_op), [K,1,1]) .* params.Q_inv ./ log(2);
    R_bar_op = log2(1 + p_op .* params.beta_0 ./ ((params.H^2 + d_op_square) * params.sigma2)) - c4;
    taylor_t = R_bar_op - (log2(exp(1)) .* p_op .* params.beta_0 .* (d_square - d_op_square)) ./ ...
        ((params.H^2+d_op_square) .* (params.sigma2 * (params.H^2+d_op_square) + p_op * params.beta_0));
    prob4.Constraints.c24 = t <= taylor_t;

    % constraint 25 and 26
    taylor_v = optimexpr(K,N-1);
    taylor_d = optimexpr(K,K,N);
    for k = 1:K
        for n = 1:N

            if n ~= 1
                q0_diff = x0.q(k,:,n) - x0.q(k,:,n-1);
                q_diff = q(k,:,n) - q(k,:,n-1); 
                v_op_square = sum(q0_diff.^2);
                grad = 2*q0_diff;
                taylor_v(k,n-1) = -v_op_square + grad*q_diff';
            end

            for k2 = 1:K
                q0_diff = x0.q(k,:,n) - x0.q(k2,:,n);
                q_diff = q(k,:,n) - q(k2,:,n); 
                dist_square = sum(q0_diff.^2);
                grad = 2*q0_diff;
                if k2 ~= k
                    taylor_d(k,k2,n) = -dist_square + grad*q_diff';
                else
                    taylor_d(k,k2,n) = params.d_min.^2;
                end
            end
        end
    end
    prob4.Constraints.c25 = v_bar.^2 * params.delta^2 <= taylor_v;
    prob4.Constraints.c26 = params.d_min.^2 <= taylor_d;
    
    % objective
    obj = sum(sum(params.beta * params.delta * (params.theta1*sqrt(veloc_square).^3 + params.theta2 ./ v_bar)));
    prob4.Objective = obj; 

    problem = prob2struct(prob4, x0, ...
        'ConstraintFunctionName', 'generatedConstraintsProb4', ...
        'ObjectiveFunctionName', 'generatedObjectiveProb4');
    problem.options = options;
    [sol,fval,exitflag,output] = fmincon(problem);

    % x0.q_var = x0.q(:,:,2:end-1);
    % [sol,fval,exitflag,output,lambda] = solve(prob4, x0, Options=options);

    q_op = reshape(sol(varindex(prob4).q), [K, 2, N]); %cat(3, params.q_I, sol.q_var, params.q_F);
    B_op = reshape(sol(varindex(prob4).B), [K, I, N]);
    t_op = reshape(sol(varindex(prob4).t), [K, I, N]);
    v_bar_op = reshape(sol(varindex(prob4).v_bar), [K, N-1]);

    d_op_square = zeros(K,I,N);
    for k = 1:K
        for i = 1:I
            for n = 1:N
                d_op_square(k,i,n) = sum((q_op(k,:,n)-params.w(i,:)).^2);
            end
        end
    end

    sol_struct = struct(q = q_op, B = B_op, t = t_op, v_bar = v_bar_op);
    prob4_output = struct('sol', sol_struct, 'fval', fval, ...
        'exitflag', exitflag, 'output', output);

    true_veloc = sqrt(evaluate(veloc_square, struct(q = q_op)));
    true_fval = sum(sum(params.beta * params.delta * (params.theta1*true_veloc.^3 + params.theta2 ./ true_veloc)));

end