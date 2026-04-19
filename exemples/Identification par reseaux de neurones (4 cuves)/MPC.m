classdef MPC < handle
    properties
        phi;             % system evolution function
        N;               % Prediction horizon
        Q;               % state stage cost weight matrix x_i^T Q x
        R;               % input stage cost weight matrix u_i^T R u
        P;               % terminal cost weight matrix u_i^T P u
        xmin;            % vector of min value for each state
        xmax;            % vector of max value for each state
        umin;            % vector of min value for each input
        umax;            % vector of max value for each input
        nx;              % size of state
        nu;              % size of input
        opti = [];       % Casadi opti object
        opti_x = [];     % (casadi.opti.variable) symbol of the states decision variable
        opti_x0 = [];    % (casadi.opti.parameter) symbol of the initial state
        opti_ref = [];   % (casadi.opti.parameter) symbol of reference target/trajectory
        opti_u = [];     % (casadi.opti.variable) symbol of input decision variable
        prev_sol = [];
        use_warm_start=true; % whether to use warm start. This is useful to toggle it internally (RTI)
    end

    methods
        function obj = MPC(phi, N, Q, R, P, xmin, xmax, umin, umax)
        %MPC Class constructor for setting up an MPC controller
        %
        % Usage:
        %
        % mpc = MPC(PHI, N, Q, R, P, XMIN, XMAX, UMIN, UMAX) creates and
        % returns an MPC object with state evolution function :
        %                    x_k+1 = PHI(x_k, u_k)
        %
        % over a control horizon of length N and quadratic cost function:
        %            V_k = (x_k-xref_k)'*Q*(x_k-xref_k) + u_k'*R*u_k.
        %
        % P is the weight on the terminal cost (x_N-xref_N)'P(x_N-xref_N).
        % XMIN, XMAX, UMIN, UMAX are respectively the constraints on the lower
        % state bound, upper state bound, lower input bound and upper input bound.

            obj = obj@handle();
            obj.phi = phi;
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.P = P;
            obj.xmin = xmin;
            obj.xmax = xmax;
            obj.umin = umin;
            obj.umax = umax;
            obj.nx = length(xmin);
            obj.nu = length(umin);
            [obj.opti, obj.opti_x, obj.opti_x0, obj.opti_ref, obj.opti_u] = obj.construct_problem();
        end

        function [opti, x, x0, ref, u] = construct_problem(obj)
        %CONSTRUCT_PROBLEM Create the MPC optimization problem
        % possibilite de generalisation en travaillant à partir de yref
        % en lieu et place de xref
        %
        % Usage :
        % [uopt, U, X] = mpc.solve(X0, REF) generates the first
        % control action UOPT, sequence of actions U and sequence of
        % predictions X over the control horizon that are produced when
        % driving the system from current state X0 to reference REF.
        %
        % You would typically want to run this function in a loop,
        % where you could compute the control action and apply it to
        % the system in the following fashion:
        %
        % for i=1:1:Tf
        %     [uopt, ~, ~] = mpc.solve(x0, ref);
        %      x = simulate_system(x, uopt);
        % end

            opti = casadi.Opti();
            opti.solver('ipopt',...
                        struct('print_time', false), ...
                        struct('print_level', 0));

            x = opti.variable(obj.nx, obj.N+1);
            u = opti.variable(obj.nu, obj.N);
            x0 = opti.parameter(obj.nx, 1); opti.subject_to(x(:, 1) == x0); % initial state
                                                                            %u0 = opti.parameter(obj.nu, 1);
            ref = opti.parameter(obj.nx, obj.N+1);

            cost = (x(:, end)-ref(:, end))'*obj.P*(x(:, end)-ref(:, end));  % Terminal Cost
                                                                            %cost = cost + (u(:, 1) - u0)'*obj.R*(u(:, 1)-u0);
            for i=1:obj.N
                if i~=1
                    cost = cost + (u(:, i)-u(:, i-1))'*obj.R*(u(:, i)-u(:, i-1));
                    %cost = cost + u(:, i)'*obj.R*u(:, i);
                end
                cost = cost + (x(:,i)-ref(:,i))'*obj.Q*(x(:,i)-ref(:, i));  % Stage cost
                opti.subject_to(x(:,i+1) == obj.phi(x(:,i), u(:,i)));  % State Propagation Constraints
                for j=1:obj.nx
                    opti.subject_to(obj.xmin(j) <= x(j, i+1) <= obj.xmax(j));  % state constraints
                end
                for j=1:obj.nu
                    opti.subject_to(obj.umin(j) <= u(j, i) <= obj.umax(j));  % input constraints
                end
            end
            opti.minimize(cost);
        end

        function [] = set_problem_parameters(obj, x0, ref)
        % Set the value of the parametrized optimization problem.
            obj.opti.set_value(obj.opti_x0, x0);
            if size(ref, 2) == 1 % reference target = constant reference trajectory
                obj.opti.set_value(obj.opti_ref, repmat(ref, 1, obj.N+1));
            else % reference trajectory
                obj.opti.set_value(obj.opti_ref, ref);
            end
        end

        function [] = warm_start(obj, x0)
        % Prepare the optimization problem by reusing the previous solution. The
        % values of the previous solution are shifted, then used as the initial
        % solution to the current problem.
            prev_u = obj.prev_sol.value(obj.opti_u); prev_x = obj.prev_sol.value(obj.opti_x);
            obj.opti.set_initial(obj.opti_u, [prev_u(:, 2:end) prev_u(:, end)]);
            obj.opti.set_initial(obj.opti_x, [prev_x(:, 2:end) prev_x(:, end)]);
        end

        function [uopt, time_exec, U, X]  = solve(obj, x0, ref, varargin)
        %SOLVE Solve the optimization problem at the current step, with the
        %provided reference XREF and current state X0.
        %
        % Usage : [UOPT, TIME_EXEC, U, X] = OBJ.SOLVE(X0, XREF) returns the
        % optimal control action to apply UOPT, as well as the future control
        % inputs in U and predicted state X over the control horizon. The
        % execution time TIME_EXEC is also returned.
            tstart=tic;

            obj.set_problem_parameters(x0, ref);

            % use warm start or best guess (repeated x0) as initialization of solution
            if ~isempty(obj.prev_sol) && obj.use_warm_start
                obj.warm_start();
            else
                obj.opti.set_initial(obj.opti_x, repmat(x0, 1, obj.N+1));
            end

            % if ~isempty(varargin)
            %     obj.opti.solver('ipopt', varargin{1}, varargin{2});
            % end

            % Solve and retrieve values
            obj.prev_sol = obj.opti.solve();
            U = obj.prev_sol.value(obj.opti_u); X = obj.prev_sol.value(obj.opti_x);
            uopt = U(:, 1);
            time_exec = toc(tstart);
        end
    end
end
