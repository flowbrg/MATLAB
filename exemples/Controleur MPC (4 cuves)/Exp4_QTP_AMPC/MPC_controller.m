classdef MPC_controller < matlab.System
    %MPC_controller Class for setting up an MPC controller
    %
    % Usage:
    %
    % mpc = MPC_controller(PHI, N, Q, R, P, XMIN, XMAX, UMIN, UMAX) creates and
    % returns an MPC object with state evolution function :
    %                    x_k+1 = PHI(x_k, u_k)
    %
    % over a control horizon of length N and quadratic cost function:
    %            V_k = (x_k-xref_k)'*Q*(x_k-xref_k) + u_k'*R*u_k.
    %
    % P is the weight on the terminal cost (x_N-xref_N)'P(x_N-xref_N).
    % XMIN, XMAX, UMIN, UMAX are respectively the constraints on the lower
    % state bound, upper state bound, lower input bound and upper input bound.
    properties (Nontunable)
        model;           % system evolution function
        N;               % Prediction horizon step
        Q;               % state stage cost weight matrix x_i^T Q x
        R;               % input stage cost weight matrix u_i^T R u
        P;               % terminal cost weight matrix x_i^T P x
        xmin;            % vector of min value for each state
        xmax;            % vector of max value for each state
        umin;            % vector of min value for each input
        umax;            % vector of max value for each input
        dt;              % sampling rate
        nx;              % size of state
        nu;              % size of input
    end

    properties(Access=private)
        opti = [];       % Casadi opti object
        opti_x = [];     % (casadi.opti.variable) symbol of the states decision variable
        opti_x0 = [];    % (casadi.opti.parameter) symbol of the initial state
        opti_ref = [];   % (casadi.opti.parameter) symbol of reference target/trajectory
        opti_u = [];     % (casadi.opti.variable) symbol of input decision variable
        opti_u0 = [];    % (casadi.opti.parameter) symbol of last input decision 
        u = [];     % previous input decision
        prev_sol = [];
        use_warm_start=true; % whether to use warm start. This is useful to toggle it internally (RTI)
    end

    methods

        function obj = MPC_controller(varargin)
            setProperties(obj, nargin, varargin{:});
        end %function (constructor MHE)



        function [opti, x, x0, u0, xref, u] = construct_problem(obj)

        %==================================================================
        % CONSTRUCT_PROBLEM initialize the MPC optimization problem
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

            opti = casadi.Opti('conic');
            opts = struct();
            opts.print_time = false;
            opts.printLevel = 'none';
            opti.solver('qpoases',opts);

            % Define decision variables
            % **Define state variable x here**
            x = opti.variable(obj.nx, obj.N + 1);   % [4 x (N+1)]
            % **Define input variable u here**
            u = opti.variable(obj.nu, obj.N);        % [2 x N]

            % Define parameter variables (Optimization parameters, not model parameters!)
            % **Define initial state here**
            x0   = opti.parameter(obj.nx, 1);        % [4 x 1]
            u0   = opti.parameter(obj.nu, 1);        % [2 x 1]

            % **Constraint your state trajectory with initial state**
            opti.subject_to( x(:, 1) == x0 );

            % **Define your state reference variable here**
            xref = opti.parameter(obj.nx, 1);        % [4 x 1]
            
            cost = 0;
            for i=1:obj.N
                % **Define your objective function here**
                % Error w.r.t. reference
                ex = x(:, i) - xref;           % state error
                eu = u(:, i);                  % input (penalise absolute value
                                               % or deviation from uref if needed)
                % Stage cost  (LQ-style)
                % J = int(x'Qx + u'Ru)
                cost = cost + ex'*obj.Q*ex + eu'*obj.R*eu;

                % **Define your model constraint here**
                % --- Equality constraint: nonlinear model dynamics ------
                % obj.model is a casadi.Function: x_next = model(x_k, u_k)
                x_next = obj.model(x(:, i), u(:, i), x0, u0);
                opti.subject_to(x(:, i+1) == x_next);

                % **Define your state limit constraint here**
                % State bounds
                opti.subject_to(obj.xmin <= x(:, i));
                opti.subject_to(x(:, i) <= obj.xmax);

                % **Define your input limit constraint here**                
                % Input bounds
                opti.subject_to(obj.umin <= u(:, i));
                opti.subject_to(u(:, i) <= obj.umax);
            end

            % **Define your terminal cost here**
            ex_N = x(:, obj.N+1) - xref;
            cost = cost + ex_N' * obj.P * ex_N;
            % **Set objective function to your optimization**
            opti.minimize(cost);
        %==================================================================
        end

        function [uopt, time_exec, U, X]  = solve(obj, x0, ref, varargin)
        %==================================================================
        % SOLVE solve the optimization problem at the current step, with the
        % provided reference XREF and current state X0.
        %
        % Usage : [UOPT, TIME_EXEC, U, X] = OBJ.SOLVE(X0, XREF) returns the
        % optimal control action to apply UOPT, as well as the future control
        % inputs in U and predicted state X over the control horizon. The
        % execution time TIME_EXEC is also returned.
            tstart=tic;
    
            obj.opti.set_initial(obj.opti_x,repmat(x0,1,obj.N+1));
            obj.opti.set_initial(obj.opti_u,zeros(obj.nu,obj.N));

            % Set optimization parameters
            % **Set feedback as initial state**
            obj.opti.set_value(obj.opti_x0,  x0(:));    % current state
            obj.opti.set_value(obj.opti_u0, obj.u(:));  % current input
            obj.opti.set_value(obj.opti_ref, ref(:));   % reference

            % **Set state reference**
            if obj.use_warm_start && ~isempty(obj.prev_sol)
                try
                    obj.opti.set_initial(obj.prev_sol.value_variables());
                catch
                    disp("WARM START FAILED - using defaut")
                end
            end

            % Solve and retrieve values
            % **Solve optimization**
            sol = obj.opti.solve();
            obj.prev_sol = sol;

            % Extract solution
            % **Extract result for each variable**
            X = sol.value(obj.opti_x);   % [nx x (N+1)]  predicted states
            U = sol.value(obj.opti_u);   % [nu x  N   ]  optimal inputs

            % **Take first time step u vector as MPC output**
            uopt = U(:, 1);              % [nu x 1]
            obj.u = uopt;
            time_exec = toc(tstart);
        %==================================================================
        end

    end

    methods (Access = protected)
        % MPC implementation
        function setupImpl(obj)
            % Matlab class execute command here before starting simulation.
            [obj.opti, obj.opti_x, obj.opti_x0, obj.opti_u0, obj.opti_ref, obj.opti_u] = obj.construct_problem();
        end

        function u = stepImpl(obj, xref, x)
            % Matlab class execute command here during simulation.
            if isempty(obj.u)
                obj.u = zeros(obj.nu, 1);  % ou obj.uref si vous l'ajoutez en propriété
            end
            u = obj.solve(x, xref);
        end
        
        % Other setup
        function fixed_out = isOutputFixedSizeImpl(obj)
            fixed_out = true;
        end%function

        function size_out = getOutputSizeImpl(obj)
            size_out = obj.nu;
        end %function

        function type_out = getOutputDataTypeImpl(obj)
            type_out = "double";
        end %function

        function complex_out = isOutputComplexImpl(obj)
            complex_out = false;
        end%function

        function resetImpl(obj)
        %    obj.counter = 0;
        end %function

        function releaseImpl(obj)
        %    obj.counter = 0;
        end

        function sts = getSampleTimeImpl(obj)
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.dt);
        end %function

        function icon = getIconImpl(obj)
            % Set appearance of simulink block
            icon = ['AMPC controller'];
        end %function
    end %protected methods

end
