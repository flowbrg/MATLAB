function phi = discretize_model(f, options)
    % Discretize f using options.METHOD and time interval options.DT.
    % f should take 2 arguments, such that xdot = f(x, u).
    arguments
        f
        options.method = 'RK4'; % default
        options.dt (1, 1) double = 0.1;
    end
    import casadi.*

    names = cellstr(name_in(f));
    x = names{1}; u = names{2};
    state = SX.sym(x, size1_in(f, x), size2_in(f, x));
    controls = SX.sym(u, size1_in(f, u), size2_in(f, u));

    switch options.method
      case 'RK4'
        k1 = f(state, controls);
        k2 = f(state+options.dt*k1/2.0,controls);
        k3 = f(state+options.dt*k2/2.0,controls);
        k4 = f(state+options.dt*k3,controls);
        states_1 = (state+options.dt*(k1+2*k2+2*k3+k4)/6.0)';
        phi = Function('phi', {state,controls},{states_1});

      case 'euler'
        states_1 = (state+options.dt*f(state, controls))';
        phi = Function('phi', {state,controls},{states_1});
    end
end