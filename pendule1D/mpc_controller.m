function u = mpc_controller(q, Ad, Bd, Q_bar, R_bar, F, G, H_qp, lb, ub, N)
    % Vecteur gradient (dépend de l'état courant)
    f_qp = (G' * Q_bar * F * q);

    % Résolution QP
    opts = optimoptions('quadprog', 'Display', 'off', ...
                        'MaxIterations', 200);
    U_opt = quadprog(H_qp, f_qp, [], [], [], [], lb, ub, [], opts);

    % On n'applique que le premier élément
    if isempty(U_opt)
        u = 0;  % repli si QP infaisable
    else
        u = U_opt(1);
    end
end