function u = mpc_blk(q)
    % Récupère les matrices pré-calculées depuis le workspace
    persistent Ad Bd Q_bar R_bar F_mat G_mat H_qp lb ub N_hor
    if isempty(Ad)
        Ad    = evalin('base', 'Ad');
        Bd    = evalin('base', 'Bd');
        Q_bar = evalin('base', 'Q_bar');
        R_bar = evalin('base', 'R_bar');
        F_mat = evalin('base', 'F');
        G_mat = evalin('base', 'G');
        H_qp  = evalin('base', 'H_qp');
        lb    = evalin('base', 'lb');
        ub    = evalin('base', 'ub');
        N_hor = evalin('base', 'N');
    end

    f_qp  = G_mat' * Q_bar * F_mat * q;
    U_opt = quadprog(H_qp, f_qp, [], [], [], [], lb, ub, [], ...
                     optimoptions('quadprog','Display','off'));
    if isempty(U_opt)
        u = 0;
    else
        u = U_opt(1);
    end
end