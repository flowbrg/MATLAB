function [xeq, ueq] = find_equilibrium(xref)
% find equilibrium of system at xref (reference for reservoirs 1 and 2)
    xref = reshape(xref, [], 1);
    import casadi.*
    f = QTP_toricelli();
    xueq = SX.sym('xueq', 4);
    rfp = struct('x', xueq, 'g', f([xref;xueq(1:2)], xueq(3:4))');
    s = casadi.rootfinder('solver', 'newton', rfp);%, struct('nlpsol', 'ipopt'));
    sol = full(s([xref; zeros(2, 1)], 0));
    xeq = [xref; sol(1:2)];
    ueq = sol(3:4);
end