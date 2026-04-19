function [U, X, time_exec, progress_window] = run_QTP_with_controller(reg, x0, xref, t, dt, uref)
%RUN_QTP_WITH_CONTROLLER Simulate the closed-loop behaviour of the four-tanks
%benchmark system with controller REG.
%
% Usage: [U, X, TIME_EXEC, PROGRESS_WINDOW] = run_QTP_with_controller(REG, X0,
% XREF, TIME_ARRAY, PREV_U) returns the series of inputs U and states X of the
% closed-loop simulation. TIME_EXEC is the average time to compute control
% inputs and PROGRESS_WINDOW is the handle to a window showing the progress of
% computation.
%
%  REG should be a function (or callable object) that accepts 2 arguments X
%  and XREF which are the current state and reference target (vector) or
%  trajectory (matrix).
%
%  X0 is a vector corresponding to the initial state of the system when
%  beginning the simulation.
%
%  XREF is the reference target (as a vector) or reference trajectory (as a
%  matrix).
%
%  T is a vector between indicating the time of start and end of simulation.
%
%  DT is the sampling time of the regulator in seconds.
%
%  (optionnal) UREF is the reference input.
import casadi.*
f = QTP_toricelli();
phi = discretize_model(f, 'dt', t(2)-t(1));

l = length(t);
n_iter = 0;
X = [x0 zeros(4, l-1)];
time_exec = 0;
msg = 'Simulation du système en boucle fermée !';
progress_window=waitbar(0,msg);
u = reg(x0, xref);
U = [u zeros(2, l-1)];
for k=1:l-1
    try
        waitbar(k/l, progress_window, msg);
    catch
        msg = 'Simulation du système en boucle fermée... Merci de ne pas me fermer !';
        progress_window=waitbar(k/l, msg);
    end

    if mod(t(k+1), dt) == 0
        try
            [u, time_compute] = reg(X(:, k), xref);
        catch
            warning('MPC failed')
            break
        end
        time_exec=time_exec+time_compute;
        n_iter=n_iter+1;
    end
    U(:, k+1) = u;
    X(:, k+1) = full(phi(X(:, k), double(U(:, k))));  % apply input to the QTP and store it
end

waitbar(1, progress_window, ['Simulation terminée, temps moyen par itération: ' num2str(time_exec/n_iter, '%f') 's']);
