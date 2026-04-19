function x = mySTksi(x,u)
    w = 1;
    K = 1;
    dt = 0.01;

    x = x + [x(2); -w^2*x(1)-2*w*x(2)*x(3)+K*w^2*u;0]*dt;
end