function [y] = rk4(limit,v0,fun)
    format long
    h = 0.02;   %Length of steps
    t = 0;
    v = v0;
    y = [];     %Containing [time, plant population]

    %Using Runge-Kutta 4 until v reaches the constant value
    while v < limit
       y(end + 1, 2) = v;
       y(end, 1) = t;
       
       %Runge-Kutta 4
       k1 = fun(t,v);
       k2 = fun(t + h/2,v + h/2*k1);
       k3 = fun(t + h/2,v + h/2*k2);
       k4 = fun(t, v + h*k3);
       v = v + h/6 * (k1 + 2*k2 + 2*k3 + k4);
       t = t + h;

    end
end