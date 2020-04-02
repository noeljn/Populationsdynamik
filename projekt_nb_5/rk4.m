function [y] = rk4(limit,v0,t0,fun,s,h)
    format long   
    t = t0;
    v = v0;
    y = []; %Containing [time, plant population]
    check = 0;

    %Using Runge-Kutta 4 until v reaches the constant value
    while check < limit
        
       y(end + 1, 2) = v;
       y(end, 1) = t;
        
       %Runge-Kutta 4
       k1 = fun(t,v);
       k2 = fun(t + h/2,v + h/2*k1);
       k3 = fun(t + h/2,v + h/2*k2);
       k4 = fun(t, v + h*k3);
       v = v + h/6 * (k1 + 2*k2 + 2*k3 + k4);
       t = t + h;
       
       if s == 1
           check = v;
       else
           check = t;
       end

    end
end