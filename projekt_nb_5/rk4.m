function [y] = rk4(limit,v0,t0,fun,s)
    format long   
    t = t0;
    v = v0;
    y = []; %Containing [time, data points]
   
    if s==1 %When calculating T1 and its unknown
        h=0.002; %Step length
        while not(v > limit)
            y = [y,[t;v]];
            v = calcRK(fun,v,h);
            t=t+h;
        end
    else    %When the time span is known
        N = 1000;   %No. of part intervals
        h = (limit-t0)/N;   %Calculating the step length
        tvec = t0:h:limit;
        y=[t0;v];
        for i = 1:N
           v = calcRK(fun,v,h);
           y = [y,[tvec(i+1);v]];
        end
    end
    
end

%Runge-Kutta 4
function [v] = calcRK(fun,v,h)
    
    k1 = fun(v);
    k2 = fun(v + h/2*k1);
    k3 = fun(v + h/2*k2);
    k4 = fun(v + h*k3);
    v = v + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

