function reliability()
    fun = @(v) 15*v - 1.7*10^(-5)*v^2;
    t = 0;
    v = 100;
    limit = 0.99*(15/(1.7*10^(-5)));
    h=0.02; %Step length
    yvec = [];
    vvec = [];
    for i = 1:3
        while  not(v > limit)
            k1 = fun(v);
            k2 = fun(v + h/2*k1);
            k3 = fun(v + h/2*k2);
            k4 = fun(v + h*k3);
            add = h/6 * (k1 + 2*k2 + 2*k3 + k4);
            v = v + h/6 * (k1 + 2*k2 + 2*k3 + k4);
            vvec = [vvec;v];
            t=t+h;
        end
    yvec = [yvec;v];
    h = h/2;
    vvec = [];
    v=100;
    end
end