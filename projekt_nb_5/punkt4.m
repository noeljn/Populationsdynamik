function [y] = punkt4(start)    
    format long
    h = 0.02;
    tMax = 1.5;
    t = start(1);
    v = [start(2);2];
    y = []; %Containing [time, plant population]
    func = @(t,v)[15.*v(1)-1.7*10^(-5).*v(1).^2-0.022.*v(1).*v(2);
             -1.9.*v(2).^(1.4)+0.088.*v(1).^(0.6).*v(2).^(0.8)];
    y(1,:) = v;


    %Uv(2)ing Runge-Kutta 4 until v(1) reachev(2) the conv(2)tant v(1)alue
    while t < tMax
       %Runge-Kutta 4
       k1 = func(t,v);
       k2 = func(t + h/2,v + h/2*k1);
       k3 = func(t + h/2,v + h/2*k2);
       k4 = func(t, v + h*k3);
       v = v + h/6 * (k1 + 2*k2 + 2*k3 + k4);
       t = t + h;
       y(end + 1,1:2) = v';
       y(end,3) = t;

    end
    %plot(y(:,3),y(:,1));
    %hold on;
    %plot(y(:,3),y(:,2));
    y = [y(end,1),y(end,2)];
end