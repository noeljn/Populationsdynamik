function [y] = rk4(v0)

format long
vMax = 15/(1.7*10^(-5));

h = 0.02;
t = 0;
v = v0;
y = [];

disp(y);
f = @(t,v) 15*v - 1.7*10^(-5)*v^2;


while vMax > v
   y(end + 1, 1) = v;
   y(end, 2) = t;
   k1 = f(t,v);
   k2 = f(t + h/2,v + h/2*k1);
   k3 = f(t + h/2,v + h/2*k2);
   k4 = f(t, v + h*k3);
   add = h/6 * (k1 + 2*k2 + 2*k3 + k4);
   v = v + add;
   t = t + h;
   
end
%plot(y(:,2),y(:,1));
end
