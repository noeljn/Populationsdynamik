
fun = @(v) 15*v - 1.7*10^(-5)*v^2;
vmax = 15/(1.7*10^(-5));
h=0.002; %Step length
time=[];
yvec = [];
for i = 1:5
    t = 0;
    v = 100;
    vvec = v;
    tvec = t;
    while not(v>vmax*0.99)
        k1 = fun(v);
        k2 = fun(v + h/2*k1);
        k3 = fun(v + h/2*k2);
        k4 = fun(v + h*k3);
        v = v + h/6 * (k1 + 2*k2 + 2*k3 + k4);
        t=t+h;
        tvec =[tvec;t];
        vvec = [vvec;v];
       
        
    end
    [~,k] = min(abs(tvec-0.8));
    yvec = [yvec;vvec(k)];
    h = h/2;
    
    ymatrix = [tvec,vvec]';
    [~,T1] = interpolT1(ymatrix,vmax*0.95);
    time = [time;T1];
end

%Error in truncation will occur due to large numbers
qvec = [];
for i = 1:3
    quotient = abs(yvec(i+1)-yvec(i))/abs(yvec(i+2)-yvec(i+1));
    qvec = [qvec;quotient];
end
disp('Order of accuracy')
disp(qvec)
disp('The error in time is:')
disp(abs(time(end)-time(1)))

