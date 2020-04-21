clear all
format long
fun = @(v) 15*v - 1.7*10^(-5)*v^2;

func1 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)];
func2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
time=[];

start = 100;
limit = 1.5;
t0 = 0;
s = 1;
h=0.002; %Step length

for n = 1:3
    hvec = [];
    yvec = [];
    N = 1000;   %No. of part intervals used in rk4()
    for i = 1:5
        t = t0;
        v = start;
        vvec = v;
        tvec = t;
        
        if s == 1   %For when finding T1
            vmax = 15/(1.7*10^(-5));
            
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
            hvec = [hvec;h];
            [~,k] = min(abs(tvec-0.8));
            yvec = [yvec;vvec(k)];  %Storing approximated V, 95% of max
            h = h/2;
            
            ymatrix = [tvec,vvec]';
            [~,T1] = interpolT1(ymatrix,vmax*0.95);
            time = [time;T1];
            
        else    %For when searching values at T2 and T3
            h = (limit-t0)/N;   %Calculating the step length
            tvec = t0:h:limit;
            y=[t0;v];
            
            for j = 1:N
                k1 = fun(v);
                k2 = fun(v + h/2*k1);
                k3 = fun(v + h/2*k2);
                k4 = fun(v + h*k3);
                v = v + h/6 * (k1 + 2*k2 + 2*k3 + k4);
                vvec = [vvec,v];
            end
            yvec = [yvec;vvec(1,end)];
            N = N*2;
            hvec = [hvec;h];
        end
        
        %Printing out the values every iteration and dividing the step size
        if n == 2
            fprintf('Approximated values are: V=%f, S=%f\n',vvec(1,end),vvec(2,end));
        elseif n == 3
            fprintf('Approximated values are: V=%f, S=%f, R=%f\n',vvec(1,end),vvec(2,end),vvec(3,end));
        end
    end
    %Prints the error in time T1 between the first and last time found
    if n == 1
        disp('The error in time is:')
        disp(abs(time(end)-time(1)))
    end
    
    %This prints out the order of accuracy that can be calculated
    %Error in truncation will occur due to large numbers
    qvec = [];
    for i = 1:3
        quotient = abs(yvec(i+1)-yvec(i))/abs(yvec(i+2)-yvec(i+1));
        qvec = [qvec;quotient];
    end
    disp('Order of accuracy')
    disp(qvec)
    
    
    %Plotting the error against step size and reference 
    %to see if RK has order of 4
    felv = abs(yvec(1:end-1)-yvec(2:end));
    figure(10+n); %Starts from 10 to not overwrite other graphs
    loglog(hvec(1:end-1),felv,'b'); 
    hold on
    loglog(hvec,hvec.^4,'k')
    hold on
    legend('Error','Reference','Location','southeast')
    hold on
    xlabel('Step size')
    hold on
    ylabel('Error')
    
    %Changing function to be analysed and initial values
    if n == 1
        fun = func1;
        start = [0.95*15/(1.7*10^(-5));2];
        t0 = time(1);
        s = 2;
    elseif n == 2
        fun = func2;
        start = [vvec(1:2,end);2];
        t0 = 1.5;
        limit = 3;
        s = 2;
    end
    
end
figure(11)
title('Analysing Runge-Kutta 4 when finding T1')
figure(12)
title('Analysing Runge-Kutta 4 for V & S between T = T1 and 1.5')
figure(13)
title('Analysing Runge-Kutta 4 for V, S & R between T = 1.5 and 3')

