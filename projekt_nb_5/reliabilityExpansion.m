function  reliabilityExpansion(dataT3)
yvec = [];
for i = 0:4
    N = 1000 *2^i;
    %By changing the tau and pest u can check how many decimals are correct
    %for each tau.
    data = expansion(dataT3,0.5,1,N);
    yvec(end + 1) = data(2,end);
    
    %You can disp the amount to check how many decimals are wrong
    %disp(data(2,end - 1));
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
    
end


%Modified expansion so u can change the N for rk4
function [data] = expansion(dataT3,tau,pest,N)
func2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];

%Collecting data between T3-T8 with Runge-Kutta 4 
t = 4;
if pest == 1
    start = [dataT3(1:2,end);dataT3(3,end)*0.3;dataT3(4,end)*0.8];
else
    start = [dataT3(1:2,end);dataT3(3,end);dataT3(4,end)];
end
data = start;

while t < 8
    
    %Collecting data between t and tau
        year = rk4Test(t + tau,start(2:4,end),start(1,end),func2,2,N);
        %Harvest
        year(2,end) = 100;
        data = [data(:,1:end-1),year];
        start = [data(1:2,end);data(3,end);data(4,end)];
        
        if tau == 0
            start = [data(1:2,end);data(3,end)*0.3;data(4,end)*0.8];
        end

    %Pest
    t = t + 1;
    if pest == 1 && tau ~= 0
        year = rk4Test(t,start(2:4,end),start(1,end),func2,2,N);
        data = [data(:,1:end-1),year];
        start = [data(1:2,end);data(3,end)*0.3;data(4,end)*0.8];
    end
    
end

end


%Modified so u can change the N 
function [y] = rk4Test(limit,v0,t0,fun,s,N)
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

