function [y,T1,v95] = interpolT1()
    vMax = 15/(1.7*10^(-5)); %Finding the constant value

    diffeq = @(v) 15*v - 1.7*10^(-5)*v^2;
    v95 = 0.95*vMax;         %Calculating 95% of the constant value

    y = rk4(vMax,100,0,diffeq,1);   %Retriving plant population with function 
                                    %for Runge-Kutta 4

    %Finding the index for value closest to v95 and interpolating
    [~,index] = min(abs(y(2,:)-v95));
    tvec = y(1,index-1:index+1)';
    A=[ones(3,1),tvec,tvec.^2];
    yint = y(2,index-1:index+1)';
    c = A\yint;

    f = @(x)c(1)+c(2).*x+c(3).*x.^2;    %Polynomial from interpolation
    
    %Solving the polynomial to find T1 and displaying the closest
    %approximated time for v95
    time1 = -c(2)/(2*c(3)) + sqrt((c(2)/(2*c(3)))^2 - (c(1)-v95)/c(3));
    time2 = -c(2)/(2*c(3)) - sqrt((c(2)/(2*c(3)))^2 - (c(1)-v95)/c(3));
    if abs(y(1,index)-time1)<abs(y(1,index)-time2)
        T1 = time1;
    else
        T1 = time2;
    end
    y = y(:,1:index);   %Extracting the datapoint between 0-T1
end