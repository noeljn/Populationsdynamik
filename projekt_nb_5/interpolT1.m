function [y,T1] = interpolT1(y,v95)
    
    %Finding the index for value closest to 95% of max
    %and interpolating second degree polynomial with three data points
    [~,index] = min(abs(y(2,:)-v95));
    tvec = y(1,index-1:index+1)';   %Time vector (independent values)
    A=[ones(3,1),tvec,tvec.^2]; %Vandermonde matrix
    yint = y(2,index-1:index+1)';   %Vector with approximated values (dependent values)
    c = A\yint;

    f = @(x)c(1)+c(2).*x+c(3).*x.^2;    %Polynomial from interpolation
    
    %Solving the polynomial with quadratic formula (pq-formlen) find T1 
    %and displaying the T1 closest to approxiamted time from rk4()
    time1 = -c(2)/(2*c(3)) + sqrt((c(2)/(2*c(3)))^2 - (c(1)-v95)/c(3));
    time2 = -c(2)/(2*c(3)) - sqrt((c(2)/(2*c(3)))^2 - (c(1)-v95)/c(3));
    if abs(y(1,index)-time1)<abs(y(1,index)-time2)
        T1 = time1;
    else
        T1 = time2;
    end
    y = y(:,1:index);   %Extracting the datapoint between 0-T1
end