%For V and S only
fun1 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)];
jac1 = @(u)[15-3.4*10^(-5).*u(1)-0.022.*u(2),-0.022.*u(1);
             0.0528.*u(2).^(0.8).*u(1).^(-0.4),-2.66.*u(2).^(0.4)+0.0704.*u(1).^(0.6).*u(2).^(-0.2)];
%For V, S and R
fun2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];
jac2 = @(u)[15-3.4*10^(-5).*u(1)-0.022.*u(2),-0.022.*u(1),0;
         0.0528.*u(2).^(0.8).*u(1).^(-0.4),-2.66.*u(2).^(0.4)+0.0704.*u(1).^(0.6).*u(2).^(-0.2)-1.4*u(3),-1.4.*u(2);
         0,0.028.*sqrt(u(3)),-1.8 + 0.014.*u(2).*1./sqrt(u(3))];

for n = 1:2
    order = -6; %Order of tolerance used in newtonsys.m
    
    if n ==1
       disp('Running for model with V and S.')
       fun = fun1;
       jac = jac1;
       start = [100000;700];
    else
        disp('Running for model with V, S and R.')
        fun = fun2;
        jac = jac2;
        start = [100000;600;20];
    end
    
    roots = start;
    tol = 10^(order);
    for i = 1:3 %Analysing for three tolerances
        h = 1;
        iter = 0;
        hvec = [];
        %Stops until tolerance is reached or iteration is 1000
        while norm(h) > tol && iter < 1000
            iter = iter+1;
            h=-jac(roots)\fun(roots);
            hvec = [hvec;norm(h)];
            roots = roots + h;
        end
        
        fprintf('Tolerance order %d\n',order)
        disp(hvec)
        
        %Printing out the final values for every tolerance
        if n ==1
            fprintf('Constants are: V=%f, S=%f\n',roots(1),roots(2));
        else
            fprintf('Constants are: V=%f, S=%f, R=%f\n',roots(1),roots(2),roots(3));
        end
        order = order-1 ;
        roots = start;
    end
end