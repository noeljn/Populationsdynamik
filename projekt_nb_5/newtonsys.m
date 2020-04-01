<<<<<<< HEAD
function [roots]=newtonsys(tol,start,fun,jac)
    roots = start;
    h = 1;
    iter = 0;
    %Stops until tolerance is reached or iteration is 1000
    while norm(h) > tol && iter < 1000
        iter = iter+1;
        h=-jac(roots(1),roots(2))\fun(roots(1),roots(2));
        roots = roots + h;
    end
end

=======
function [roots]=newtonsys(tol,start,fun,jac)

fun = @(V,S)[15.*V-1.7*10^(-5).*V.^2-0.022.*V.*S;
             -1.9.*S.^(1.4)+0.088.*V.^(0.6).*S.^(0.8)];
jac = @(V,S)[15-3.4*10^(-5).*V-0.022.*S,-0.022.*V;
             0.0528.*S.^(0.8).*V.^(-0.4),-2.66.*S.^(0.4)+0.0704.*V.^(0.6).*S.^(-0.2)];

    roots = [880000;700];
    tol = 10^(-7);
    h = 1;
    iter = 0;
    while norm(h) > tol && iter < 5000
        h=-jac(roots(1),roots(2))\fun(roots(1),roots(2));
        roots = roots + h;
        iter = iter+1;
    end
    disp('roots')
    disp(roots)
end
>>>>>>> 4c477505fc197b078db1edbea31ae8d645388d51
