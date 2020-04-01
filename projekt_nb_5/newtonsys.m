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
