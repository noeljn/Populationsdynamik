function [roots]=newtonsys(start,fun,jac)
    roots = start;
    h = 1;
    iter = 0;
    %Stops until tolerance is reached or iteration is 1000
    while norm(h) > 10^(-6) && iter < 1000
        iter = iter+1;
        h=-jac(roots)\fun(roots);
        roots = roots + h;
    end
end
