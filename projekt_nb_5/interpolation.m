v95 = 0.95*15/(1.7*10^(-5));

y = rk4(100);
disp(y);


[c,index] = min(abs(y(:,1)-v95));
tvec = y(index-1:index+1,2)
A=[ones(3,1),y(index-1:index+1,2),y(index-1:index+1,2).^2];
yint = y(index-1:index+1,1);

disp(A)
disp(yint)

c = A\yint;

f = @(x)c(1)+c(2).*x+c(3).*x.^2;
xvec=linspace(tvec(1),tvec(end));
plot(f(xvec))
hold on 
plot(y(index-1:index+1,2),y(index-1:index+1,1))