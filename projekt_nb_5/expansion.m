function expansion(dataT3,tau)

func2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];

%Collecting data between T3-T4 with Runge-Kutta 4         
data = rk4(4,[dataT3(2:3,end);2],dataT3(1,end),func2,2);



figure(1)
plot(data(1,:),data(2,:),'color',[0.6666,0.2,0.3490]);
hold on
figure(2)
plot(data(1,:),data(3,:),'color',[0.6666,0.2,0.3490]);
hold on
figure(3)
plot(data(1,:),data(4,:),'color',[0.6666,0.2,0.3490]);
hold on
end