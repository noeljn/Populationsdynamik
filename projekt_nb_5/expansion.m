function harvest = expansion(dataT3,tau)

func2 = @(u)[15.*u(1)-1.7*10^(-5).*u(1).^2-0.022.*u(1).*u(2);
             -1.9.*u(2).^(1.4)+0.088.*u(1).^(0.6).*u(2).^(0.8)-1.4.*u(2).*u(3);
             -1.8.*u(3) + 0.028.*u(2).*sqrt(u(3))];

%Collecting data between T3-T4 with Runge-Kutta 4 
t = 3;
start = [dataT3(1:2,end);dataT3(3,end)*0.3;dataT3(4,end)*0.8];
data = start;
while t < 8
    
    %Collecting data between t and tau
    year = rk4(t + tau,start(2:4,end),start(1,end),func2,2);
    %Harvest
    year(2,end) = 100;
    data = [data(:,1:end-1),year];
    start = [data(1:2,end);data(3,end);data(4,end)];

    %Pest
    t = t + 1;
    year = rk4(t,start(2:4,end),start(1,end),func2,2);
    data = [data(:,1:end-1),year];
    start = [data(1:2,end);data(3,end)*0.3;data(4,end)*0.8];
    
    
end


%{
figure(1)
plot(data(1,:),data(2,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
figure(2)
plot(data(1,:),data(3,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
figure(3)
plot(data(1,:),data(4,:),'color',[0.9568, 0.5019, 0.1254]);
hold on
%}
%Returns harvest for year 7
index = find(data(1,:) == (7 + tau));
harvest = data(2,index - 1) - 100;

end