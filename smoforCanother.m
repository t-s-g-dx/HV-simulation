close all
clear
clc
% define initial conditions and time range
tspan = 0:0.01:10; % time range from t=0 to t=10
x0 = [1; 1; 3; 1; 1; 1];% initial condition 
z0=[0;0;0;0;0;0;0;0;0];
lam=[3;1;1;3;1;1;3;1;1];
lam=1.*lam;
M=5;
I=0.5;
options = odeset('RelTol', 1e-6);
% options = odeset('MaxStep', 1e-4);
% Solve the differential equation
[tx, x] = ode45(@(tx,x) myODE(tx,x,M,I), tspan, x0,options);
tic
[t,z] = ode45(@(t,z) odefcn(t,z,tx,x,lam), tspan,z0,options);
elapsed_time = toc;
disp(['program runing time：', num2str(elapsed_time), 's']);

F=2*sin(1*tx);
T=0.5*sin(1*tx);
Fx=z(:,3).*M./cos(atan2(x(:,2),x(:,1)));
Fy=z(:,6).*M./sin(atan2(x(:,2),x(:,1)));
Tt=z(:,9).*I;
Fe=z(:,3).*M./cos(atan2(x(:,2),x(:,1)))-F;
Fe2=z(:,6).*M./sin(atan2(x(:,2),x(:,1)))-F;
Te=z(:,9).*I-T;
meanF=mean(Fe.^2)

e=zeros(size(x));
e(:,1)=z(:,1)-x(:, 1);
e(:,2)=z(:,4)-x(:, 2);
e(:,3)=z(:,7)-x(:, 3);
e(:,4)=z(:,2)-x(:, 4);
e(:,5)=z(:,5)-x(:, 5);
e(:,6)=z(:,8)-x(:, 6);



% draw the result
plot(tx, x(:, 1), tx, x(:, 2),tx,x(:,3),tx,x(:,4),tx,x(:,5),tx,x(:,6));
legend('x1(t)', 'x2(t)','x3(t)','x4(t)','x5(t)','x6(t)');
xlabel('t');
ylabel('x');
figure
plot(t, z(:, 1), t, z(:, 4),t,z(:,7),t,z(:,2),t,z(:,5),t,z(:,8));
legend('z1(t)', 'z4(t)','z7(t)','z2(t)','z5(t)','z8(t)');
xlabel('t');
ylabel('z');
figure
plot(t,Fe,t,Te,t,Fe2);
figure
plot(t, e(:, 1), t, e(:, 2),t,e(:,3),t,e(:,4),t,e(:,5),t,e(:,6));
legend('e1(t)', 'e4(t)','e7(t)','e2(t)','e5(t)','e8(t)');
xlabel('t');
ylabel('z');


% Define the differential equation function
function dxdt = myODE(tx,x,M,I)
    F=2*sin(1*tx);
    T=0.5*sin(1*tx);
%     F=1;
%     T=0;
    dxdt = zeros(size(x));
    dxdt(1) = x(4);
    dxdt(2) = x(5);
    dxdt(3)=x(6);
    dxdt(4)=cos(atan2(x(2),x(1)))/M*F;
    dxdt(5)=sin(atan2(x(2),x(1)))/M*F;
    dxdt(6)=T/I;
end
function dzdt = odefcn(t,z,tx,x,lam)
  x = interp1(tx,x,t);
  dzdt = zeros(size(z));
  v1=-lam(1)*abs(z(1)-x(1)).^(2/3)*sign(z(1)-x(1))+z(2);
  dzdt(1) = v1;
  v2=-lam(2)*abs(z(1)-x(1)).^(1/3)*sign(z(1)-x(1))+z(3);
  dzdt(2) = v2;
  dzdt(3)=-lam(3)*sign(z(1)-x(1));

  v3=-lam(4)*abs(z(4)-x(2)).^(2/3)*sign(z(4)-x(2))+z(5);
  dzdt(4) = v3;
  v4=-lam(5)*abs(z(4)-x(2)).^(1/3)*sign(z(4)-x(2))+z(6);
  dzdt(5) = v4;
  dzdt(6)=-lam(6)*sign(z(4)-x(2));

  v5=-lam(7)*abs(z(7)-x(3)).^(2/3)*sign(z(7)-x(3))+z(8);
  dzdt(7) = v5;
  v6=-lam(8)*abs(z(7)-x(3)).^(1/3)*sign(z(7)-x(3))+z(9);
  dzdt(8) = v6;
  dzdt(9)=-lam(9)*sign(z(7)-x(3));
 
end