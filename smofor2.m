close all
clear 
tic
n=1;%number of cars
Nx=(n+1)*6;%number of states
Nz=(n)*6+(n)*3+6;%number of estimate states
tspan = 0:0.01:10; % time range from t=0 to t=10

% initial condition 
index=1;
car=[1 1 1 0 0 0; 2 2 2 0 0 0];
x0 = zeros(Nx,1);
for i=1:6:Nx
    c=car(:,index);
    x0(i:i+5,1)=car(index,:);
    index=index+1;
end
z0=zeros(Nz,1);
lam=ones(Nz,1);
lam(1)=3;
lam(2)=1;
lam(3)=1;
lam(4)=3;
lam(5)=1;
lam(6)=1;
lam(7)=3;
lam(8)=1;
lam(9)=1;
lam(Nz-5)=3;
lam(Nz-4)=5;
lam(Nz-3)=3;
lam(Nz-2)=5;
lam(Nz-1)=3;
lam(Nz)=5;
M=5;
I=0.5;
options = odeset('RelTol', 1e-4);
% options = odeset('MaxStep', 1e-4);
% Solve the differential equation
[tx, x] = ode45(@(tx,x) myODE(tx,x,M,I), tspan, x0,options);

[t,z] = ode45(@(t,z) odefcn(t,z,tx,x,lam), tspan,z0,options);


F=2*sin(1*tx);
T=0.5*sin(1*tx);
Fx=(z(:,3)+z(:,Nz-4)).*M./cos(atan2(x(:,2),x(:,1)));
Fy=(z(:,6)+z(:,Nz-2)).*M./sin(atan2(x(:,2),x(:,1)));
Tt=(z(:,9)+z(:,Nz)).*I;
e1=x(:,1)-x(:,Nx-5);
e2=z(:,1);
e=e2-e1;
e_1=(x(:,3)-x(:,Nx-3))-z(:,7);
e_2=(x(:,2)-x(:,Nx-4))-z(:,4);
Fe=(z(:,3)+z(:,Nz-4)).*M./cos(atan2(x(:,2),x(:,1)))-F;
Te=(z(:,9)+z(:,Nz)).*I-T;
meanF=mean(Fe.^2)
% draw the result
figure
hold on
for i=1:Nx
    plot(tx,x(:,i))
    legendInfo{i}=['x',num2str(i)];
end
legend(legendInfo);
xlabel('t');
ylabel('x');
figure
plot(t,e,t,e_1,t,e_2)
figure
plot(t,Fe,t,Te)
elapsed_time = toc;
disp(['program runing time：', num2str(elapsed_time), 's']);
function dxdt = myODE(tx,x,M,I)
     F=2*sin(1*tx);
    T=0.5*sin(1*tx);
%      F=1;
%      T=0.1;
    nx=size(x);
    dxdt = zeros(nx);
    for i=1:6:nx
    dxdt(i)=x(i+3);
    dxdt(i+1)=x(i+4);
    dxdt(i+2)=x(i+5);
    dxdt(i+3)=cos(atan2(x(i+1),x(i)))/M*F;
    dxdt(i+4)=sin(atan2(x(i+1),x(i)))/M*F;
    dxdt(i+5)=T/I;
    end
end

function dzdt = odefcn(t,z,tx,x,lam)
  x = interp1(tx,x,t);
  dzdt = zeros(size(z));
  nz=size(z,1);
  nx=size(x,2);

  %C Acceleration Estimation
  v1=-lam(nz-5)*abs(z(nz-5)-x(nx-2)).^(1/2)*sign(z(nz-5)-x(nx-2))+z(nz-4);
  dzdt(nz-5) = v1;
  dzdt(nz-4)=-lam(nz-4)*sign(z(nz-4)-v1);
  v2=-lam(nz-3)*abs(z(nz-3)-x(nx-1)).^(1/2)*sign(z(nz-3)-x(nx-1))+z(nz-2);
  dzdt(nz-3) = v2;
  dzdt(nz-2)=-lam(nz-2)*sign(z(nz-2)-v2);
  v3=-lam(nz-1)*abs(z(nz-1)-x(nx)).^(1/2)*sign(z(nz-1)-x(nx))+z(nz);
  dzdt(nz-1)=v3;
  dzdt(nz)=-lam(nz)*sign(z(nz)-v3);
  i=1;
  j=1;
  while i<(nz-6)&&j<(nx-5)
  v1=-lam(i)*abs(z(i)-(x(j)-x(nx-5))).^(2/3)*sign(z(i)-(x(j)-x(nx-5)))+z(i+1);
  dzdt(i) = v1;
  v2=-lam(i+1)*abs(z(i+1)-v1).^(1/2)*sign(z(i+1)-v1)+z(i+2);
  dzdt(i+1) = v2;
  dzdt(i+2)=-lam(i+2)*sign(z(i+2)-v2);

  v3=-lam(i+3)*abs(z(i+3)-(x(j+1)-x(nx-4))).^(2/3)*sign(z(i+3)-(x(j+1)-x(nx-4)))+z(i+4);
  dzdt(i+3) = v3;
  v4=-lam(i+4)*abs(z(i+4)-v3).^(1/2)*sign(z(i+4)-v3)+z(i+5);
  dzdt(i+4) = v4;
  dzdt(i+5)=-lam(i+5)*sign(z(i+5)-v4);

  v5=-lam(i+6)*abs(z(i+6)-(x(j+2)-x(nx-3))).^(2/3)*sign(z(i+6)-(x(j+2)-x(nx-3)))+z(i+7);
  dzdt(i+6) = v5;
  v6=-lam(i+7)*abs(z(i+7)-v5).^(1/2)*sign(z(i+7)-v5)+z(i+8);
  dzdt(i+7) = v6;
  dzdt(i+8)=-lam(i+8)*sign(z(i+8)-v6);
  i=i+9;
  j=j+6;
  end
end