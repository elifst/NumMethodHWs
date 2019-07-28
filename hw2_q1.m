method = input('Please choose the method you want to use( 1 for FTCS, 2 for Crank Nicholson )  :  ');
switch method
    case 1
%FTCS expilicit method
clear ,clc
t=40; %time in seconds
dt=0.04; %time step
v=1/(pi^2);
M=11; %number of nodes
x0=0; xL=1; x=xL-x0;
dx=x/(M-1);
N=t/dt+1; %number of nodes in time axis
T=zeros(N,M);
%Tex=zeros(N,M);
%initial conditions:
for j=1:M
    x=(j-1)*dx;
    T(N,j)=cos(pi*(x-0.5));
end
%boundary conditions:
for i=1:N
    T(i,1)=0;
    T(i,M)=0;
end
%numeric solutions:
for j=N:-1:2
    for i=2:M-1
        T(j-1,i)=T(j,i)+(v*dt/(dx^2))*(T(j,i+1)-2*T(j,i)+T(j,i-1));
    end
end
%exact solutions:
Tex=zeros(N,M);
for i=0:M-1
   x=i*dx;
   for j=0:N-1
       tx=j*dt;
       Tex(N-j,i+1)=exp(-tx)*cos(pi*(x-0.5));
   end
end
figure
X=0:dx:1;
for i=N:-1:1
    t=dt*(N-i);
    if t==0 || t==0.12 || t==0.2 || t==0.32 || t==0.4 || t==4 || t==40
    hold on
    plot(X,T(i,:),'color','black','linewidth',1);
    plot(X,Tex(i,:),'or','linewidth',2);
    t=num2str(t);
    str={['t=', t]};
    text(0.55,T(i,(M+1)/2),str,'color','blue');
    t=double(t);
    title('FTCS Method');
    end
end
legend('Tnumeric','Texact');
N=num2str(M-1); xlabel('x'); ylabel('T'); title(['FTCS Method, N=', N]);

    case 2

%Crank-Nicolson method
clear ,clc
t=40; %time in seconds
dt=0.01; %time step
v=1/(pi^2);
M=101; %number of nodes
x0=0; xL=1; x=xL-x0; dx=x/(M-1);
N=t/dt+1; %number of nodes in time axis
T=zeros(N,M);
%initial conditions:
for j=1:M
    x=(j-1)*dx;
    T(N,j)=cos(pi*(x-0.5));
end
%boundary conditions:
for i=1:N
    T(i,1)=0;
    T(i,M)=0;
end
%y0=0; yN=l; F0=0; delta0=0;
a=-(2+2*(dx^2)/(v*dt)); b=1; c=1; d=zeros(N,M);
F=zeros(N,M); delta=zeros(N,M);
%TDMA
for j=N:-1:2
    for i=2:M-1 %forward elemination
        d(j,i)=(-2*(dx^2)/(v*dt))*T(j,i)-T(j,i+1)+2*T(j,i)-T(j,i-1);
        F(j-1,i)=-b/(a+c*F(j-1,i-1));
        delta(j-1,i)=(d(j,i)-c*delta(j-1,i-1))/(a+c*F(j-1,i-1));
    end
    for i=M:-1:3 %backward substitution
        T(j-1,i-1)=T(j-1,i)*F(j-1,i-1)+delta(j-1,i-1);
    end
end
%exact solutions:
Tex=zeros(N,M);
for i=0:M-1
   x=i*dx;
   for j=0:N-1
       tx=j*dt;
       Tex(N-j,i+1)=exp(-tx)*cos(pi*(x-0.5));
   end
end
figure ; X=0:dx:1;
for i=N:-1:1
    t=dt*(N-i);
    if t==0 || t==0.1 || t==0.2 || t==0.3 || t==0.4 || t==4 || t==40
    hold on
    plot(X,T(i,:),'color','black','linewidth',3);
    plot(X,Tex(i,:),'or','linewidth',1); t=num2str(t);
    str={['t=', t]};
    text(0.55,T(i,(M+1)/2),str,'color','blue'); t=double(t);
    end
end
legend('Tnumeric','Texact');
N=num2str(M-1); xlabel('x'); ylabel('T'); title(['Crank Nicholson Method, N=', N]);
end
