clear all
clc

delta_x = 0.2;
delta_y = 0.2;
N=(1-0)/delta_x; %horizontal axis
M=(2-0)/delta_y; %vertical axis
x=0:delta_x:1;
y=0:delta_y:2;
Gama=(delta_x/delta_y)^2;
% Re = 1 5 10
Re=1 ; %Reynolds number

u(1,1:M+1)=0;    %Boundary condition
u(1:N+1,1)=0;    %Boundary condition
u(N+1,1:M+1)=0;  %Boundary condition
u(1:N+1,M+1)=1;  %Boundary condition
u(2:N,2:M)=0.5;  %Initial guess for internal points

w=0.5; %under relaxation parameter
epsilon=1;
k=0; %iteration counter
while epsilon>0.0001
    k=k+1
    for i=2:N
        for j=2:M
            A1=0.5*Re*delta_x;
            A2=2*(1+Gama);
            B1=1-A1*u(i,j);
            B3=2-B1;
            B2=Gama;
            B4=Gama;
            u_new(i,j)=(u(i+1,j)*B1+u(i-1,j)*B3+u(i,j+1)*B2+u(i,j-1)*B4)/A2;  % 2-D Burger equation expanded with
                                                                              % Finite Difference Method 
            u_new(i,j)=(1-w)*u(i,j) + w*u_new(i,j);
            epsilon = (u_new(i,j)-u( i,j))./u(i,j);  %Convergence test
            u(i,j)=u_new(i,j);
        end
    end
end
% Title name should be changed according to Re number
u=[u'];
[X,Y] = meshgrid(x,y);
[C,h] = contour(X,Y,u,15);
xlabel('x'), ylabel('y'), title('Velocity Seperation on The Wall ( u\infty = 1 m/s and Re = 1 )');
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)


