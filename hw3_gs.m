clc
clear all
h=pi./10;
k=pi./20;
gama=(h./k).^2;
y=transpose((pi./2):-k:0);       % y values %
x=0:h:pi;                        % x values %
u_0y=cos(y);                     % BC-1 u(0,y)=cos(y) %
u_piy=-cos(y);                   % BC-2 u(pi,y)=-cos(y) %
u_x0=cos(x);                     % BC-3 u(x,0)=cos(x) %
u_xpi2=linspace(0,0,length(x));  % BC-4 u(x,0)=0 %
u=ones(11).*0.5;                      % initial guess for all u is 0.5 %
u(1:11,1)=u_0y;
u(1:11,11)=u_piy;
u(1,:)=u_xpi2;
u(11,:)=u_x0;
% RHS of the governing equation %
for j=1:11
    for i=1:11
        f(j,i)=-1.*(cos(x(i)+y(j))+cos(x(i)-y(j)));
    end
end
for k=1:1000
    for j=2:10
        for i=2:10
    u(1:11,1,k+1)=u_0y;
    u(1:11,11,k+1)=u_piy;
    u(1,:,k+1)=u_xpi2;
    u(11,:,k+1)=u_x0;
    u(i,j,k+1)=(1/(2+2.*gama)).*(-1.*h.*h.*f(i,j)+u(i,j+1,k)+u(i,j-1,k+1)+gama.*(u(i+1,j,k)+u(i-1,j,k+1)));

        end
    end
    A(:,:,k)=abs(1-u(2:10,2:10,k)./u(2:10,2:10,k+1));
    B=max(max(A));
    if B>=0.0001
        continue
    else
        break
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%5
fig=figure('Name','Gauss Seidel Iteration');
[x,y] = meshgrid(x,y);
surf(x,y,u(:,:,k+1))
title('Gauss Seidel Iteration','FontSize',20)
xlabel('x axes','FontSize',16)
ylabel('y axes','FontSize',16)
zlabel('u values','FontSize',16)
c=colorbar;
c.Label.String = 'U values';
grid on
orient(fig,'landscape')
print('Gauss Seidel Iteration', '-dpdf','-fillpage');
