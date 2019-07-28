clear;
clc;
var=5;
n=20;
k=20;
zt=2/n;
et=1/k;
ze=0:zt:2;
e=0:et:1;
for i=1:length(ze)
    x(i,1)=ze(i);
    if x(i,1)>1
        y(i,1)=0.5*x(i,1)-0.5;
    else
        y(i,1)=0;
    end
    x(i,length(e)) = (i-1)*et*0.5;
    y(i,length(e)) = 1;
end
for j = 1:length(e)
    if j > 0.75*(length(e)-1)
        x(length(e), j) = 0.5;
    else
        x(length(e),j) = 2 - (j-1)*(2/(length(e)-1));
    end
    if 1.5<x(length(e),j) && x(length(e),j)<=2;
        y(length(e),j) = -1*x(length(e),j) + 2.5;
    end
    if 1<x(length(e),j) && x(length(e),j)<=1.5;
        y(length(e),j) = x(length(e),j) - 0.5;
    end
    if 0.50<x(length(e),j) && x(length(e),j)<=1;
        y(length(e),j) = 0.5;
    end
    if x(length(e),j) == 0.5
        y(length(e),j) = 0.5 + (j-0.75*(length(e)-1)-1)*zt;
    end
    x(1,j) = 0;
    y(1,j) = (j-1)*et;
end
for l=1:2000
    for i=2:n
        for j=2:k
            a=((x(i+1,j)-x(i-1,j))/(2*et))^2+((y(i+1,j)-y(i-1,j))/(2*et))^2;
            b=((x(i+1,j)-x(i-1,j))/(2*et))*((x(i,j+1)-x(i,j-1))/(2*zt))+((y(i+1,j)-y(i-1,j))/(2*et))*((y(i,j+1)-y(i,j-1))/(2*zt));
            l=((x(i,j+1)-x(i,j-1))/(2*zt))^2+((y(i,j+1)-y(i,j-1))/(2*zt))^2;
            x(i,j)=((a/zt^2)*(x(i,j+1)+x(i,j-1))+(l/et^2)*(x(i+1,j)+x(i-1,j))-b/(2*zt*et)*(x(i+1,j+1)-x(i-1,j+1)+x(i-1,j-1)-x(i+1,j-1)))/(2*(a/zt^2+l/et^2));
            y(i,j)=((a/zt^2)*(y(i,j+1)+y(i,j-1))+(l/et^2)*(y(i+1,j)+y(i-1,j))-b/(2*zt*et)*(y(i+1,j+1)-y(i-1,j+1)+y(i-1,j-1)-y(i+1,j-1)))/(2*(a/zt^2+l/et^2));
        end
    end
end
hold on
for i=1:n+1
    plot(x(i,1:n+1),y(i,1:n+1))
end
for j=1:k+1
    plot(x(1:k+1,j),y(1:k+1,j))
    title('Mesh')
    xlabel('ze')
    ylabel('e')
end
vel=zeros(length(ze),length(e));
vel(1,:)=1;
vel(:,length(e))=1;
vel(:,1)=1;
for i=1:((length(e)-1)*0.25)+1
    vel(length(ze),i)=1;
end
for l=1:2000
    for t=2:((length(e)-1)*0.25)+1
        vel(length(ze),t)=0.25*(vel(length(ze),t-1)+vel(length(ze),t+1)+2*vel(length(ze)-1,t));
    end
    for t=2:length(ze)-1
        vel(t,length(ze))=0.25*(vel(length(ze),t-1)+vel(length(ze),t+1)+2*vel(length(ze)-1,t));
    end
    for i=2:n
        for j=2:k
            vel(i,j)=0.25*(vel(i-1,j)+vel(i+1,j)+vel(i,j+1)+vel(i,j-1)+zt*zt*var);
        end
    end
end
figure
[C,h]=contourf(x,y,vel(1:n+1,1:k+1));
clabel(C,h)
title('Velocity')
xlabel('x')
ylabel('y')

