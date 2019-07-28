%%to find the solution for
%%30x10 make the c=1
%%90x30 make the c=3
%%180x60 make the c=6
clear, clc
b=10; a=3*10; c=6;
N=b/(1/c)+1; M=a/(1/c)+1;
dx=a/(M-1); dt=(dx^2)/4;
T=18*ones(N,M); T(N,:)=150;
T(N,1)=(150+18)/2; T(N,M)=(150+18)/2;
n=1; w=1; T2=T;
err=1; tol=10^-2;
while err>tol %five point formula
    t=dt*(n-1); 
    T3(1:N,:)=T(N:-1:1,:);
    if t<=5
    if rem(t,1)==0
        figure
        contour(T3,1000);
        colormap(flipud(hot))
        t=num2str(t);
        title(['t=' t 's']); xlabel(['a','(', num2str(M-1) ,')']); 
        ylabel(['b','(', num2str(N-1) ,')']);     
        t=str2double(t);
    end
    else 
        if rem(t,10)==0
        figure
        contour(T3,1000);
        colormap(flipud(hot))
        t=num2str(t);
        title(['t=' t 's']); xlabel(['a','(', num2str(M-1) ,')']); 
        ylabel(['b','(', num2str(N-1) ,')']);     
        t=str2double(t);
        end
    end
    n=n+1;
    for j=N-1:-1:2
        for i=2:M-1
            T(j,i)=0.25*(T(j,i+1)+T(j,i-1)+T(j+1,i)+T(j-1,i));
        end
    end
    Tm=T((N+1)/2,(M+1)/2);
    if Tm>=30
        while w==1
        tm=t;
        disp(t);
        disp('midpoint of the plate reached 30 C.');
        w=2;
        end
    end
    err=max(max(T-T2)); T2=T;
end
disp(t)
disp('Result converged in sec')

