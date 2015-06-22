clear all;close all; clc
%1
%a
dt=0.1;
f=@(t,y)(2-2*y*t)/(t.^2+1);
y(1)=1;
t(1)=0;
T=1;
y(2)=y(1)+dt*f(t(1),y(1));
t(2)=t(1)+dt;
for i=1:T/dt
    t(i+1)=t(i)+dt;
    y(i+1)=y(i)+dt*f(t(i),y(i));
end
A1=y';
save -ascii A1.dat A1

%b
f=@(t,y)(2-2*y*t)/(t.^2+1) ;
h=0.1
for j=1:10
k1=f(t(j) ,y(j) );
k2=f(t(j)+h/2,y(j)+h/2*k1);
k3=f(t(j)+h/2,y(j)+h/2*k2);
k4=f(t(j)+h ,y(j)+h *k3);
y(j+1)=y(j)+h/6*(k1+2*k2+2*k3+k4);
t(j+1)=t(j)+h;
end 
A2=y';
save -ascii A2.dat A2

%c
[t45,y]=ode45(f,[0:1],1);
A3=[t45,y];
save -ascii A3.dat A3


%2a
clear all;close all;clc
r1=2;
r2=5;
v1=110;
y(:,1)=v1;
y(:,2)=-100;
f=@(r,y)[y(2);-2/r*y(2)];
[t,y]=ode45(f,[r1 r2],[v1 -100]);
A4 =[t,y(:,1)];
save -ascii A4.dat A4

%b
r1=2;
r2=5;
v1=110;
y(:,1)=v1;
y(:,2)=-90;
f=@(r,y)[y(2);-2/r*y(2)];
[t,y]=ode45(f,[r1 r2],[v1 -90]);
A5=[t y(:,1)];
save -ascii A5.dat A5

%c
A6(:,1)=(-100:1:-90)';
for i=-100:1:-90
    y(:,1)=v1;
    y(:,2)=i;
    f=@(r,y)[y(2);-2/r*y(2)];
    [t,y]=ode45(f,[r1 r2],[v1 i]);
    A6(101+i,2)=y(end,1);
end
save -ascii A6.dat A6

%d 
x=A6(:,2);
y=A6(:,1);
y0=interp1(x,y,0)
save -ascii A7.dat y0

%problem3
%a
clear all; close all;clc
lamda=-41;dt=0.05;
ya(1)=1;ta(1)=0;
f=@(t,y) lamda*(y-t)+1;

for i=1:1/0.05
    ya(i+1)=ya(i)+dt*f(ta(i),ya(i));
    ta(i+1)=dt+ta(i);
end
A8=ya';
save -ascii A8.dat A8


%b
f=@(t,y) lamda*(y-t)+1;
yb(1)=1;
tb(1)=0;
h=0.05;
for j=1:1/0.05
    tb(j+1)=tb(j)+h;
    yb(j+1)=fzero(@(z)z-yb(j)-h*f(tb(j+1),z),yb(j));
    
end 
A9=yb';
save -ascii A9.dat A9;


%c
[t,y]=ode23tb(f,[0 1],[1]);
A10=[t,y];
save -ascii A10.dat A10


%problem 4
%a
clear all;close all;clc
sigma=10;
beta=8/3;
rho=28;
h=0.01;
lorenz=@(t,Y) [sigma*(Y(2)-Y(1)); Y(1)*(rho-Y(3))-Y(2);Y(1)*Y(2)-beta*Y(3)]
[t,Y]=ode45(lorenz,[0:0.01:25],[2 3 14]);
%plot3(Y(:,1),Y(:,2),Y(:,3)); view(2);
A11=[t,Y];
save -ascii A11.dat A11

%b
[t,Yb]=ode45(lorenz,[0:0.01:25],[2 3 (14+10^-9)]);
A12=[t,Yb];
save -ascii A12.dat A12

%c
deltat=[];
for i=1:length(Y(:,1))
    deltat(i,:)=sqrt((Y(i,1)-Yb(i,1))^2+(Y(i,2)-Yb(i,2))^2+(Y(i,3)-Yb(i,3))^2);
    lndeltat(i,:)=log(deltat(i,:));
end

p=polyfit((0:0.01:25)',lndeltat,1)
save -ascii A13.dat p

