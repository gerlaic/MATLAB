%%
% This is my homework for AMATH301 in 2015 WINTER. Current homework can be find in 
% http://courses.washington.edu/am301/
% They are not exactly the same but they are very similar. (Basically same question with different parameters)
%%

close all; clear all;clc;
%problem 1 a 
A1 = [];
n = 2;
i = 1;
while(cond(hilb(n),inf)< 10^16&&n < 2000)
    A1(i,:)= cond(hilb(n),inf);
    i = i + 1;
    n = n + 1;
end
save -ascii A1.dat A1

%b
A2 = [];
for i=2:1:(n-1)
    z=ones(i,1);
    h=hilb(i);
    A2(i-1,1)=norm(z-h\(h*z),Inf);
    A2(i-1,2)=norm(z-h*(h\z),Inf);
end
save -ascii A2.dat A2


%c
A3 = [];
for i=2:1:(n-1)
    z=ones(i,1);
    h=hilb(i);
    ih=invhilb(i);
    A3(i-1,1)=norm(z-ih*(h*z),Inf);
    A3(i-1,2)=norm(z-h*(ih*z),Inf);
end
A3
save -ascii A3.dat A3

%problem2 a
P=[4 -1 -1 0 0 0 0 0 0 0;
    -1 5 -1 -1 -1 0 0 0 0 0;  
    -1 -1 5 0 -1 -1 0 0 0 0;
    0 -1 0 5 -1 0 -1 -1 0 0;
    0 -1 -1 -1 6 -1 0 -1 -1 0;
    0 0 -1 0 -1 5 0 0 -1 -1;
    0 0 0 -1 0 0 4 -1 0 0;
    0 0 0 -1 -1 0 -1 5 -1 0;
    0 0 0 0 -1 -1 0 -1 5 -1;
    0 0 0 0 0 -1 0 0 -1 4];
b=[0 0 0 0 0 0 1 1 1 1].';
A4=P\b
save -ascii A4.dat A4;

%b

x0=[1; zeros(9,1)];
A5=[x0];
D=diag(diag(P));
T=P-D;
for i=1:1:100
    A5(:,i+1)=D\(b-T*A5(:,i));
    if norm((A5(:,i+1)-A5(:,i)),Inf)<=(1e-4)
        break
    end 
end
save -ascii A5.dat A5


%b
x0=[1; zeros(9,1)];
A6=[x0];
S = tril(P);
T=P-S;
for i=1:1:100
    A6(:,i+1)=S\(b-T*A6(:,i));
    if norm((A6(:,i+1)-A6(:,i)),Inf)<=(1e-4)
        break
    end 
end
save -ascii A6.dat A6

%c
[A7,a1,a2,n]=bicg(P,b,1e-4)
A7=[A7;n]
save -ascii A7.dat A7


%c
M =[0 1/4 1/4 0 0 0 0 0 0 0;
    1/2 0 1/4 1/4 1/6 0 0 0 0 0;
1/2 1/4 0 0 1/6 1/4 0 0 0 0;
0 1/4 0 0 1/6 0 1/2 1/4 0 0;
0 1/4 1/4 1/4 0 1/4 0 1/4 1/4 0;
0 0 1/4 0 1/6 0 0 0 1/4 1/2;
0 0 0 1/4 0 0 0 1/4 0 0;
0 0 0 1/4 1/6 0 1/2 0 1/4 0;
0 0 0 0 1/6 1/4 0 1/4 0 1/2;
0 0 0 0 0 1/4 0 0 1/4 0];
sum(M,1);
[V,D]=eig(M)
A8=V(:,1);
save -ascii A8.dat A8

%problem3
%a
N = [ 7.24; 9.64; 12.87; 17.07; 23.19; ...
 31.44; 38.56; 50.19; 62.98; 76.21; ...
 92.23;106.02;123.20;132.16;151.33; ...
 179.32;203.30;226.54;248.71;281.42; ...
 307.75];
x=(-90:10:110)';
y=log(500./N -1)
A9=polyfit(x,y,1)
save -ascii A9.dat A9

%a ii
A=[x ones(21,1)];
[X,a1,a2,n]=lsqr(A,y);
X

A10=[A9-X' n]
save -ascii A10.dat A10

%iii
x1=120:10:300;
Y = polyval(A9,x1')
A11 = 500./(1+exp(Y))
save A11.dat A11 -ascii

%b i
t=(-90:10:110)';
x0=[0.01;1;500];
f=@(x,t) (x(3))./(1+x(2).*exp(-x(1).*t));
A12 = (lsqcurvefit(@(x,t) f(x,t),x0,t,N))'
save -ascii A12.dat A12

%ii
t2=(120:10:300)';
A13=f(A12,t2)
save -ascii A13.dat A13



%problem 4 a
pcoef=polyfit(t,N,length(t));
tp=(-90:1:110)';
Np=polyval(pcoef,tp);
plot(tp,Np)
save -ascii A14.dat Np

%b
Nint=interp1(t,N,tp)
save -ascii A15.dat Nint

%c
A16=spline(t,N,tp)
save -ascii A16.dat A16



