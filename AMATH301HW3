%%
% This is my homework for AMATH301 in 2015 WINTER. Current homework can be find in 
% http://courses.washington.edu/am301/
% They are not exactly the same but they are very similar. (Basically same question with different parameters)
%%

clear all;close all;clc
%problem1
%a 
X =(0:0.01:5)';
K=2;
n=3;
b0=10;
p=2;
bx=b0.*X.^n./(K^n+X.^n);
cx=p.*X;
fx=bx-cx;
A1=[bx cx fx];
save -ascii A1.dat A1

%b
[xopt,fval]=fminsearch('q1b',5)
A2=[xopt -fval]';
save -ascii A2.dat A2

%c


m=200; % number of generations
n=50; % number of trials
n2=10; % number of trials to be kept 

rng(0);
% initial (random) guesses
A=5+randn(n,1);
A3=[];
E=[];
for jgen=1:m
 for j=1:n % evaluate objective function
     X=A(j);
    E(j)= p.*X-b0.*X.^3./(K^3+X.^3);
 end
 [Es,Ej]=sort(E); % sort from small to large
 Ak1=A(Ej(1:n2)); % best 10 solutions

 Ak2=Ak1+randn(n2,1)/jgen; % 10 new mutations
 
 Ak3=Ak1+randn(n2,1)/jgen; % 10 new mutations
 
 Ak4=Ak1+randn(n2,1)/jgen; % 10 new mutations

 Ak5=Ak1+randn(n2,1)/jgen; % 10 new mutations
 A3(jgen,:)=[A(Ej(1)) -Es(1)];
 A=[Ak1; Ak2; Ak3; Ak4; Ak5];% group new 50
end
A3
save -ascii A3.dat A3

%d
rng(0);
[xopt1,fval]=ga(@q1b,1,[],[],[],[],0,[],[],[])
A4=[xopt1;-fval]
save -ascii A4.dat A4


%p2
%clear all; close all;clc
%a
A=[1000 600 301;10 5 3;4 2 3;2 1 1];
Aeq=[];
b=[2800000;25000;10000;6000];
beq=[];
A5=[A b];
save -ascii A5.dat A5

%b
p=[38;22;12];
Lb=[0;0;0];
A6=[p Lb];
save -ascii A6.dat A6

%c
[x,P]=linprog(-p,A,b,[],[],Lb,[]);
Nr=x(1);
Nd=x(2);
Nh=x(3);
A7=[round(x);-P];
save -ascii A7.dat A7

%p3
clear all; close all; clc;
N = [  7.24;  9.64; 12.87; 17.07; 23.19; 
            31.44; 38.56; 50.19; 62.98; 76.21; 
            92.23;106.02;123.20;132.16;151.33; 
           179.32;203.30;226.54;248.71;281.42; 
           307.75];
t = (-90:10:110)';
h = 10;
n = size(t,1);
A8 = zeros(n,1);
A8(1) = (-N(3) + 4*N(2) - 3*N(1))/(2*h);
A8(end) = (3*N(end) - 4*N(end-1) + N(end-2))/(2*h);
A8(2:end-1) = (N(3:end) - N(1:end-2))/(2*h);
save -ascii A8.dat A8

%p4
clear all; close all; clc;
%4.a
re = 0.308;
r0 = 0.478;
theta = 0.7051;
h = 0.017;
r=[.308 .325 .342 .359 .376 .393 .410 .427 .444 .461 .478];
T=[640 794 885 943 1034 1064 1114 1152 1204 1222 1239];
y1 = r.*T.*theta;
y2 = r.*theta;
num = (h/3).*sum(y1(1:2:end-2)+4.*y1(2:2:end-1)+y1(3:2:end));
denom = (h/3).*sum(y2(1:2:end-2)+4.*y2(2:2:end-1)+y2(3:2:end));
T1 = num/denom;

%4.b
T2 = trapz(r,y1)/trapz(r,y2);
A9 = [T1 T2].';
save -ascii A9.dat A9


%p5
m=10;
v0=10;
vt=5;
A10=quad(@(x) m./(-x.*sqrt(x)),10,5,1.e-4)
        
save -ascii A10.dat A10
        
