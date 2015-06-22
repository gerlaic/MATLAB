%%
% This is my homework for AMATH301 in 2015 WINTER. Current homework can be find in 
% http://courses.washington.edu/am301/
% They are not exactly the same but they are very similar. (Basically same question with different parameters)
%%

% init
close all; clear all; clc;

% Problem 1
A = [1 -6; -2 5];
B = [0 7; -7 0];
C = [1 -9 5; 2 0 -5];
D = [-1 0; 0 4; 3 0];
x = [1; 0];
y = [0; 1];
z = [-1; 2; -3];

A1 = A+B;
A2 = 2*x+3*y;
A3 = A*y;
A4 = B*(2*x-3*y);
A5 = D.'* z;
A6 = C.'*x + 3.*z;
A7 = C*D;
A8 = B*A;
A9 = A - D'*C';

A10 = D(1,:);
A11 = C(1,2);
A12 = D(2:3,2);
A13 = D.^2;

for i=1:13
    str = strcat('A', int2str(i), '.dat ');
    var = strcat('A', int2str(i));
    save(str,var,'-ascii');
end

%Problem 2
x1 = 200;
for i = 1:2000
    x1 = (x1-0.1);
end
x1 = abs(x1);

x2 = 200;
for i = 1:1600
    x2 = (x2-0.125);
end
x2 = abs(x2);

x3 = 200;
for i = 1:1000
    x3 = (x3-0.2);
end
x3= abs(x3);

x4 = 200;
for i = 1:800
    x4 = x4-0.25;
end
x4 = abs(x4);


x5 = 200;
for i = 1:500
    x5 = x5-0.4;
end
x5=abs(x5);

x6 = 200;
for i = 1:400
    x6 = x6-0.5;
end
x6 = abs(x6);

A14 = [x1 x2 x3 x4 x5 x6];
save -ascii A14.dat A14;

%Problem 3
v_0 = 36;
t = 4;
g = 9.8;
c_d = 0.25;

mr = 200;
ml = 100;
A15 = 0:1:50;
for i = 1:100
    mc = (ml+mr)/2;
    fc = sqrt(mc .*g./c_d).*tanh(t.*sqrt(c_d.*g./mc)) - v_0;
    A15(1,i) = mc;
    if fc > 0
        mr = mc;
    else 
        ml = mc;
    end
    
    if abs(fc)< 10^(-4)
        break
    end
end

A15 = A15(1,1:i);
fc = fzero(@(mc) sqrt(mc .*g./c_d).*tanh(t.*sqrt(c_d.*g./mc)) - v_0, 100, optimset('TolX', 1e-6));
A15(1,i+1) = fc;
A15 = A15.';
save -ascii A15.dat A15;

%Problem 4
R1 = 20;
R2 = 5;
R3 = 30;
R4 = 25;
R5 = 10;
R6 = 15;
V2 = 0;
V3 = 50;

A = [R6+R1+R2 -R1 -R2; -R1 R3+R4+R1 -R4; -R2 -R4 R5+R2+R4];
[L,U,P] = lu(A);
save -ascii A16.dat A P L U;

A17=[];
for V1=49:2:99
    z = L\(P*[V1;0;50]);
    x = U\z;
    A17 = [A17 x];
end
save -ascii A17.dat A17;

i = 1;
for V1=49:2:99
    A18(:,i) =inv(A) * [V1;0;50];
    i=i+1;
end
A18 = A18-A17;
save -ascii A18.dat A18;
