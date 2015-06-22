function fx=q1b(x0)
X=x0;
K=2;
n=3;
b0=10;
p=2;
bx=b0*X^n/(K^n+X^n);
cx=p*X;
fx=-(bx-cx);

