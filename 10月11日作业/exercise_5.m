%生成矩阵A,b1,B,b2过程
dimension=[10:10:1000];
forwarder=zeros(100,1);
res=zeros(100,1);
for k=1:100
    n=dimension(k);
v1=ones(n-1,1);
v2=ones(n-1,1).*6;
v3=ones(n-1,1).*8;
v4=ones(n,1).*6;
v5=ones(n,1).*8;
A1=diag(v5)+diag(v1,1)+diag(v2,-1);
B1=diag(v4)+diag(v1,1)+diag(v3,-1);
b1=zeros(n,1);
b2=b1;
b1(1,1)=9;
for i=1:n-1
    b1(i,1)=15;
end
b1(n,1)=14;
b2(1,1)=7;
for i=1:n-1
    b2(i,1)=15;
end
b2(n,1)=14;

%% 不用选主元的方法求解x，并计算其向前误差与残差|Ax-b1|
b=b1;
A=A1;
for i=1:n-1
    t=A(i+1,i)/A(i,i);
    A(i+1,i:i+1)=A(i+1,i:i+1)-t*A(i,i:i+1);
    b(i+1,1)=b(i+1,1)-t*b(i,1);
end
for i=n:-1:2
    b(i,1)=b(i,1)/A(i,i);
    A(i,i)=1;
    b(i-1,1)=b(i-1,1)-b(i,1)*A(i-1,i);
    A(i-1,i)=0;
end
x=b;
ferror=max(abs(x-ones(n,1)));%取无穷范数作为向前误差
r=max(abs(A1*x-b1));%取无穷范数作为残差
forwarder(k,1)=ferror;
res(k,1)=r;

%% 用选主元的方法求解x,并计算其向前误差与残差|Ax-b1|
% A=B1;
% b=b2;
% [l,u,p]=lu(A);
% y=l\(p*b);
% x=u\y;
% x=A\b;
% 
% ferror=max(abs(x-ones(n,1)));%取无穷范数作为向前误差
% r=max(abs(B1*x-b2));%取无穷范数作为残差
% forwarder(k,1)=ferror;
% res(k,1)=r;

end
figure
scatter(dimension,log(forwarder));%向前误差图表
figure
scatter(dimension,log(res));%残差图表