m=4;
n=4;
% %方法一：采用partial pivoting，用matlab自带的LU分解即可实现
A=randn(m,n);
% [L1,U1,P1]=lu(A);

%方法二：直接进行LU分解，用自己写的程序完成
last=min(m,n);
for j=1:last
        A(j+1:m,j)=A(j+1:m,j)./A(j,j);%把第j列变成各行之间的"系数之比"
        A(j+1:m,j+1:n)=A(j+1:m,j+1:n)-A(j+1:m,j)*A(j,j+1:n);%对Schur complement进行消元操作 
end
v=ones(max(m,n));
%需要保证自己写的L2和U2矩阵的大小和Matlab给出的L1与U1的大小是一致的，方便后面计算F范数
if n>=m
    U2=zeros(m,n);
    L2=zeros(m,m);
end
if n<=m
    U2=zeros(n,n);
    L2=zeros(m,n);
end

for i=1:m
    for j=1:n
        if j>=i
            U2(i,j)=A(i,j);
        end
        if j<i
            L2(i,j)=A(i,j);
        end
        if i==j
            L2(i,j)=1;%L的对角线上的元素设置为1
        end
    end
end
U2
L2

A
% norm(L1-L2)%计算L1-L2的F范数
% norm(U1-U2)%计算U1-U2的F范数
