%% 生成一个tridiagonal矩阵
n=10;
% while 1
A=zeros(n,n);
A(1,1:2)=rand(1,2)+[2,0];%保证对角元最大
for i=2:n-1
    A(i,i-1:i+1)=rand(1,3)+[0,2,0];%保证对角元最大
end
A(n,n-1:n)=rand(1:2)+[0,10];%保证对角元最大
% if det(A)~=0 
%     break;%保证A是非奇异的，这样之后才能解出对应的b
% end
% 期末复习时的批注：由圆盘定理，对角占优阵一定是非奇异的，故这一步不需要！
% end
%再生成一个b用于测试
b0 = rand(n,1);
b = b0; %保留原始数据
A0 = A;
%% 然后从上到下进行Gaussian elimination
for i=2:n
    A(i,i) = A(i,i)-A(i-1,i).*A(i,i-1)./A(i-1,i-1);
    b(i,1) = b(i,1)-b(i-1,1).*A(i,i-1)./A(i-1,i-1);
end
%然后从下到上进行Gaussian elimination，解出x并储存在向量b中，这个过程中不需要再去改变矩阵A
for i=n:-1:2
    b(i,1) = b(i,1)./A(i,i);
    b(i-1,1) = b(i-1,1) - b(i,1).*A(i-1,i);
end
b(1,1) = b(1,1)./A(1,1);

b %打印出b
norm(b-A0\b0)