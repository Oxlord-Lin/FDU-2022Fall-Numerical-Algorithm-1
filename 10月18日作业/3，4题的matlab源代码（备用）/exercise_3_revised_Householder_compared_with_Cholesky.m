m = 100; %矩阵大小
n = 50; 
start = 0;
ending = 14;
step = 0.01;
f_error = zeros((ending-start)/step+1,1);
f_error2 = f_error; %Matlab使用深拷贝
count = 1;
for k = start:step:ending
    ex = 10^k;
    singular = 1 + (ex-1)*rand(n,1);
    singular(1,1) = ex;
    singular(n,1) = 1;
    V = orth(randn(n));
    U = orth(randn(m));
    A = U*[diag(singular);zeros(m-n,n)]*V';
    b = randn(m,1);
    standard_sol = pinv(A)*b; %用最小二乘解，精度最高
    x1 = sol_Hou(A,b); %Householder
    x2 = sol_Cho(A,b); %Cholesky
    f_error(count,1) = norm(x1-standard_sol); %向前误差，Householder
    f_error2(count,1) = norm(x2-standard_sol);%向前误差，Cholesky
    count = count + 1;
end
% 作图版本一
plot([start:step:ending],log10(f_error),'--');
hold on
plot([start:step:ending],log10(f_error2));
xlabel('log10(kapa)');
ylabel('向前误差(对10取对数)');
% title('用Householder作QR从而求解法方程时向前误差与条件数之间的关系（取对数）');
title('通过Cholesky与Householder求解法方程时的精确度随矩阵条件数的变化');

% % 作图版本二
% figure
% scatter([start:step:ending],f_error);
% title('用Householder作QR从而求解法方程时向前误差与条件数之间的关系');
% xlabel('log10(kapa)');
% ylabel('向前误差');

%% 使用householder方法进行QR分解,进而求解法方程
function x = sol_Hou(A,b)
temp = size(A);
m = temp(1,1);
n = temp(1,2);
for j =1:n
    v=A(j:m,j);
    len=norm(v);
    if real(v(1,1))>=0 %此处相当于对v与一个长度为norm(v)且只有第一个分量不为零的向量作差
        v(1,1)=v(1,1)+len; %避免舍入误差
    else
        v(1,1)=v(1,1)-len;%避免舍入误差
    end
    H=eye(m-j+1)-2*(v*v')./(v'*v); %基本上找到处理第j列的Householder反射子，但还要调整相位，以保证相乘为实数
    A(j:m,j:n)=H*A(j:m,j:n); %Householder反射子只需要作用于schur补即可
    b(j:m,1) = H*b(j:m,1); %Householder还要同时作用于右侧的向量b
end
x = U_sol(A(1:n,1:n),b(1:n,1));
end

%% 通过Cholesky分解求法方程
function x = sol_Cho(A,b)
R = Cholesky(A'*A);
y = L_sol(R',A'*b);
x = U_sol(R,y);
end

%% 计算Cholesky因子
function R = Cholesky(S)
n = size(S);
n = n(1,1);
for i = 1:n
    for j = i+1:n
        S(j,j:n) = S(j,j:n) - S(i,j:n)*S(i,j)'/S(i,i); % 此处勿忘转置S（i，j）
    end
    S(i,i:n) = S(i,i:n)./sqrt(S(i,i));
end
R = triu(S);
end

%% 求解下三角线性方程组
function x = L_sol(L,b)
n = size(L);
n = n(1,1);
for i = 1:n-1
    b(i,1) = b(i,1)/L(i,i);
    for j = i+1:n
        b(j,1) = b(j,1) - L(j,i)*b(i,1);
    end
end
b(n,1) = b(n,1)/L(n,n);
x = b;
end

%% 求解上三角线性方程组
function x = U_sol(U,b)
n = size(U);
n = n(1,1);
for i = n:-1:2
    b(i,1) = b(i,1)/U(i,i);
    for j = i-1:-1:1
        b(j,1) = b(j,1) - U(j,i)*b(i,1);
    end
end
b(1,1) = b(1,1)/U(1,1);
x = b;
end

