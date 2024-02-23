m = 100; %矩阵大小
n = 50; 
start = 0;
ending = 14;
step = 0.01;
f_error = zeros((ending-start)/step+1,1);
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
    x = sol_Cho(A,b);
    standard_sol = pinv(A)*b; %用最小二乘解，精度最高
    f_error(count,1) = norm(x-standard_sol); %向前误差
    count = count + 1;
end
% 作图版本一
scatter([start:step:ending],log10(f_error));
title('用Cholesky分解求解法方程时向前误差与条件数之间的关系（取对数）');
xlabel('log10(kapa)');
ylabel('向前误差(对10取对数)');

% 作图版本二
figure
scatter([start:step:ending],f_error);
title('用Cholesky分解求解法方程时向前误差与条件数之间的关系');
xlabel('log10(kapa)');
ylabel('向前误差');

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

%% 通过Cholesky分解求法方程
function x = sol_Cho(A,b)
R = Cholesky(A'*A);
y = L_sol(R',A'*b);
x = U_sol(R,y);
end
