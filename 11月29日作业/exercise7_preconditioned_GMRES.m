%% GMRES Solver for nonsymmetric coefficient matrix
n = 250;
residuals = zeros(n-1,1); %用于记录残差的变化情况
L0 = diag(randn(1,n)) + diag(randn(1,n-1),-1);
U0 = diag(randn(1,n)) + diag(randn(1,n-1),+1);
M = L0*U0;
L = L0 + (1e-3)*randn(n,1)*randn(1,n);
L = tril(L);
U = U0 + (1e-3)*randn(n,1)*randn(1,n);
U = triu(U);
A = L*U;
b = randn(n,1);
%%
% % 这样使用预条件其实不太对，M应该体现在每一次迭代中；时间不够用了，临时这样处理
A0 = A;
A = A/M;
%%
x_real = A0\b; %真实的解x_real
x_ori = randn(n,1); %初始向量
r = b-A*x_ori;%初始残差
residuals(1,1) = norm(r);
Q = zeros(n,n); %用于存储Krylov子空间的正交基
H = zeros(n,n); %用于存储A*q_k的在q_1至q_k+1的投影
Q(1:n,1) = r./norm(r);
flg = 0; %用于标记是否因H(i+1,i)足够小而提前终止
for i = 1:n-1 %由于用于表示A^-1的A的多项式次数不会超过n-1，故外循环次数最多n-1次
    y = A*Q(1:n,i); %y=A*q_i
    for j = 1:i
        H(j,i) = Q(1:n,j)'*y; %向q_1至q_i进行投影
        y = y - H(j,i)*Q(1:n,j); %减去相应的分量
    end
    H(i+1,i) = norm(y); 
    if abs(H(i+1,i)) < 1e-16 %若H(i+1,i)足够小
        flg = 1;
        disp("提前终止，因为H(i,i+1)足够小，认为A^-1*r已经完全落在i维Krylov子空间中")
        break
    end
    Q(1:n,i+1) = y./H(i+1,i); %生成第i+1个正交基向量
    r_2 = zeros(i+1,1);
    r_2(1,1) = norm(r);
    [c,res] = least_square(H(1:i+1,1:i),r_2); %最小二乘问题的局部函数（详后）
    %注意：尽管此步骤已经把q_i+1生成出来，但此处的c是对q_1到q_i的线性组合
    %如果没有提前终止Arnoldi过程，对于最后一次循环，c是对q_1到q_n-1的线性组合，即仍然不是精确解
    %精确解应该是x_ori加上q_1到q_n的某个线性组合
    residuals(i,1) = res; %res是残差
end
if flg ==0 
    figure
    plot([1:n-1],log10(residuals));
    xlabel("Arnoldi过程的迭代次数");
    ylabel("残差，取对数");
    title("非对称矩阵在不带GMRES过程中残差的变化情况(阶数n=250)")
     %这里补充一个步骤，通过对Q的各列进行线性组合得到x，并与通过x=A\b得到的x_real进行比较
    y= A*Q(1:n,n);
    for k = 1:n
    H(k,n) = Q(1:n,k)'*y; %向q_1至q_n进行投影
    end
    r_2 = zeros(n,1);
    r_2(1,1) = norm(r);
    c = H\r_2;
    x_add = Q*c;
    x = x_ori + x_add;
    er = norm(M\x-x_real)
end
if flg==1
    %如果提前结束，那么q_1至q_i的某个线性组合就能得到x_add
    figure()
    plot([1:i-1],log10(residuals(1:i-1,1)));
    xlabel("迭代次数");
    ylabel("残差");
    title("非对称矩阵在GMRES过程中残差的变化情况(阶数n=500)")
end

%%
% 用于处理最小二乘问题并返回最优解与残差
function [x,res] = least_square(A,b)
x = pinv(A)*b;
res = norm(A*x-b);
end

