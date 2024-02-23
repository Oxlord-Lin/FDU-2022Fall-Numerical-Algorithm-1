%% Paige-Saunders bidiagonalization
m = 1e3;
n = 50;
A = randn(m,n); % 最小二乘问题中，m>n
B = zeros(m,n);
main_diag = zeros(1,n);
sub_diag = zeros(1,n); % 第一个元素不在此对角线上，但为了编程方便加上它
x0 = zeros(n,1);
b = randn(m,1);
residual = b - A*x0;
uc = residual./norm(residual); %uc是U的第一列
U = zeros(m);
U(1:m,1) = uc;
p = uc;
sub_diag(1,1) = 1;
V = zeros(n);
X=zeros(n,2); %用于储存产生的向量
X(1:n,1) = x0;
error = zeros(1,n); % 用于储存残差
error(1,1) = norm(residual);
k = 1;

while(sub_diag(1,k)>1e-16 && k<=n+1)
    % 首先进行二对角化
    U(1:m,k) = p./sub_diag(1,k);
    if k==n+1
        break;
    end
    if k==1
        r = A'*U(1:m,k);
    else
    r = A'*U(1:m,k)- sub_diag(1,k)*V(1:n,k-1);
    end
    main_diag(1,k) = norm(r);
    V(1:n,k) = r./main_diag(1,k);
    p = A*V(1:n,k) - main_diag(1,k)*U(1:m,k);
    sub_diag(1,k+1) = norm(p);
    B(k,k) = main_diag(1,k);
    B(k+1,k) = sub_diag(1,k+1);
    k = k+1;
end
norm(U'*A*V-B)


