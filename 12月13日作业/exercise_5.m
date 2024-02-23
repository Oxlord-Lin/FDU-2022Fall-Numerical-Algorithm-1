%% LSQR
m = 1e3;
n = 60;
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
    
    % 然后进行LSQR过程
    Bk = B(1:k+1,1:k);
%     imagesc(Bk);
%     pause(1);
    v = norm(residual)*[1,zeros(1,k)]';
    % 对Bk和v同时做Givens旋转
    for j = 1:k
        alpha = Bk(j,j);
        beta = Bk(j+1,j);
        c = alpha/sqrt(alpha^2+beta^2);
        s = beta/sqrt(alpha^2+beta^2);
        G = eye(k+1);
        G(j,j) = c;
        G(j+1,j+1) = c;
        G(j,j+1) = s;
        G(j+1,j) = -s;
        Bk = G*Bk;
        v = G*v;
    end
    Rk = Bk(1:k,1:k);
    y = Rk\v(1:k,1);
    xk = x0 + V(1:n,1:k)*y;
    X(1:n,k+1) = xk;
    error(1,k+1) = norm(Bk*y-v);
    k = k+1;
end
norm(U'*A*V-B)
figure()
plot(error)
figure
plot(log10(abs(error-(norm(b-A*(A\b)))*ones(size(error)))))

