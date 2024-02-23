%% 生成两个矩阵，其中A的条件数10^6,B的条件数10^16
m = 100;
n = 50;
U = orth(randn(m)); %U = eye(m);
V = orth(randn(n)); %V = eye(n);
% % 奇异值服从几何分布
% sig1 = linspace(1,6,n);
% sig2 = linspace(1,16,n);
% for i = 1:n
%     sig1(1,i) = 10^sig1(1,i);
%     sig2(1,i) = 10^sig2(1,i);
% end

% % 奇异值服从均匀分布
sig1 = linspace(1,1e6,n);
sig2 = linspace(1,1e16,n);
b1 = U*[sig1';ones(m-n,1)]; %右侧向量1
b2 = U*[sig2';ones(m-n,1)]; %右侧向量2
x_star = V*[ones(n,1)]; %精确解 
sig1 = [diag(sig1);zeros(m-n,n)]; %奇异值1
sig2 = [diag(sig2);zeros(m-n,n)]; %奇异值2
A1 = U*sig1*V'; %矩阵1
A2 = U*sig2*V'; %矩阵2

x1_CGS = CGS(A1,b1);
x1_MGS = MGS(A1,b1);
x1_Hou = Hou(A1,b1);
error_A1 = log10([norm(x1_CGS-x_star),norm(x1_MGS-x_star),norm(x1_Hou-x_star)]);
x2_CGS = CGS(A2,b2);
x2_MGS = MGS(A2,b2);
x2_Hou = Hou(A2,b2);
error_A2 = log10([norm(x2_CGS-x_star),norm(x2_MGS-x_star),norm(x2_Hou-x_star)]);
X = categorical({'CGS','MGS','Householder'});
X = reordercats(X,{'CGS','MGS','Householder'});
Y = [error_A1',error_A2'];
figure
bar(X,Y)
legend('条件数10^6','条件数10^{16}');
ylabel('log10(norm(x-x_*))');
% title('奇异值服从几何分布');
title('奇异值服从均匀分布');



%% 用经典GS求解LS问题
function x = CGS(A,b)
[m,n] = size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j = 1:n
    y=A(1:m,j);
    for i = 1:j-1
        R(i,j)=Q(1:m,i)'*A(1:m,j);
        y=y-Q(1:m,i).*R(i,j);%经典的G-S方法
    end
    R(j,j)=norm(y);
    Q(1:m,j)=y/R(j,j);
end
b2 = Q'*b;
x = U_sol(R,b2);
end

%% 改进GS，不进行重正交化
function x = MGS(A,b)
[m,n] = size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j = 1:n
    y=A(1:m,j);
    for i = 1:j-1
        R(i,j)=Q(1:m,i)'*y;%改进的G-S方法
        y=y-Q(1:m,i).*R(i,j);
    end
    R(j,j)=norm(y);
    Q(1:m,j)=y/R(j,j);
end
b2 = Q'*b;
x = U_sol(R,b2);
end

%% 使用householder方法进行QR分解,进而求解法方程
function x = Hou(A,b)
[m,n] = size(A);
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