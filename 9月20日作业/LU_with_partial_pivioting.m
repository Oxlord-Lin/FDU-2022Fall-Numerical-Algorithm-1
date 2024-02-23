%% 复习时发现未做选主元的LU分解，故补充如下
% 以长方形矩阵作为例子（m<n）
m = 5;
n = 5;
A = randn(m,n);
A0 = A;
P = eye(m);
for i = 1:m-1
    [M,index] = max(abs(A(i:m,i)));
    index = index + i - 1;
    P_temp = eye(m); 
    P_temp(i,i) = 0; P_temp(index,index) = 0; P_temp(i,index) = 1; P_temp(index,i) = 1;
    P = P_temp * P; % 记录置换阵
    v = A(i,:);
    A(i,:) = A(index,:);
    A(index,:) = v; % 交换行，注意是整行交换，而不是在Schur补内交换
    A(i+1:m,i) = A(i+1:m,i)./A(i,i);
    A(i+1:m,i+1:n) = A(i+1:m,i+1:n) - A(i+1:m,i)*A(i,i+1:n); % 修正Schur补
end

L = zeros(m,m);
U = zeros(m,n);
for i = 1:m
    for j = 1:n
        if i>j
            L(i,j) = A(i,j);
        else 
            U(i,j) = A(i,j);
        end
    end
end
L = L + eye(m);


