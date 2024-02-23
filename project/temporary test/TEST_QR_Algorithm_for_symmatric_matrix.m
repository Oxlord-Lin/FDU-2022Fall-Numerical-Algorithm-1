%% 对称矩阵的QR算法(辅助自己理解单步位移迭代，toy algorithm)
% 此部分属于算法主体结构
n = 4; %n为矩阵阶数
% 生成对称矩阵A
A = zeros(n);
for i = 1:n
    for j = 1:i
        A(j,i) = rand();
        A(i,j) = A(j,i);
    end
end
%% 三对角化
[A,U] = upper_hessenberg(A); %对A进行上海森伯格化，由于A是对称的，故实际上为三对角化

%% test
%A = diag([1:5]);

%% 带Wilkinson位移的隐式QR算法
x = A(1,1)-Wilkinson_shift(A);
y = A(2,1);
for k = 1:n-1
    %生成Givens旋转中的角度
    temp = sqrt(x^2+y^2);
    c = x/temp;
    s = y/temp;
    %生成Givens旋转矩阵
    G = eye(n); G(k,k) = c; G(k+1,k+1) = c; G(k,k+1) = -s; G(k+1,k) = s;
    %驱赶"气泡"
    A = G'*A*G
    %为下一次迭代的Givens矩阵所需参数做准备
    if k<n-1
    x = A(k+1,k);
    y = A(k+2,k); 
    end
end



