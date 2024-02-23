f_er_A_no_piv = zeros(100,1);
res_A_no_piv = zeros(100,1);

f_er_B_piv = zeros(100,1);
res_B_piv = zeros(100,1);

f_er_B_no_piv = zeros(100,1);
res_B_no_piv = zeros(100,1);

for k = 1:100
    n = 10*k;

    [f_er,res] = func_A_no_piv(n);
    f_er_A_no_piv(k,1) = f_er;
    res_A_no_piv(k,1) = res;

    [f_er2,res2] = func_B_piv(n);
    f_er_B_piv(k,1) = f_er2;
    res_B_piv(k,1) = res2;

    [f_er3,res3] = func_B_no_piv(n);
    f_er_B_no_piv(k,1) = f_er3;
    res_B_no_piv(k,1) = res3;
end
figure()
plot([1:100],log10(f_er_A_no_piv),'--o');
hold on;
plot([1:100],log10(res_A_no_piv),'-*');
title('Ax=b 不选主元 向前误差（圆圈）与残差（星号），取对数');
hold off;

figure()
plot([1:100],log10(f_er_B_no_piv),'--o');
hold on;
plot([1:100],log10(res_B_no_piv),'-*');
title('Bx=b 不选主元 向前误差（圆圈）与残差（星号），取对数');
hold off;

figure()
plot([1:100],log10(f_er_B_piv),'--o');
hold on;
plot([1:100],log10(res_B_piv),'-*');
title('Bx=b 列主元 向前误差（圆圈）与残差（星号），取对数');
hold off;







%% 计算Ax=b，列主元高斯消元（理论上无需选主元）
function [f_error, residual] = func_A_piv(n)
% 首先创建A与b
A = diag(8*ones(n,1)) + diag(ones(n-1,1),1) + diag(6*ones(n-1,1),-1);
b = 15*ones(n,1);
b(1,1) = 9; 
b(n,1) = 14;
% 以下Gauss消元
count = 0;
for i = 1:n-1
    [~,index] = max(abs(A(i:i+1,i))); %三对角，保结构，只考虑i和i+1行
    index = i+index-1;
    if index ~= i %说明需要发生行交换
        count = count +1; %可以用于查看是否发生过行交换（理论上不需要）
        v = A(i,:);
        A(i,:) = A(index,:);
        A(index,:) = v; % 交换行，注意是整行交换，而不是在Schur补内交换
        temp = b(i,1); %对右侧向量b也要记得行交换
        b(i,1) = b(i+1,1);
        b(i+1,1) = temp;
    end
    A(i+1,i) = A(i+1,i)/A(i,i); %消元只局限于两行之间
    A(i+1,i+1:n) = A(i+1,i+1:n) - A(i+1,i)*A(i,i+1:n); % 修正Schur补
    b(i+1,1) = b(i+1,1) - b(i,1)*A(i+1,i); %A(i+1,i)已经被修正过，乃消元时系数
end
% 以下进行上三角线性系统Ux=b的求解
for i = n:-1:2
    b(i,1) = b(i,1)/A(i,i);
    b(i-1,1) = b(i-1,1) - b(i,1)*A(i-1,i);
end
b(1,1) = b(1,1)/A(1,1);
% 此时b储存的是Ax=b的解
standard_solution = ones(n,1); %真实的解
f_error = standard_solution - b; %向前误差
residual = norm(A*(f_error)); %残差的2范数
f_error = norm(f_error); %向前误差的2范数
end

%% 计算Ax=b，无选主元的高斯消元（理论上无需选主元）
function [f_error, residual] = func_A_no_piv(n)
% 首先创建A与b
A = diag(8*ones(n,1)) + diag(ones(n-1,1),1) + diag(6*ones(n-1,1),-1);
b = 15*ones(n,1);
b(1,1) = 9; 
b(n,1) = 14;
% 以下Gauss消元
for i = 1:n-1
    A(i+1,i) = A(i+1,i)/A(i,i); %消元只局限于两行之间
    A(i+1,i+1:n) = A(i+1,i+1:n) - A(i+1,i)*A(i,i+1:n); % 修正Schur补
    b(i+1,1) = b(i+1,1) - b(i,1)*A(i+1,i); %A(i+1,i)已经被修正过，乃消元时系数
end
% 以下进行上三角线性系统Ux=b的求解
for i = n:-1:2
    b(i,1) = b(i,1)/A(i,i);
    b(i-1,1) = b(i-1,1) - b(i,1)*A(i-1,i);
end
b(1,1) = b(1,1)/A(1,1);
% 此时b储存的是Ax=b的解
standard_solution = ones(n,1); %真实的解
f_error = standard_solution - b; %向前误差
residual = norm(A*(f_error)); %残差的2范数
f_error = norm(f_error); %向前误差的2范数
end

%% 计算By=b，列主元高斯消元
function [f_error, residual] = func_B_piv(n)
% 首先创建B与b
A = diag(6*ones(n,1)) + diag(ones(n-1,1),1) + diag(8*ones(n-1,1),-1);
b = 15*ones(n,1);
b(1,1) = 7; 
b(n,1) = 14;
% 以下Gauss消元
count = 0;
for i = 1:n-1
    [~,index] = max(abs(A(i:i+1,i))); %三对角，保结构，只考虑i和i+1行
    index = i+index-1;
    if index ~= i %说明需要发生行交换
        count = count +1; %可以用于查看是否发生过行交换（理论上不需要）
        v = A(i,:);
        A(i,:) = A(index,:);
        A(index,:) = v; % 交换行，注意是整行交换，而不是在Schur补内交换
        temp = b(i,1); %对右侧向量b也要记得行交换
        b(i,1) = b(i+1,1);
        b(i+1,1) = temp;
    end
    A(i+1,i) = A(i+1,i)/A(i,i); %消元只局限于两行之间
    A(i+1,i+1:n) = A(i+1,i+1:n) - A(i+1,i)*A(i,i+1:n); % 修正Schur补
    b(i+1,1) = b(i+1,1) - b(i,1)*A(i+1,i); %A(i+1,i)已经被修正过，乃消元时系数
end
% 以下进行上三角线性系统Ux=b的求解
for i = n:-1:2
    b(i,1) = b(i,1)/A(i,i);
    b(i-1,1) = b(i-1,1) - b(i,1)*A(i-1,i);
end
b(1,1) = b(1,1)/A(1,1);
% 此时b储存的是Bx=b的解
standard_solution = ones(n,1); %真实的解
f_error = standard_solution - b; %向前误差
residual = norm(A*(f_error)); %残差的2范数
f_error = norm(f_error); %向前误差的2范数
% count
end
%% 计算By=b，不选主元的高斯消元
function [f_error, residual] = func_B_no_piv(n)
% 首先创建B与b
A = diag(6*ones(n,1)) + diag(ones(n-1,1),1) + diag(8*ones(n-1,1),-1);
b = 15*ones(n,1);
b(1,1) = 7; 
b(n,1) = 14;
% 以下Gauss消元
for i = 1:n-1
    A(i+1,i) = A(i+1,i)/A(i,i); %消元只局限于两行之间
    A(i+1,i+1:n) = A(i+1,i+1:n) - A(i+1,i)*A(i,i+1:n); % 修正Schur补
    b(i+1,1) = b(i+1,1) - b(i,1)*A(i+1,i); %A(i+1,i)已经被修正过，乃消元时系数
end
% 以下进行上三角线性系统Ux=b的求解
for i = n:-1:2
    b(i,1) = b(i,1)/A(i,i);
    b(i-1,1) = b(i-1,1) - b(i,1)*A(i-1,i);
end
b(1,1) = b(1,1)/A(1,1);
% 此时b储存的是Bx=b的解
standard_solution = ones(n,1); %真实的解
f_error = standard_solution - b; %向前误差
residual = norm(A*(f_error)); %残差的2范数
f_error = norm(f_error); %向前误差的2范数
end