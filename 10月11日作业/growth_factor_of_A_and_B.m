rho_A = zeros(100,1);
rho_B_piv = zeros(100,1);
rho_B_no_piv = zeros(100,1);
for k = 1:100
    n = 10*k;
   rho_A(k,1) = func_A_piv(n);
   rho_B_piv(k,1) = func_B_piv(n);
   rho_B_no_piv(k,1) = func_B_no_piv(n);
end
plot([1:100],(rho_A));
title('A（不选主元）的增长因子');
figure
plot([1:100],(rho_B_piv));
title('B（列主元）的增长因子');
figure
plot([1:100],(rho_B_no_piv));
title('B（不选主元）的增长因子');




%% 计算Ax=b，不选主元的高斯消元（理论上无需选主元），并返回增长因子
function rho = func_A_piv(n)
% 首先创建A
A = diag(8*ones(n,1)) + diag(ones(n-1,1),1) + diag(6*ones(n-1,1),-1);
A_max = max(max(A));
U_max = A_max;
% 以下Gauss消元
for i = 1:n-1
    A(i+1,i) = A(i+1,i)/A(i,i); %消元只局限于两行之间
    A(i+1,i+1:n) = A(i+1,i+1:n) - A(i+1,i)*A(i,i+1:n); % 修正Schur补
    A(i+1,i) = 0;
    temp = max(max(A));
    if U_max < temp
       U_max = temp;
    end
end
rho = U_max/A_max;
end


%% 计算By=b，列主元高斯消元，并返回增长因子
function rho = func_B_piv(n)
% 首先创建B
A = diag(6*ones(n,1)) + diag(ones(n-1,1),1) + diag(8*ones(n-1,1),-1);
A_max = max(max(A));
U_max = A_max;
% 以下Gauss消元
for i = 1:n-1
    [~,index] = max(abs(A(i:i+1,i))); %三对角，保结构，只考虑i和i+1行
    index = i+index-1;
    if index ~= i %说明需要发生行交换
        v = A(i,:);
        A(i,:) = A(index,:);
        A(index,:) = v; % 交换行，注意是整行交换，而不是在Schur补内交换
    end
    A(i+1,i) = A(i+1,i)/A(i,i); %消元只局限于两行之间
    A(i+1,i+1:n) = A(i+1,i+1:n) - A(i+1,i)*A(i,i+1:n); % 修正Schur补
    A(i+1,i) = 0;
    temp = max(max(A));
    if U_max < temp
       U_max = temp;
    end
end
rho = U_max/A_max;
end

%% 计算By=b，不选主元高斯消元，并返回增长因子
function rho = func_B_no_piv(n)
% 首先创建B
A = diag(6*ones(n,1)) + diag(ones(n-1,1),1) + diag(8*ones(n-1,1),-1);
A_max = max(max(A));
U_max = A_max;
% 以下Gauss消元
for i = 1:n-1
    A(i+1,i) = A(i+1,i)/A(i,i); %消元只局限于两行之间
    A(i+1,i+1:n) = A(i+1,i+1:n) - A(i+1,i)*A(i,i+1:n); % 修正Schur补
    A(i+1,i) = 0;
    temp = max(max(A));
    if U_max < temp
       U_max = temp;
    end
end
rho = U_max/A_max;
end