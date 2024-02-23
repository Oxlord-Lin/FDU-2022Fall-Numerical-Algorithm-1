%% 用于计算对称矩阵A的特征值与特征向量的Lanczos算法
% 第一种生成A的方式
n = 2e3;
D1 = diag(rand(1,n));
U1= orth(randn(n));
A = U1*D1*U1'; 
%for test
% A = D1;

% 另一种生成A的方式
% for i = 1:n
%     for j=1:i
%         A(i,j) = randn();
%         A(j,i) = A(i,j);
%     end
% end
% [U1,D1] = eig(A);
% [~,ind] = sort(diag(D1));
% D1 = D1(ind,ind);
% U1 = U1(:,ind);

figure()
plot([1:n],diag(D1));
title('eigenvalues of A');
xlabel('index');
Max_Rits_value = zeros(1,n);
min_Rits_value = zeros(1,n);
middle_Rits_value = zeros(1,n);
Max_Rits_vector = zeros(n,n);
min_Rits_vector = zeros(n,n);

orthogonality_of_Lanczos_vector = zeros(1,n); %用于记录Lanczos向量的正交性

Q = zeros(n,n);
T = zeros(n,n);
b = randn(n,1);
Q(1:n,1) = b./norm(b);

for i = 1:n-1
    y = A*Q(1:n,i); % y=A*q_i
    T(i,i) = Q(1:n,i)'*y;
    y = y - T(i,i)*Q(1:n,i);
    if i-1>=1
        y = y - T(i,i-1)*Q(1:n,i-1);
    end
    T(i+1,i) = norm(y);
    if  T(i+1,i)<1e-16
        disp('提前结束');
        break
    end
    T(i,i+1) = T(i+1,i);
    Q(1:n,i+1) = y./T(i+1,i); %生成第i+1个正交基向量

%   计算Tj的特征值（Rits value），特征向量(Rits vector)，并计算Lancsoz向量的正交性
    
    orthogonality_of_Lanczos_vector(1,i) = norm(Q(1:n,1:i)'*Q(1:n,1:i)-eye(i));
    [U2,D2]=eig(T(1:i,1:i));
%     for t = 1:i
%         if i>30
%             break
%         end
%         scatter(i,D2(t,t));
%         pause(0.1)
%         hold on
%     end
    D2 = diag(D2);
    D2 = sort(D2);
    a = max(D2);
    Max_Rits_value(1,i) = a;
    b = min(D2);
    min_Rits_value(1,i) = b;
    if i~=1
        c = D2(round(i/2));
        middle_Rits_value(1,i) = c;
    end
    Rits_vector = (Q(1:n,1:i)*U2);
    Max_Rits_vector(1:n,i) = Rits_vector(1:n,1);
    min_Rits_vector(1:n,i) = Rits_vector(1:n,i);
    
end

%由于i最大到n-1，故最后一次正交性得单独算
orthogonality_of_Lanczos_vector(1,n) = norm(Q(1:n,1:n)'*Q(1:n,1:n)-eye(n)); 

figure()
bar([1:n],Max_Rits_value);
figure()
bar([1:n],min_Rits_value);
figure()
bar([1:n],middle_Rits_value);

% lanczos向量的正交性
figure()
bar([1:n],orthogonality_of_Lanczos_vector,10);
title('orthogonality of Lanczos vectors');
xlabel('iter');
ylabel('Lanczos vectors的正交性');

% 排序后的Rits values(即T的特征值，作为A的特征值的近似)
figure();
[U2,D2] = eig(T);
[d,ind] = sort(diag(D2));
D2_sorted = D2(ind,ind);
U2_sorted = U2(:,ind);
plot(diag(D2_sorted));
title('Rits values(sorted)')
xlabel('index')

% Rits value 的收敛性
figure()
plot([1:n],diag(abs(D2-D1))); % Rits value的绝对误差
title('convergence of Rits value')
xlabel('index')
ylabel('absolute error of Rits value')

%Rits vector 的收敛性
figure()
Rits_vector = Q*U2_sorted;
convergence_of_Rits_vector = zeros(1,n);
for i = 1:n
    convergence_of_Rits_vector(1,i) = norm(Rits_vector(1:n,i)-U1(1:n,i));
end
plot([1:n], convergence_of_Rits_vector);
title('convergence of Rits vector');

%Rits vector收敛性 （个人感觉更直观，只考虑最大和最小Rits value对应的Rits vector）
%最大
count=1;
Max_Rits_vector_2 = zeros(n,1);
for i = 1:n
    if norm(Max_Rits_vector(1:n,i))>1e-16
        Max_Rits_vector_2(1:n,count) = Max_Rits_vector(1:n,i);
        count = count+1;
    end
end

convergence_of_Max_Rits_vector = Max_Rits_vector_2(1:5,1:end);

figure()
for i =1:5
    plot(abs(convergence_of_Max_Rits_vector(i,1:end)));
    hold on
end
title('最大Rits vector收敛过程');

%最小
count=1;
min_Rits_vector_2 = zeros(n,1);
for i = 1:n
    if norm(min_Rits_vector(1:n,i))>1e-16
        min_Rits_vector_2(1:n,count) = min_Rits_vector(1:n,i);
        count = count+1;
    end
end
convergence_of_min_Rits_vector = min_Rits_vector_2(end-4:end,1:end);

figure()
for i =1:5
    plot(abs(convergence_of_min_Rits_vector(i,1:end)));
    hold on
end
title('最小Rits vector收敛过程');