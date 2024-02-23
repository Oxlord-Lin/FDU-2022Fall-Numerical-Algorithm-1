%% 探索Arnodi过程中产生的向量的正交性
n = 250;
departure = zeros(n,1);
flg = 0; %用于标记是否因H(i+1,i)足够小而提前终止
A = randn(n);
r = rand(n,1);
Q = zeros(n,n); %用于存储Krylov子空间中的正交基向量
H = zeros(n,n);
Q(1:n,1) = r./norm(r);
departure(1,1) = abs(Q(:,1)'*Q(:,1) - 1);
for i = 1:n-1
    y = A*Q(1:n,i); % y=A*q_i
    for j = 1:i
        H(j,i) = Q(1:n,j)'*y; % 向前i个正交基向量上投影
        y = y - H(j,i)*Q(1:n,j);
    end
    H(i+1,i) = norm(y);
    norm(y);
    if abs(H(i+1,i)) < 1e-15 %若H(i,i+1)足够小
        flg = 1;
        disp("提前终止，因为H(i,i+1)足够小，认为r已经完全落在i维Krylov子空间中")
        break
    end
    Q(1:n,i+1) = y./H(i+1,i); %生成第i+1个正交基向量
    departure(i+1,1) = norm(Q(:,1:i+1)'*Q(:,1:i+1)-eye(i+1));
end
plot(log10(departure));
xlabel('迭代次数');
ylabel("已经生成向量的正交程度，取对数");
title("探索Arnoldi过程产生的向量的正交性");