%% 探索不同阶数的Arnodi procedure产生的向量的正交性
begin = 5;
ending = 250;
depature = zeros((ending - begin + 1),1);
for n = begin:ending
    flg = 0; %用于标记是否因H(i+1,i)足够小而提前终止
    A = rand(n);
    r = rand(n,1);
    Q = zeros(n,n); %用于存储Krylov子空间中的正交基向量
    H = zeros(n,n);
    Q(1:n,1) = r./norm(r);
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
    end
    if flg ==1 %如果是提前终止的
        depature(n-begin+1,1) = norm(Q(1:n,1:i)'*Q(1:n,1:i)-eye(i));
    else
        depature(n-begin+1,1) = norm(Q'*Q-eye(n));
    end
end
figure();
plot([begin:ending],log10(depature));
xlabel("矩阵A的阶数");
ylabel("正交程度，通过log10(norm(Q'*Q-eye(n))进行衡量");
title("探索Arnoldi过程产生的正交向量的正交性与阶数n之间的关系");