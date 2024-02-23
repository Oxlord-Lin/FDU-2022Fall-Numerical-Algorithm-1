%% 预条件反幂法
n = 2e3;
D = diag(randn(1,n));
U = orth(randn(n));
A = U*D*U';
X1 = zeros(n,1); %存储带预条件的反幂法产生的向量
X2 = zeros(n,1); %存储常规反幂法产生的向量
X3 = zeros(n,1); %存储利用Rayleigh-Rits进行改造后的反幂法产生的向量
%% 带预条件
theta = 1;
v = randn(n,1);
v = v./norm(v);
X1(1:n,1) = v;
for count = 2:1000 %迭代次数保护
    rayleigh = (v'*A*v)/(v'*v);
    r = A*v-v*rayleigh;
    v2 = v - (A-theta*eye(n))\r;
    v2 = v2./norm(v2);
    X1(1:n,count) = v2;
    v = v2;
    if norm(X1(1:n,count-1)-X1(1:n,count)) < n*1e-16
        break
    end
end
for j = 1:100:n
    plot(X1(j,1:count))
    hold on
end
xlabel('iter');
title('带预条件反幂法产生的向量的部分分量收敛过程');

%% 常规反幂法
w = randn(n,1);
w = w./norm(w);
X2(1:n,1) = w;
for count = 2:1000 %迭代次数保护
    w2 = A\w;
    w2 = w2./norm(w2);
    X2(1:n,count) = w2;
    w = w2;
    if norm(X2(1:n,count-1)-X2(1:n,count)) < n*1e-16
        break
    end
end
figure()
for j = 1:100:n
    plot(X2(j,1:count))
    hold on
end
xlabel('iter');
title('常规反幂法产生的向量的部分分量收敛过程');

%% 利用Rayleigh-Rits过程来优化预条件反幂法
theta = 1;
x = randn(n,1);
x = x./norm(x);
X3(1:n,1) = x;
for count = 2:1000 %迭代次数保护
    rayleigh = (x'*A*x)/(x'*x);
    r = A*x-x*rayleigh;
    x2 = (A-theta*eye(n))\r;
    % 进行Gram-Schmidt正交化
    x2 = x2 - x'*x2;
    x2 = x2./norm(x2);
    % Rayleigh-Rits过程
    S = [x,x2]; % S包含两个正交向量
    T = S'*A*S;
    [U2,D2] = eig(T);
    [~,ind] = sort(diag(D2)); %默认升序排列
    Ds = D2(ind,ind);
    Us = U2(:,ind);
    Rits_vector = S*Us;
    X3(1:n,count) = Rits_vector(1:n,1);
    if norm(X3(1:n,count-1)-X3(1:n,count)) < n*1e-16
        break
    end
end
figure()
for j = 1:100:n
    plot(X3(j,1:count))
    hold on
end
xlabel('iter');
title('利用Rayleigh-Rits进行改造后的带预条件反幂法产生的向量的部分分量收敛过程');



