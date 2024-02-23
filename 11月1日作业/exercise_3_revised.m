%% 先随机产生一个实对称矩阵
Q = orth(rand(200,200));
D = diag(linspace(2,400,200));
A = Q*D*Q';
%% Rayleigh Quotient Iteration
flg = 0;
NA = norm(A);
lambda = 100+0.1;
x = rand(200,1);
x = x/norm(x);%归一化
lambda_list=zeros(10^2,1);
count = 1;
lambda_list(count,1) = x'*A*x;
for i = 1:(10^2)%迭代保护
    count = count + 1;
    B = A - lambda*eye(200);
    y = B\x;
    lambda = y'*A*y/(y'*y); %Rayleight Quotient
    lambda_list(count,1) = lambda;
    r = A*x - lambda*x; %计算残差
    if norm(r) <= (NA + abs(lambda))*1e-15
        flg = 1;
        disp("残差r足够小，结束反幂法迭代")
        break
    end
    x = y/norm(y); % 归一化，准备进入下一轮循环
end

%% 绘图
plot([1:count],lambda_list(1:count,1),'b--o');
title(['使用Rayleigh商迭代计算特征值的收敛过程,选取lambda=100']);
xlabel("迭代次数")
ylabel("特征值的近似值")


