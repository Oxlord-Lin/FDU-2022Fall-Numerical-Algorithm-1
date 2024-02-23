%% 循环Jacobi对角化方法，并展示非对角元范数收敛过程，以及按元素收敛过程
n=50;
lambda = diag(randn(n,1));
Q = orth(randn(n));
A_ori = Q*lambda*Q';
% disp('初始时的A如下')
% imagesc(A_ori);

%% 正常的扫描次序
A = A_ori;
figure()
imagesc(A);
colorbar
title('初始状态');
count=0; % 用于统计扫描次数
er=[F_norm_of_offdiag_entries(A)];
while (1)
    count=count+1; % 扫描次数加一
    % 以下是"一轮扫描"的过程
    for j = 1:n-1 %外面的大循环是列循环
        for i = j+1:n
            Q = Jacobi_rotation(A,i,j);
            A = Q'*A*Q; %进行一次Jacobi旋转
        end
    end
%     figure()
%     imagesc(A);
%     colorbar
%     title(['扫描次数为:',num2str(count)]);
    E=F_norm_of_offdiag_entries(A) % 计算非对角元的F范数
    er(1,count+1) = E;
    if E<(1e-16)
        break
    end
end
figure()
plot([0:count],log10(er),'o--')
xlabel('扫描次数')
ylabel('log(非对角元F范数)')

%% 【不正常的】的扫描次序（一）(外循环沿着行走，内循环沿着列走)(也能成功)
A = A_ori;
count=0; % 用于统计扫描次数
er=[F_norm_of_offdiag_entries(A)];
while (1)
    count=count+1; % 扫描次数加一
    % 以下是"一轮扫描"的过程
    for i = 2:n %外面的大循环是行循环
        for j = 1:(i-1)
            Q = Jacobi_rotation(A,i,j);
            A = Q'*A*Q; %进行一次Jacobi旋转
        end
    end
    E=F_norm_of_offdiag_entries(A); % 计算非对角元的F范数
    er(1,count+1) = E;
    if E<(1e-16)
        break
    end
end
figure()
plot([0:count],log10(er),'o--')
xlabel('扫描次数')
ylabel('log(非对角元F范数)')


%% 【不正常的】的扫描次序（二）(随机选取p,q)(试了几个例子，能够成功，相当于扫描16轮)
A = A_ori;
% figure()
% imagesc(A);
% colorbar
% title('初始状态');
% pause(1)
count=0; % 用于统计旋转次数
count2 = 1;
er=[F_norm_of_offdiag_entries(A)];
flg = 0;
while (count<1e5)
    count=count+1; % 扫描次数加一
    p = randi(n);
    q = randi(n);
    Q = Jacobi_rotation(A,p,q);
    A = Q'*A*Q; %进行一次Jacobi旋转
    if mod(count,10) == 0
        count2 = count2 + 1;
        E = F_norm_of_offdiag_entries(A); % 计算非对角元的F范数
        er(1,count2+1) = E;
        if E<1e-16
            flg = 1;
        end
    end
    if flg == 1
        disp('达到收敛条件')
        disp(count)
        break 
    end

end
figure()
plot(log10(er))
ylabel('log10(非对角元F范数)')

%% 【不正常的】的扫描次序（二）（在一轮扫描中顺序随机，但不重不漏）(试验了几个例子，也能收敛)
A = A_ori;
figure()
imagesc(A);
colorbar
title('初始状态');
E=F_norm_of_offdiag_entries(A)
pause(1)
count=0; % 用于统计扫描次数
er=[F_norm_of_offdiag_entries(A)];
while (1)
    count=count+1; % 扫描次数加一
    % 以下是"一轮扫描"的过程
    row = randperm(n);
    col = randperm(n);
    for i = 1:n 
        for j = 1:n
            p = row(1,i);
            q = col(1,j);
            if i>=j % 只需要考虑下三角部分即可
                continue
            end
            Q = Jacobi_rotation(A,p,q);
            A = Q'*A*Q; %进行一次Jacobi旋转
        end
    end
    E = F_norm_of_offdiag_entries(A) % 计算非对角元的F范数
    er(1,count+1) = E;

    if E<(1e-16)
        break
    end
end
figure()
plot([0:count],log10(er),'o--')
xlabel('扫描次数')
ylabel('log10(非对角元F范数)')

%% 局部函数如下 
function Jacobi_Q = Jacobi_rotation(A,p,q) % 返回一个对应于（p,q）的Jacobi旋转阵
n = size(A);
n = n(1,1);
ep = zeros(n,1);
ep(p,1) = 1;
eq = zeros(n,1);
eq(q,1) = 1;
tao = (A(q,q)-A(p,p)) / (2*A(p,q));
t = sign(tao)/(abs(tao) + sqrt(1+tao*tao));
c = 1 / sqrt(1+t*t);
s = t*c;
Jacobi_Q = eye(n)+(c-1)*(ep*ep'+eq*eq')+s*(ep*eq'-eq*ep');
end

function E = F_norm_of_offdiag_entries(A) % 返回所有非对角元的F范数
E_2 = norm(A,"fro")^2 - norm(diag(A),"fro")^2;
E = sqrt(E_2);
end