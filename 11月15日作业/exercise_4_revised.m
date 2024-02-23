% 提示：所需的函数作为局部函数定义在后面
% 下面是正式的试验部分
%% 首先对5阶，10阶，15阶的随机生成的对称矩阵使用Jacobi对角化，并展示其运行结果
for i = 1:5
    n = 5*i; % n是阶数
    % 首先生成n阶对称实矩阵
    lambda=diag([1:n]);
    temp=rand(n,n);
    temp=orth(temp);
    A=temp*lambda*temp';
    % 此处临时加一个随机生成对称矩阵的部分
    for j=1:n
        for k = j:n
            A(j,k)=10*rand()*cos(rand());
            A(k,j)=A(j,k);
        end
    end
%     disp('阶数:');
%     n;
%     disp('初始时A为：')
%     A; %展示A
    % 接下来进行Jacobi对角化过程
    count = 0;
    while(1)
        count = count + 1;
        [p,q,c,s] = find_max_offdiag_and_theta(A);
        Q = Jacobi_rotation(A,p,q,c,s);
        A_2 = Q'*A*Q;
        if F_norm_of_offdiag_entries(A_2) <= (n*(n-1))*(1e-16)
            A = A_2;
            break
        end
        A = A_2;
    end
%     disp('Jacobi对角化的结果为：');
    figure
    imagesc(A)
    title(['阶数:',num2str(n),' 迭代次数：',num2str(count)]);
end
%% 下面对阶数n=20的收敛过程进行可视化展示（用off-diagnal的平方和表示）
disp('以下通过观察非对角元F范数的收敛过程的方式来展现一个阶数n=20实对称矩阵的Jacobi对角化过程')
n = 20;
lambda=10*diag(rand(n,1));
temp=rand(n,n);
temp=orth(temp);
A=temp*lambda*temp;
% 此处临时加一个随机生成对称矩阵的部分
for j=1:n
    for k = j:n
        A(j,k)=10*rand()*cos(rand());
        A(k,j)=A(j,k);
    end
end
% disp('初始时的A如下')
% A
er = [F_norm_of_offdiag_entries(A)];
count = 1;
while(1)
    count = count + 1;
    [p,q,c,s] = find_max_offdiag_and_theta(A);
    Q = Jacobi_rotation(A,p,q,c,s);
    A_2 = Q'*A*Q;
    E = F_norm_of_offdiag_entries(A_2);
    er(1,count) = E;
    if E <= (n*(n-1)/2)*(1e-15)
        A = A_2;
        break
    end
    A = A_2;
end
% disp('Jacobi对角化的结果为：');
% A
disp('非对角元的F范数收敛过程如下')
figure();
plot(er);
title(['非对角元F范数收敛过程',' 阶数',num2str(n),' 迭代次数',num2str(count-1)]);
xlabel('迭代次数')
ylabel('非对角元F范数')

%% 下面对阶数n=40的收敛过程进行可视化展示（entry_wise）
disp('以下通过可视化每个元素的变化的方式来可视化Jacobi对角化过程')
p = 7;
q = 8;
entry_pq = zeros(1,1); % 用于跟踪某个元素的变化情况
n = 40;
lambda=10*diag(randn(n,1));
temp=randn(n,n);
temp=orth(temp);
A=temp*lambda*temp';
% 此处临时加一个随机生成对称矩阵的部分
% for j=1:n
%     for k = j:n
%         A(j,k)=10*rand()*cos(rand());
%         A(k,j)=A(j,k);
%     end
% end
% disp('初始时的A如下')
% A
figure()
imagesc(abs(A))
colorbar
title("初始时A的元素分布如下")
count = 0;
entry_pq(count+1,1) = abs(A(p,q));
while(1)
    count = count + 1;
    [p,q,c,s] = find_max_offdiag_and_theta(A);
    Q = Jacobi_rotation(A,p,q,c,s);
    A_2 = Q'*A*Q;
    E = F_norm_of_offdiag_entries(A_2);
    entry_pq(count+1,1) = abs(A_2(p,q));
    if E <= (n*(n-1)/2)*(1e-16)
        A = A_2;
        break
    end
    if count==1 || count==3 || count==5 || count==7 || count==10 || count==15 || count==20
        figure()
        imagesc(abs(A_2))
        colorbar
        [~,~] = title('迭代次数为：',count);

    end
    A = A_2;
end
% disp('Jacobi对角化的结果为：');
% A
figure()
imagesc(abs(A))
colorbar
[~,~] = title('Jacobi对角化最终结果的分布情况',['迭代次数：',num2str(count)]);
%% 
figure
plot(log10(entry_pq(1:100,1)));
title(['(',num2str(p),',',num2str(q),')','位置上元素绝对值（取对数）的变化'])

%% 正所谓，工欲善其事，必先利其器，首先准备好需要的函数，如下
function Jacobi_Q = Jacobi_rotation(A,p,q,c,s) % 返回一个Jacobi旋转阵
n = size(A);
n = n(1,1);
ep = zeros(n,1);
ep(p,1) = 1;
eq = zeros(n,1);
eq(q,1) = 1;
Jacobi_Q = eye(n)+(c-1)*(ep*ep'+eq*eq')+s*(ep*eq'-eq*ep');
end

function [p,q,c,s] = find_max_offdiag_and_theta(A) % 找到最大的非对角元并返回cos(theta)&sin(theta)
n = size(A);
n = n(1,1);
p=1;
q=2;
max = abs(A(1,2));
for i =  1:n-1
    for j = i+1:n
    if abs(A(i,j))>abs(max)
        max = A(i,j);
        p = i;
        q = j;
    end
    end
end
tao = (A(q,q)-A(p,p)) / (2*A(p,q));
t = sign(tao)/(abs(tao) + sqrt(1+tao*tao));
c = 1 / sqrt(1+t*t);
s = t*c;
end

function E = F_norm_of_offdiag_entries(A) % 返回所有非对角元的平方和
E_2 = norm(A,"fro")^2 - norm(diag(A),"fro")^2;
E = sqrt(E_2);
end