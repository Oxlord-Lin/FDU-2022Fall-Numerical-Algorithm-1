%% QR_Algorithm_for_symmatric_matrix

%% 准备工作：生成对称矩阵A
n = 20; %n为矩阵阶数
A = zeros(n);
for i = 1:n
    for j = 1:i
        A(j,i)=-2+4*rand();
        A(i,j) = A(j,i);
    end
end

% M用于记录T33阶数的变化（即已经算出的特征值的个数）
M=zeros(n,1);

%% 准备工作：三对角化（调用自己写的函数）
    % T与U满足：T=U'AU
[T,U] = upper_hessenberg(A); %对A进行上海森伯格化，由于A是对称的，故实际上为三对角化

%% QR算法的主体部分
tol = 1e-14;
Q = eye(n); %用于储存所有的Givens旋转阵
max = 1000; %迭代次数保护
m=0; % m是T33阶数，即已经算出的特征值的个数
for count = 1:max
    % 一、收敛性判定
    % 把足够小的次对角元变成0
    for i = 1:n-m-1
        abs(T(i+1,i))
        if abs(T(i+1,i))<=(abs(T(i,i))+abs(T(i+1,i+1)))*tol
            T(i+1,i) = 0;
            T(i,i+1) = 0;
        end
    end
    
    % 二、对T进行分块（调用自己写的函数）
        % 从左上角到右下角分别记作T11，T22，T33，其中T33是对角阵，T22是最大的不可约三对角阵
        % 分块的目的以及细节可参见北大教材《数值线性代数》对称QR算法有关章节
    [l,m] = block_size(T); %l是T11阶数，m是T33阶数
    M(count,1) = m; %记录已经得到的特征值的个数
    if m==n  % 此时已经是对角阵，得到全部特征值
        break;
    end

    % 三、对T22进行一次Wilkinson位移隐式QR迭代（调用自己写的函数）
    [T(l+1:n-m,l+1:n-m),G] = one_step_QR_with_Wilkinson_shift(T(l+1:n-m,l+1:n-m));
    
    % 四、储存旋转阵
    Q = Q*blkdiag(eye(l),G,eye(m)); %Q用于存储Givens旋转阵
    
end %最外层循环结束（有迭代次数保护的那个）

Q = U*Q; %把三对角化时的U也包含进来，此处的Q满足 A=QDQ',其中D=T为对角阵，元素为特征值，得到A的谱分解

%% 对称QR算法成果展示
plot([1:count],M,"-o"); %迭代次数和T33（已经算出的特征值）大小的关系
xlabel("rounds of iteration")
ylabel("the number of obtained eigenvalue")
title({['QR algorithm for symmatric matrix'];['matrix size:',num2str(n),' by ',num2str(n)];['tol = ',num2str(tol)]});

disp(["Average iter for an eigenvalue:",count/n])

disp(["Depature:",norm(Q*T*Q'-A)])