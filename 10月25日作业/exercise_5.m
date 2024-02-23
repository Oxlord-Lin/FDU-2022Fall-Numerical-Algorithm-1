% %% 初始化（一）
m=35;
n=25;
U =orth(rand(m,m));%生成一组标准正交基u
V =orth(rand(n,n));%生成另一组标准正交基v
% 初始化（二）
step=0.1;
last=60;
cond=zeros(1,((last-1)/step+1)); %cond是行向量
dep=[cond;cond;cond;cond;cond;cond]; %一共六行，记录CGS,CGS_reorth,MGS,MGS_reorth,Cho,Hou的正交程度
residual = dep; %记录残差 （Matlab使用深拷贝）
count=1;
%% 计算部分
for k= 1:step:last
    Max_ex = 10^k;
    %最大的奇异值从10增加到10^last，按照指数的方式增长,指数上步长为step
    t = linspace(0,k,n);
    for j = 1:n
        t(1,j) = 10^t(1,j); %奇异值服从几何分布，属于病态矩阵
    end
    sig = diag(t);
%     sig = diag(1 + (Max_ex-1)*rand(n,1));
    sig(1,1) = 1;
    sig(n,n) = 10^k;
    cond(1,count)=sig(n,n);
    A = U(:,1:n)*sig*V';
    [dep(1,count),residual(1,count)] = CGS(A);
    [dep(2,count),residual(2,count)] = CGS_reorth(A);
    [dep(3,count),residual(3,count)] = MGS(A);
    [dep(4,count),residual(4,count)] = MGS_reorth(A);
    [dep(5,count),residual(5,count)] = QR_by_Cho(A);
    [dep(6,count),residual(6,count)] = QR_by_Hou(A);
    count=count+1;
end
%% 绘图（y轴取对数）
% 首先绘制departure图
x = cond;
y1 = log10(dep(5,:));
y2 = log10(dep(1,:));
y3 = log10(dep(3,:));
y4 = log10(dep(2,:));
y5 = log10(dep(4,:));
y6 = log10(dep(6,:));
plot(x,y1,'r');
hold on;
plot(x,y2,'y');
% legend('Cholesky','CGS')
plot(x,y3,'g');
% legend('MGS')
plot(x,y4,'c');
% legend('CGS_reorth')
plot(x,y5,'b');
% legend('MGS_reorth')
plot(x,y6,'k');
% legend('Householder')
legend('Cholesky','CGS','MGS','CGS-reorth','MGS-reorth','Householder');
% figure
% hold on
% plot(cond,log10(dep(5,:)),'Color','r'); %Cho
% plot(cond,log10(dep(1,:)),'Color','y'); %CGS
% legend('Cho','CGS');
% plot(cond,log10(dep(3,:)),'Color','g'); %MGS
% plot(cond,log10(dep(2,:)),'Color','c'); %CGS_reorth
% plot(cond,log10(dep(4,:)),'Color','b'); %MGS_reorth
% plot(cond,log10(dep(6,:)),'Color','k'); %Hou
% % legend('Cholesky','CGS','MGS','CGS-reorth','MGS-reorth','Householder')
title('departure');
xlabel('条件数k（A）');
ylabel('log10(departure), departure=norm(Q*Q-eye(n),fro)')

%% 然后绘制residual图
figure
hold on
y1 = log10(residual(5,:));
y2 = log10(residual(1,:));
y3 = log10(residual(3,:));
y4 = log10(residual(2,:));
y5 = log10(residual(4,:));
y6 = log10(residual(6,:));
plot(x,y1,'r');
hold on;
plot(x,y2,'y');
plot(x,y3,'g');
plot(x,y4,'c');
plot(x,y5,'b');
plot(x,y6,'k');
title('residual');
xlabel('条件数k（A）');
ylabel('log10(residual), departure=norm(Q*Q-eye(n),fro)');
legend('Cholesky','CGS','MGS','CGS-reorth','MGS-reorth','Householder');
%% 绘图（都取对数）
temp = [5,1,3,2,4,6]';
display_name = ['Cho';'CGS';'MGS';'CRO';'MRO';'Hou'];
figure
for k = 1:6
    index = temp(k,1);
    name = display_name(k,:);
    plot(log10(cond),log10((dep(index,:))));
    legend(name);
    hold on
end
title('departure');
xlabel('log10(条件数k(A))');
ylabel('log10(departure), departure=norm(Q*Q-eye(n),fro)');
legend('Cholesky','CGS','MGS','CGS-reorth','MGS-reorth','Householder');
%%
figure
for k = 1:6
    index = temp(k,1);
    plot(log10(cond),log10((residual(index,:))));
    hold on
end
title('residual');
xlabel('log10(条件数k(A))');
ylabel('log10(residual), departure=norm(Q*Q-eye(n),fro)');
legend('Cholesky','CGS','MGS','CGS-reorth','MGS-reorth','Householder');
%% 用经典GS，不进行重正交化
function [dep,residual] = CGS(A)
[m,n] = size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j = 1:n
    y=A(1:m,j);
    for i = 1:j-1
        R(i,j)=Q(1:m,i)'*A(1:m,j);
        y=y-Q(1:m,i).*R(i,j);%经典的G-S方法
    end
    R(j,j)=norm(y);
    Q(1:m,j)=y/R(j,j);
end
%正交化完成，计算残差、depature等
[dep,residual] = dep_and_res(A,Q,R);
end

% scatter(cond2,dep);
% title("经典GS,不用重正交化,departure of orthogonality")
% xlabel("K2(A)")
% f2=figure(2);
% scatter(cond2,residual);
% title("经典GS,不用重正交化,residual");
% xlabel("K2(A)")


%% 用经典GS，进行重正交化
function [dep,residual] = CGS_reorth(A)
[m,n] = size(A);
Q=zeros(m,n);
R=zeros(n,n);
R2=R;
for j = 1:n
    y=A(1:m,j);
    for i = 1:j-1
        R(i,j)=Q(1:m,i)'*A(1:m,j);%经典的G-S方法
        y=y-Q(1:m,i).*R(i,j);
        R2(i,j)=Q(1:m,i)'*y;%重正交时得对y进行投影
        y=y-Q(1:m,i).*R2(i,j);
    end
    R(j,j)=norm(y);
    Q(1:m,j)=y/R(j,j);
end
%正交化完成，计算残差、depature等
[dep,residual] = dep_and_res(A,Q,R);
end

%% 改进GS，不进行重正交化
function [dep,residual] = MGS(A)
[m,n] = size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j = 1:n
    y=A(1:m,j);
    for i = 1:j-1
        R(i,j)=Q(1:m,i)'*y;%改进的G-S方法
        y=y-Q(1:m,i).*R(i,j);
    end
    R(j,j)=norm(y);
    Q(1:m,j)=y/R(j,j);
end
%正交化完成，计算残差、depature等
[dep,residual] = dep_and_res(A,Q,R);
% figure(5);
% scatter(cond2,dep);
% title("改进GS,不进行重正交化,departure of orthogonality")
% xlabel("K2(A)")
% figure(6);
% scatter(cond2,residual);
% title("改进GS,不进行重正交化,residual");
% xlabel("K2(A)")
end
%% 改进GS，进行重正交化
function [dep,residual] = MGS_reorth(A)
[m,n] = size(A);
Q=zeros(m,n);
R=zeros(n,n);
R2=R;
for j = 1:n
    y=A(1:m,j);
    for i = 1:j-1
        R(i,j)=Q(1:m,i)'*y;%改进的G-S方法
        y=y-Q(1:m,i).*R(i,j);
        R2(i,j)=Q(1:m,i)'*y;%重正交化
        y=y-Q(1:m,i).*R2(i,j);
    end
    R(j,j)=norm(y);
    Q(1:m,j)=y/R(j,j);
end
%正交化完成，计算残差、depature等
[dep,residual] = dep_and_res(A,Q,R);
% figure(7);
% scatter(cond2,dep);
% title("改进GS,进行重正交化,departure of orthogonality")
% xlabel("K2(A)")
% figure(8);
% scatter(cond2,residual);
% title("改进GS,进行重正交化,residual");
% xlabel("K2(A)");
end
%% 使用cholesky分解，不进行重正交化(我也不知道怎么reorthogonalization)
function [dep,residual] = QR_by_Cho(A)
[~,n] = size(A);
S = A'*A;
for i = 1:n
    for j = i+1:n
        S(j,j:n) = S(j,j:n) - S(i,j:n)*S(i,j)'/S(i,i); % 此处勿忘转置S（i，j）
    end
    S(i,i:n) = S(i,i:n)./sqrt(S(i,i));
end
R = triu(S);
Q = A/R;
%正交化完成，计算残差、depature等
[dep,residual] = dep_and_res(A,Q,R);
end
% figure(9);
% scatter(cond2,dep);
% title("Cholesky分解,不进行重正交化,departure of orthogonality")
% xlabel("K2(A)")
% figure(10);
% scatter(cond2,residual);
% title("Cholesky分解，不进行重正交化,residual");
% xlabel("K2(A)");

%% 使用cholesky分解，进行重正交化
% Q=zeros(m,n);
% R=zeros(n,n);
% R2=R;
% count=1;
% for k= 1:step:last
%     sig(1,1)=10^k;%最大的奇异值从10增加到10^last，按照指数的方式增长,指数上步长为step
%     A = U*sig*V';%生成A
%     R=chol(A'*A);
%     %此处存疑，cholesky分解的重正交化怎么做？
%     R=chol(R'*R);
%     Q=A/R;
%     
% 
% %正交化完成，计算残差、depature等
% dep(1,count)=norm((Q'*Q-eye(n)),'fro');
% residual(1,count)=norm((A-Q*R),'fro');
% count=count+1;
% end
% figure(11);
% scatter(cond2,dep);
% title("Cholesky分解,进行重正交化,departure of orthogonality")
% xlabel("K2(A)")
% figure(12);
% scatter(cond2,residual);
% title("Cholesky分解，进行重正交化,residual");
% xlabel("K2(A)");

%% 使用householder方法进行QR分解
function [dep,residual] = QR_by_Hou(A)
[m,n] = size(A);
A0 = A;
Qt = eye(m);
for j =1:n
    v = A(j:m,j);
    len = norm(v);
    if v(1,1)>=0%此处相当于对v与一个长度为||v||且只有第一个分量不为零的向量作差
        v(1,1)=v(1,1)+len;%避免舍入误差
    else
        v(1,1)=v(1,1)-len;%避免舍入误差
    end
    v = v/norm(v);
    A(j:m,j:n) = A(j:m,j:n) - 2*v*(v'*A(j:m,j:n));
    Qt(j:m,:) = Qt(j:m,:) - 2*v*(v'*Qt(j:m,:));
end
Q = Qt';
Q = Q(:,1:n);
R = A(1:n,1:n);
%正交化完成，计算残差、depature等
[dep,residual] = dep_and_res(A0,Q,R);
end

%% 用于计算正交程度和残差的辅助函数
function [dep,residual] = dep_and_res(A,Q,R)
[~,n] = size(A);
dep = norm((Q'*Q-eye(n)),'fro');
residual = norm((A-Q*R),'fro');
end
