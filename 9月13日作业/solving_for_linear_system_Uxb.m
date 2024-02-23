%% 【准备部分】
% 本题中，n可以取任何值，在本程序中，不妨令n=10
n=10;
%以下部分用于随机生成一个非奇异的上三角阵U
while true;
    A=rand(n,n);
    [L,U]=lu(A);%将A分解为LU，实际上只需要使用U，也就是题目中的Ux=b中的U，此时U应该是非奇异的。
    if abs(det(U)-0)>1.0e-06;
        break %防止出现奇异阵
    end
end
% 以下部分用于取一个b
b0 = rand(n,1);%由于U是非奇异的，所以b可以任意取

%% 【方法一】从第n行开始回代
b = b0;
for j=n:-1:1
    b(j,1)=b(j,1)/U(j,j);%得出的xj储存在b中
    for i = 1:j-1
        b(i,1)=b(i,1)-U(i,j)*b(j,1);%逐行消元
    end
end
b1 = b; %最终的结果储存在向量b中，这一步用于打印出结果

%% 【方法二】直接用克莱默法则暴算
b = b0;
Dn=det(U);
for j=1:n
    Ubj=U;
    Ubj(1:n,j)=b(1:n,1);
    Dj=det(Ubj);
    b(j,1)=Dj/Dn;
end
b2 = b; %最终的结果储存在向量b中，这一步用于打印出结果
