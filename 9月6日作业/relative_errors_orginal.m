k=2.^10;
x=[rand();rand()];%随机生成向量x
y=x';%转置
%对绝对误差和相对误差e与er（都是k+1维数组）进行初始化
e=zeros(k+1,1);
er= zeros(k+1,1);

for n =0:k;
    theta=2*pi*(n/k);
    A=[cos(theta),sin(theta);-sin(theta),cos(theta)];
    tr=(y*x)*cos(theta); %真实值（或者说尽可能接近真实的值），按照向量的投影后的几何关系计算
    fl=single(y)*single(A)*single(x); %机器按照矩阵乘法计算后所得浮点数
    e(n+1,1)=abs(fl-tr);%绝对误差
    er(n+1,1)=e(n+1,1)./(abs(tr)); %相对误差，此处为了监控计算过程，让er输出，没有打分号
end
plot(0:k,er)%绘制图表，观察er随着i的变化情况
xlabel("n")
ylabel("relative error")