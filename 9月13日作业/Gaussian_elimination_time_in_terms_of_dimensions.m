% 本程序用于对非奇异阵进行高斯消元来得到上三角阵U
Max=600;%矩阵的维度上限
Min=250;%矩阵维度的下限。这是因为，如果维度太少，消元的运行时间非常短，接近0，会出现Inf。
t=zeros(Max-Min+1,1);

for dim=Min:Max %矩阵的维度变化范围
    A=randn(dim,dim);%随机生成(dim*dim)方阵
    %以下部分对测试矩阵A进行【高斯消元】变成上三角矩阵并记录所需时间
    t0=clock;%本轮高斯消元开始的时间
    for j=1:dim
        A(j+1:dim,j)=A(j+1:dim,j)./A(j,j);%把第j列变成各行之间的"系数之比"
        A(j+1:dim,j+1:dim)=A(j+1:dim,j+1:dim)-A(j+1:dim,j)*A(j,j+1:dim);%对Schur complement进行消元操作
        duration=etime(clock,t0);%本轮高斯消元结束的时间
        t(dim-Min+1,1)=duration;
    end
end

scatter(log(Min:Max),log(t))
xlabel("log(dimension)")
ylabel("log(execution time)")