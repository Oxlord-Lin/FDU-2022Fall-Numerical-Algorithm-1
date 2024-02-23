%先随机产生一个实对称矩阵
Q = orth(rand(200,200));
D = diag([1:200]);
A = Q*D*Q';
%选择lambda=69作为试验对象，对应的特征向量就是Q的第69列
v=Q(1:200,69) + 0.01*(rand(200,1));%对特征向量加扰动
%v=20*rand(200,1);
flg=0;
count=0;
lambda_list=zeros(100,1);
for i = 1:(100)%迭代保护
    count=count+1;
    lambda=(v'*A*v)/(v'*v)%瑞利商
    lambda_list(count,1)=lambda;
    y=(A-lambda*eye(200))\v;
    r=y-lambda*v;%计算残差
    v=y/norm(y);%归一化，以备进入下一轮循环
    if count>1 && abs(lambda_list(count,1)-lambda_list(count-1,1))<(1e-16)
        disp("两次计算的瑞利商足够接近，结束瑞利商迭代")
        flg=1;
        break
    end
end
if flg==0
    disp("触发迭代保护")
end
if flg==1
    plot([1:count],lambda_list(1:count,1),'b--o');
    title("使用瑞利商计算特征值的收敛过程")
    xlabel("迭代次数")
    ylabel("瑞利商近似值")
end

