format long;
A=rand(1000,1000);
x=rand(1000,1);
norm_of_A=norm(A);
x=x/norm(x);%输入的x进行归一化
flg=0;
count=0;
lambda_list=zeros(10^3,1);
for i=1:(10^3) %迭代次数保护
    count=count+1;
    y=A*x;
    L_square=x'*x;
    lambda=(x'*y)/L_square; %瑞利商
    lambda_list(count,1)=lambda;%记录一下本次获得的特征值的近似值
    r=y-lambda*x;%残差
    x=y/norm(y);%使用2-范数进行归一化，得到单位化的x进入下一轮循环
    if norm(r)<=(norm_of_A+abs(lambda))*(10^-16)%
        % 10^-16是机器精度，这个tolerance是相对A与lambda而变化的
        norm(r) %显示一下r的2-范数
        disp("r足够小，结束迭代")
        flg=1;
        break;
    end
end
if flg==0
    disp("触发迭代保护而退出")
end
plot([1:count],lambda_list(1:count,1),'b--o');
title("用乘幂法求A的谱半径")
xlabel("迭代次数")
ylabel("瑞利商（特征值的近似值）")


