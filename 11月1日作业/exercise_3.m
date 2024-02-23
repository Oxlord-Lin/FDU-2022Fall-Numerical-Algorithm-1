%% 先随机产生一个实对称矩阵
Q = orth(rand(200,200));
D = diag([1:200]);
A = Q*D*Q';
%% 选择lambda=69作为试验对象
offset=0.01;
[l,u,p] = lu(A-(69+offset)*eye(200));%距离69+offset最近的特征值是69
N=norm(inv(A-(69+offset)*eye(200)));
x = rand(200,1);
x = x/norm(x);%归一化
flg=0;
count=0;
lambda_list=zeros(10^4,1);
for i = 1:(10^4)%迭代保护
    count=count+1;
    y = u\(l\(p'*x));
   %y=(A-(69+offset)*eye(200))\x;
    t=(y'*x)/(x'*x);%瑞利商
    lambda=1/t+(69+offset);
    lambda_list(count,1)=lambda;
    r=y-t*x;%计算残差
    x=y/norm(y);%归一化，以备进入下一轮循环
    if norm(r)<=(N+abs(t))*(1e-16)
        disp("残差r足够小，结束反幂法迭代")
        flg=1;
        break
    end
end

%% 绘图
% if flg==0
%     disp("触发迭代保护")
% end
% if flg==1
%     plot([1:count],lambda_list(1:count,1),'b--o');
%     title("使用反幂法计算特征值的收敛过程")
%     xlabel("迭代次数")
%     ylabel("特征值的近似值")
% end

