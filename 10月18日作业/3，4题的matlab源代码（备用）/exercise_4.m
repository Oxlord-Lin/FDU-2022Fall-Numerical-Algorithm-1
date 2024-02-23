n=1:6;
y=(exp(-n/2))';
%绘制散点图
scatter(n,y,"filled");
hold on
%用最小二乘法求斜率和截距
A=[n',ones(6,1)];
x=pinv(A)*y;
k=x(1,1);
b=x(2,1);
%绘制对应直线
Y=k*n+b;
plot(n,Y);

