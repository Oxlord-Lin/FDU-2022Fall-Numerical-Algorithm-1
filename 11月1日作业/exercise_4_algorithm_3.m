n = 20;
M=diag(ones(n,1))+diag(ones(n-1,1),1);
M(n,1)=1;
M=0.5*M;
for j =1:n
    tao(j,1)=2*pi/n*(j-1);
end
c=sqrt(2/n)*cos(tao);
s=sqrt(2/n)*sin(tao);
theta1=10*rand();
theta2=10*rand();
x = cos(theta1)*c + sin(theta1)*s;
y = cos(theta2)*c + sin(theta2)*s;
figure();%初始状态
plot(x,y,'b--o');
lim=0.4;
xlim([-lim,lim]);
ylim([-lim,lim]);
hold on
for k = 1:10
    x=M*x;
    x=x/norm(x);
    y=M*y;
    y=y/norm(y);
    if mod(k,2)==0
        scatter(x,y,'b','filled')
        hold on
    end
    if mod(k,2)==1
        scatter(x,y,'r','filled')
        hold on
    end
end
title("蓝色是偶数次迭代结果，红色是奇数次迭代结果")