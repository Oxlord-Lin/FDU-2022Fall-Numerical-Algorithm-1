n = 40;
M=diag(ones(n,1))+diag(ones(n-1,1),1);
M(n,1)=1;
M=0.5*M;
x = rand(n,1);
x = x/norm(x);
y = rand(n,1);
y = y/norm(y);
x=x-sum(x)/n;
y=y-sum(y)/n;
figure();%初始状态
plot(x,y,'b-o');
lim=0.5;
xlim([-lim,lim]);
ylim([-lim,lim]);
title(["迭代次数:",0]);
for k = 1:800
    x=M*x;
    x=x/norm(x);
    y=M*y;
    y=y/norm(y);
    if k==5 || k==20 || k==100 ||k==200||k==400||k==800
        figure();
        plot(x,y,'b-o')
        xlim([-lim,lim]);
        ylim([-lim,lim]);
        title(['迭代次数:',num2str(k)]);
    end
end