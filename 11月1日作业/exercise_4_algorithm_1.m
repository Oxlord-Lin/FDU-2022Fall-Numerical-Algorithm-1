n = 20;
M=diag(ones(n,1))+diag(ones(n-1,1),1);
M(n,1)=1;
M=0.5*M;
x = -0.5+rand(n,1);
x = x/norm(x);
y = -0.5+rand(n,1);
y = y/norm(y);
figure();%初始状态
plot(x,y,'b-o');
lim=0.5;
xlim([-lim,lim]);
ylim([-lim,lim]);
title(["迭代次数:",0]);
for k = 1:200
    x=M*x;
    y=M*y;
    if k==5 || k==20 || k==100 ||k==200
        figure();
        plot(x,y,'b-o')
        xlim([-lim,lim]);
        ylim([-lim,lim]);
        title(["迭代次数:",k]);
    end
end