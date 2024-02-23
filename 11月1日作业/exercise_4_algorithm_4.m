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
u = cos(theta1)*c + sin(theta1)*s;
v = cos(theta2)*c + sin(theta2)*s;
u = u/norm(u,1);
v = v/norm(v,"inf");

figure();%初始状态
plot(u,v,'b:o');
lim=0.4;
hold on
for k = 1:5
    u=M*u;
    v=M*v;
    u = u/norm(u,1);
    v = v/norm(v,"inf");
    if mod(k,2)==0
        scatter(u,v,'b','filled')
    end
    if mod(k,2)==1
        scatter(u,v,'r','filled')
    end
end
title("蓝色是偶数次迭代结果，红色是奇数次迭代结果，u使用1范数，v使用无穷范数")