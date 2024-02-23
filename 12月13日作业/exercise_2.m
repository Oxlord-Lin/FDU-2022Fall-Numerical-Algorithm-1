%% 用于求v'*f(A)*v的Lanczos算法
n = 1e3;
D1 = diag(randn(1,n));
U1= orth(randn(n));
A = U1*D1*U1';
v = randn(n,1);
length = norm(v);
Q = zeros(n,n);
Q(1:n,1) = v./length;  %v单位化后作为Q的第一列
T = zeros(n,n);
convergence_History_1 = zeros(1,n+1);
% convergence_History_2 = zeros(1,n+1);
convergence_History_3 = zeros(1,n+1);
convergence_History_4 = zeros(1,n+1);
convergence_History_1(1,n+1) = v'*expm(A)*v;
% convergence_History_2(1,n+1) = v'*logm(A)*v; 
convergence_History_3(1,n+1) = v'*poly(A)*v; 
convergence_History_4(1,n+1) = v'*funm(A, @sin)*v; 

for i = 1:n-1
    y = A*Q(1:n,i); % y=A*q_i
    T(i,i) = Q(1:n,i)'*y;
    y = y - T(i,i)*Q(1:n,i);
    if i-1>=1
        y = y - T(i,i-1)*Q(1:n,i-1);
    end
    T(i+1,i) = norm(y);
    if  T(i+1,i)<1e-16
        disp('提前结束');
        break
    end
    T(i,i+1) = T(i+1,i);
    Q(1:n,i+1) = y./T(i+1,i); %生成第i+1个正交基向量
%   计算v'*f(T_i阶子阵))*v作为v'*f(A)*v的近似
    e = [1,zeros(1,i-1)]';
    convergence_History_1(1,i) = length^2*e'*expm(T(1:i,1:i))*e;
%     convergence_History_2(1,i) = length^2*e'*logm(T(1:i,1:i))*e;
    convergence_History_3(1,i) = length^2*e'*poly(T(1:i,1:i))*e; 
    convergence_History_4(1,i) = length^2*e'*funm(T(1:i,1:i), @sin)*e; 

end
i = n;
e = [1,zeros(1,i-1)]';
convergence_History_1(1,i) = length^2*e'*expm(T(1:i,1:i))*e;
% convergence_History_2(1,i) = length^2*e'*logm(T(1:i,1:i))*e;
convergence_History_3(1,i) = length^2*e'*poly(T(1:i,1:i))*e; 
convergence_History_4(1,i) = length^2*e'*funm(T(1:i,1:i), @sin)*e; 

figure()
plot(convergence_History_1);
title("convergence history of length^2*e'*expm(T(1:i,1:i))*e");
xlabel('iter');
ylabel('value');
% figure()
% plot(convergence_History_2);
figure()
plot(convergence_History_3);
title("convergence history of length^2*e'*poly(T(1:i,1:i))*e");
xlabel('iter');
ylabel('value');
figure()
plot(convergence_History_4);
title("convergence history of length^2*e'*funm(T(1:i,1:i), @sin)*e");
xlabel('iter');
ylabel('value');

function B = poly(A)
B = 2*A*A*A-8.5*A*A+(2/3)*A;
end
