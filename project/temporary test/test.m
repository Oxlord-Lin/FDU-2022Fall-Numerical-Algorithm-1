T=[-1.0194    0.0017    0.0000    0.0000   -0.0000    0.0000
    0.0017    1.1086    0.0050    0.0000    0.0000    0.0000
   -0.0000    0.0050   -0.9234    0.0000   -0.0000   -0.0000
   -0.0000    0.0000    0.0000    0.8280    0.0000   -0.0000
   -0.0000   -0.0000   -0.0000    0.0000    0.4482    0.0009
    0.0000   -0.0000    0.0000   -0.0000    0.0009    0.5223
];
%%
[A,G] = one_step_QR_with_Wilkinson_shift(A);
A;
%%
% 返回A与G_accumulated，    满足关系：A(output) = G_accumulated' * A(input) * G_accumulated

%only for test
% n = size(A);
% n = n(1,1);
% [q,r] = qr(A-A(n,n)*eye(n));
% A = r*q +A(n,n)*eye(n)


%% 带Wilkinson位移的隐式QR迭代
function [A,G_accumulated] = one_step_QR_with_Wilkinson_shift(A)
% 返回A与G_accumulated，满足关系：A(output) = G_accumulated' * A(input) * G_accumulated
n = size(A);
n = n(1,1);
sigma = (A(n-1,n-1)-A(n,n))/2;
Wilkinson_shift = A(n,n) - (A(n,n-1)^2) / (sigma + sign(sigma) * sqrt(sigma^2 + A(n,n-1)^2));
x = A(1,1)-Wilkinson_shift;
y = A(2,1);
G_accumulated = eye(n);
for k = 1:n-1
    %生成Givens旋转中的角度
    %[c,s] = givens(x,y);
    if y == 0
        c = 1; s = 0;
    else
        if abs(y)>abs(x)
            tao = x/y; s = 1/sqrt(1+tao^2); c = s*tao;
        else
            tao = y/x; c=1/sqrt(1+tao^2); s=c*tao;
        end
    end
    %生成Givens旋转矩阵
    G = eye(n); G(k,k) = c; G(k+1,k+1) = c; G(k,k+1) = -s; G(k+1,k) = s;
    %驱赶"气泡"
    A = G'*A*G
    G_accumulated = G_accumulated * G;
    %为下一次迭代的Givens矩阵所需参数做准备
    if k<n-1
    x = A(k+1,k);
    y = A(k+2,k); 
    end
end
end

