%% the function "upper_hessenberg" using householder orthogonalization to make matrix A an upper Hessenberg matrix 
% 最后返回A，U，满足：H=U'AU，其中H是上海森伯格阵
function [A,U] = upper_hessenberg(A)
n = size(A);
n = n(1,1);
U = eye(n);
for i = 1:n-2
    u = A(i+1:n,i);
    u(1,1) = u(1,1) + sign(u(1,1)) * norm(u);
    if norm(u)==0
        continue
    end
    w = [zeros(i,1);u];
    Q = eye(n)-(2/(w'*w))*(w*w');  %generate the Householder orthogonalizer
    A = Q*A*Q';
    U = U*Q';
end
end

