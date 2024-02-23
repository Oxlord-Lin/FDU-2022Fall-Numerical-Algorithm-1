n = 10;
% 构造爪形矩阵A
A = diag(randn(n,1));
for i = 1:n
    A(1,i) = randn();
    A(i,1) = randn();
end
A
% form the Hessenberg form of A using Givens rotation
for i = n:-1:3
    c = A(i-1,1);
    s = A(i,1);
    len = norm([c,s]);
    c = c/len;
    s = s/len;
    Q = [c,s;-s,c];
    A(i-1:i,:) = Q*A(i-1:i,:);
end
A
for i = 1:n-1
    c = A(i,i);
    s = A(i+1,i);
    len = norm([c,s]);
    c = c/len;
    s = s/len;
    Q = [c,s;-s,c];
    A(i:i+1,i:n) = Q*A(i:i+1,i:n);
end


