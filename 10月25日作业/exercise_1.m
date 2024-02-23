%% 构造矩阵A
m = 80;
n = 50;
lambda = randn();
A = [zeros(m,n);lambda*eye(n)];
Q = eye(m+n);
for i = 1:n-1
    A(i,i) = randn();
    A(i,i+1) = randn();
end
A(n,n) = randn();
for i = 1:n-1
    c = A(i,i);
    s = A(m+i,i);
    len = norm([c,s]);
    c = c/len;
    s = s/len;
    Q_temp1 = Givens(m+n,i,m+i,c,s);
    A = Q_temp1*A; 
    Q = Q_temp1*Q;
    c2 = A(m+i,i+1);
    s2 = A(m+i+1,i+1);
    len2 = norm([c2,s2]);
    c2 = c2/len2;
    s2 = s2/len2;
    Q_temp2 = Givens(m+n,m+i,m+i+1,s2,-c2); % 注意这一步是要把下半部分多出来的那个元素消掉
    A = Q_temp2*A;
    Q = Q_temp2*Q;
end
% 然后要把最右下角的元素消掉
c = A(n,n);
s = A(m+n,n);
len = norm([c,s]);
c = c/len;
s = s/len;
Q_temp1 = Givens(m+n,n,m+n,c,s);
A = Q_temp1*A; 
Q = Q_temp1*Q;

%%
function Q = Givens(p,index1,index2,c,s)
Q = eye(p);
Q(index1,index1) = c;
Q(index2,index2) = c;
Q(index1,index2) = s;
Q(index2,index1) = -s;
end