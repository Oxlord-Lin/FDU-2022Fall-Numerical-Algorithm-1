n = 10;
R = randn(n);
for i = 1:n
    for j = 1:i-1
        R(i,j) = 0;
    end
end
u = rand(n,1);
v = rand(n,1);
R0 = R + u*v';
% 从下到上通过n-1个Givens rotation将u变成±norm(u)*e1
for i = n-1:-1:1
    c = u(i,1);
    s = u(i+1,1);
    d = norm([c,s]);
    c = c/d;
    s = s/d;
    Q = [c,s;-s,c];
    u(i:i+1,1) = Q*u(i:i+1,1);
    R(i:i+1,i:n) = Q*R(i:i+1,i:n);
end
R = R + u*v';
for i = 1:n-1
    c = R(i,i);
    s = R(i+1,i);
    d = norm([c,s]);
    c = c/d;
    s = s/d;
    Q = [c,s;-s,c];
    R(i:i+1,i:n) = Q*R(i:i+1,i:n);
end
imagesc(abs(R))
colorbar
