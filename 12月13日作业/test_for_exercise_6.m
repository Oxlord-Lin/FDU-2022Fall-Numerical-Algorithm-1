n = 1e3;
S = diag([1:n]);
for i = 1:n-1
    for j = i+1:n
        S(i,j) = 5000*rand();
    end
end

V = eye(n);
V(n,1) = 1e-6;
U = orth(V);
A = U*S*U';
[U,D] = eig(A);
[~,ind] = sort((diag(D)));
D = D(ind,ind);
U = U(:,ind);
plot(abs(diag(D)));
figure()
plot(real(diag(D)));
figure()
R = real(diag(D));
I = imag(diag(D));
scatter3(R,I,[1:n]);