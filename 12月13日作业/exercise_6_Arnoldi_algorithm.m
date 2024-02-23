n = 1e3;
%生成非对称的A
% lambda = sort(randn(1,n));
% S = diag(lambda);
% for i = 1:n-1
%     for j = i+1:n
%         S(i,j) = (1e-1)*rand();
%     end
% end
S = randn(n);
S = triu(S);
lambda = sort(diag(S));
plot(lambda);
title('eigenvalues of A')
U = orth(randn(n));
A = U*S*U';
Q = zeros(n);
H = zeros(n);
v = rand(n,1);
v = v./norm(v);
Q(1:n,1) = v;
for i = 1:n-1
    y = A*Q(1:n,i); % y=A*q_i
    for j = 1:i
        H(j,i) = Q(1:n,j)'*y; % 向前i个正交基向量上投影
        y = y - H(j,i)*Q(1:n,j);
    end
    H(i+1,i) = norm(y);
    if abs(H(i+1,i)) < 1e-16 %若H(i,i+1)足够小
        flg = 1;
        disp("提前终止，因为H(i,i+1)足够小，认为r已经完全落在i维Krylov子空间中")
        break
    end
    Q(1:n,i+1) = y./H(i+1,i); %生成第i+1个正交基向量
end
% H的最后一列单独生成
y = A*Q(1:n,n-1); % y=A*q_i
for j = 1:n
    H(j,n) = Q(1:n,j)'*y; % 向前i个正交基向量上投影
    y = y - H(j,n)*Q(1:n,j);
end


[V,D] = eig(H);
eigenvalues = sort(diag(D));

figure()
plot(lambda,'--*')
hold on
plot(sort(real(eigenvalues)),'r');
title("用Arnoldi过程计算得到的Heseenberg阵的特征值的实部");

% figure()
% plot(abs([(eigenvalues)-[1:n]']))
% title("abs([(eigenvalues of H)-(eigenvalues of A)])");

% figure()
% plot(real([(eigenvalues)-[1:n]']));
% title("real([(eigenvalues of H)-(eigenvalues of A)])");
%%
figure()
eigenvalues = sort(eigenvalues,'ComparisonMethod','real') ;
R = real(eigenvalues);
I = imag(eigenvalues);
scatter3(lambda,zeros(1,n),[1:n],'b','*')
hold on
scatter3(R,I,[1:n],'r');
xlabel('Real')
ylabel('Image')
zlabel('index')
title("用Arnoldi过程计算得到的Hessenberg阵的特征值的实部与虚部");