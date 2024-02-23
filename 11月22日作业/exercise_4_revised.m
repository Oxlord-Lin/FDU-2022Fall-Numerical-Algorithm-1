%% test the algorithm
max_n=100;
depature = zeros(1,max_n);
for n = 1:max_n
    U = orth(randn(n,n));
    lambda = diag(rand(1,n)); % 这样的条件数是比较良态，而且各阶矩阵均的kapa相近，适合作图比较
%     lambda = diag([1:n]); % 这个的条件数会迅速增大，即便用到了前5项截尾Taylor，误差仍然巨大
    A = U*lambda*U';
    exp_A_by_function = matrix_exponential(A);
    exp_A = U*expm(lambda)*U'; %注意，这里要用expm，而不能用exp
    depature(1,n) = norm(exp_A_by_function-exp_A);
end
figure()
plot([1:max_n],log10(depature),'-*');
xlabel("矩阵阶数")
ylabel("log10(depature)")
title('误差（取对数） Taylor展开保留前两项，特征值为[1:n]');
%% implement the scaling_and_squaring algorithm for computing the matrix exponential
% (combined with truncated Taylor series)
function exp_M = matrix_exponential(M)
n = size(M);
n = n(1,1);
e = eig(M);
max_eigenvalue = max(e);
k = ceil(log2(1e5*max_eigenvalue)); 
% 1e3是试验下来与前5项截尾泰勒比较相合的一个缩放因子
% 如果只用前三项截尾泰勒，那么取1e5会更好
% 这里的缩放系数越大，理论上用taylor来逼近也会更准确
% 但是，由于之后还要求exp_M_scaled的幂，因而也不能缩的太小，否则之后会由于矩阵乘法过程中的舍入误差而得不偿失
M_scaled = M/(2^k);
exp_M_scaled = eye(n) + M_scaled + (1/2)*M_scaled*M_scaled %+ (1/6)*M_scaled*M_scaled*M_scaled + (1/24)*M_scaled*M_scaled*M_scaled*M_scaled;% truncated Taylor series
% 这个截尾泰勒展开一共五项，用前三项效果已经不错
exp_M = exp_M_scaled;
for i = 1:k
    exp_M = exp_M*exp_M;
end
% exp_M = eye(n) + M_scaled + (1/2)*M*M + (1/6)*M*M*M;
end