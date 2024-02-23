%% test the algorithm
max_n=50;
depature = zeros(1,max_n);
for n = 1:max_n
    U = orth(rand(n,n));
    A = U*diag(10*rand(1,n))*U';
    [u1,t1] = schur(A);
    exp_A_by_function = u1*Schur_Parlett_algorithm_for_exp_A*u1';
    exp_A = U*exp(diag(1:n))*U';
    depature(1,n) = norm(exp_A_by_function-exp_A,'fro');
end
figure()
plot([1:max_n],log10(depature),'--o');
xlabel("矩阵阶数")
ylabel("log10(depature)")
[t,s] = title('误差（取对数）','矩阵阶数：1~50')

%% Schur-Parlett algorithm (简化版，主要练习求解Sylvester方程)
function exp_A = Schur_Parlett_algorithm_for_exp_A(A)

end









end