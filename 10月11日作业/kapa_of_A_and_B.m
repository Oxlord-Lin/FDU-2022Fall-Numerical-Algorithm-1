kapa_A = zeros(100,1);
kapa_B = zeros(100,1);
for k = 1:100
    n = 10*k;
    A = generate_A(n);
%     kA = Kapa(A) 
    kapa_A(k,1) = Kapa(A);
    B = generate_B(n);
%     kB = Kapa(B)
    kapa_B(k,1) = Kapa(B);
end
plot([1:100],log10(kapa_A));
title('矩阵A的2-条件数（取对数）')
figure
plot([1:100],log10(kapa_B));
title('矩阵B的2-条件数（取对数）')


function A = generate_A(n)
A = diag(8*ones(n,1)) + diag(ones(n-1,1),1) + diag(6*ones(n-1,1),-1);
end

function B = generate_B(n)
B = diag(6*ones(n,1)) + diag(ones(n-1,1),1) + diag(8*ones(n-1,1),-1);
end

function k = Kapa(M)
k = max(svd(M))/min(svd(M));
end