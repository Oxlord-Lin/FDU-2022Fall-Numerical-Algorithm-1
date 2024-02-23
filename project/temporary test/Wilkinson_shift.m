function shift = Wilkinson_shift(A)
n = size(A);
n = n(1,1);
sigma = (A(n-1,n-1)-A(n,n))/2;
shift = A(n,n) - (A(n,n-1)^2) / (sigma + sign(sigma) * sqrt(sigma^2 + A(n,n-1)^2));
end