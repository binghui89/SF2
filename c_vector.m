function [c1, c2] = c_vector(copula_discrete, k)
% Construct vector multiplier for piecewise convolution. The returned c
% vectors are used in the convolution over interval 
% (x^c_1 + (k-1) ?x, x^c_1 + k ?x)
% Input:
% copula_discrete is the discretized copula, size: n x m (m >= n). 1 <= k <= (m + n).
% Note m is length of the first dimension (horizontal axis), and n is 
% length of the second dimension (vertical axis). 
% Output:
% c1 is a (m x 1) vector, c2 is a (n x 1) vector.

if size(copula_discrete, 2) < size(copula_discrete, 1)
    copula_discrete = copula_discrete';
    flag_t = true;
else
    flag_t = false;
end

[n, m] = size(copula_discrete);

c1 = zeros(m, 1); % First dimension, horizontal direction
c2 = zeros(n, 1); % Second dimension, vertical direction

if k <= n
    c1(1) = 1;
    c2(k) = copula_discrete(k, 1);
    for i = 2:k
        c1(i) = c1(i-1)*copula_discrete(k-i+1, i)/copula_discrete(k-i+1, i-1);
        c2(k-i+1) = c2(k-i+2)*copula_discrete(k-i+1, i-1)/copula_discrete(k-i+2, i-1);
    end
elseif (k > n) && (k <= m)
    c1(k-n) = 1;
    c2(n) = copula_discrete(n, k-n);
    for i = 1: n-1
        c1(k-n+i) = c1(k-n+i-1)*copula_discrete(n-i+1, k-n+i)/copula_discrete(n-i+1, k-n+i-1);
        c2(n-i) = c2(n-i+1)*copula_discrete(n-i, k-n+i)/copula_discrete(n-i+1, k-n+i);
    end
    i = n;
    c1(k-n+i) = c1(k-n+i-1)*copula_discrete(n-i+1, k-n+i)/copula_discrete(n-i+1, k-n+i-1);
else % k > m
    c1(k-n) = 1;
    c2(n) = copula_discrete(n, k-n);
    for i = 1: m+n-k
        c1(k-n+i) = c1(k-n+i-1)*copula_discrete(n-i+1, k-n+i)/copula_discrete(n-i+1, k-n+i-1);
        c2(n-i) = c2(n-i+1)*copula_discrete(n-i, k-n+i)/copula_discrete(n-i+1, k-n+i);
    end
end

if flag_t
    tmp = c1;
    c1 = c2;
    c2 = tmp;
end

end