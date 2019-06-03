function cdf = cdf_poly(a, t)
[I, M] = size(a);
a_int = cat(2, zeros(I, 1), a./repmat(1:M, I, 1)); % Coefficients of integrated function, constant term is 0
cdf = zeros(size(t));
for i = 2: size(t, 1)
    x1 = t(i-1);
    x2 = t(i);
    cdf(i) = cdf(i-1) + polyval(fliplr(a_int(i-1, :)), x2) - polyval(fliplr(a_int(i-1, :)), x1);
end
cdf(end) = 1;
cdf(cdf>1)=1;
end