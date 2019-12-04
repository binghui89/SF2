function [anew, tnew] = cleanup_poly(a, t)
% To clean up trailing ones and leading zeros in the CDF
[I, M] = size(a);
a_int = cat(2, zeros(I, 1), a./repmat(1:M, I, 1)); % Coefficients of integrated function, constant term is 0
cdf_results = zeros(size(t));
for i = 2: size(t, 1)
    x1 = t(i-1);
    x2 = t(i);
    cdf_results(i) = cdf_results(i-1) + polyval(fliplr(a_int(i-1, :)), x2) - polyval(fliplr(a_int(i-1, :)), x1);
end

tmp1 = cumsum(abs(cdf_results));
ind_tmp1_eq_0 = find(tmp1==0);
istart = ind_tmp1_eq_0(end); % The last zero is the starting point

cdf_diff = diff(cdf_results);
tmp2 = cumsum(abs(flipud(cdf_diff)));
ind_tmp2_eq_0 = find(tmp2==0);
iend = I - ind_tmp2_eq_0(end); % The end point is where CDF becomes constant

anew = a(istart:iend, :);
tnew = t(istart:iend, :);

end