function c = c_by_a_and_t(a, t)
[I, M] = size(a);
[m, l] = meshgrid(1:M, 1:M);
l_minus_m = tril(l - m);
e1_part = tril(...
    factorial(l-1)./factorial(l_minus_m) ...
    );
c = nan(size(a));
for i = 1: I
    if i == 1
        c(i, :) =  a(i, :)*tril(e1_part.*t(i).^l_minus_m);
    else
        c(i, :) =  (a(i, :) - a(i-1, :))*tril(e1_part.*t(i).^l_minus_m);
    end
end
end