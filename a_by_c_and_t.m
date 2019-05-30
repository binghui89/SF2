function a = a_by_c_and_t(c, t)
[I, M] = size(c);
[m, l] = meshgrid(1:M, 1:M);
l_minus_m = tril(l - m);
e2_part = tril(1./...
    (factorial(m-1).*factorial(l_minus_m)) ...
    );
a = nan(size(c));
for i = 1: I
    if i == 1
        a(i, :) =  c(i, :)*tril(e2_part.*(-t(i)).^l_minus_m);
    else
        a(i, :) =  a(i-1, :) + c(i, :)*tril(e2_part.*(-t(i)).^l_minus_m);
    end
end
end