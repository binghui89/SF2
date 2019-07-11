function [a, t] = reformulate_poly(a0, t0)
% Reformulate spline PDF to reduce numerical errors
% Now we just convert any degree of spline PDF into 0 degree spline PDF
cdf0 = cdf_poly(a0, t0);
if max(cdf0) < 1
    warning('Maximum CDF < 1!');
end

diff_cdf0 = diff(cdf0);
diff_t0   = diff(t0);
if any(diff_cdf0 < 0)
    warning('Negative pdf detected!');
end

% a = diff_cdf0(diff_cdf0>0)./diff_t0(diff_cdf0>0);

i_start = find(cdf0<1e-4, 1, 'last' );
i_end   = find(abs(cdf0-1)<1e-4, 1, 'first');

cdf_valid = cdf0(i_start: i_end);
t_valid   = t0(i_start: i_end);

a = diff(cdf_valid)./diff(t_valid);
a = [a;0];
t = t_valid;
end