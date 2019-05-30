function [a_convcorr, t_convcorr] = conv_poly_corr(a1, t1, a2, t2, copula_discrete, delta_t)
% copula_discrete's x direction corresponds to a1 and t1, y direction 
% corresponds to a2 and t2

tol = 1e-6; % Tolerance

% With correlation
% deltat = unique(diff(t1)); % Samples must be equally placed
deltat = delta_t; % Samples must be equally placed
tc_1 = t1(1)+t2(1); % Starting point of (t1 + t2)
tc_e = t1(end)+t2(end); % Ending point of (t1 + t2)
K = round((tc_e - tc_1)/deltat); % Number of discretized intervals

t_convcorr = nan(K+1, 1);
a_convcorr = nan(K+1, size(a1, 2) + size(a2, 2));
t_convcorr(1) = tc_1;
a_convcorr(end, :) = 0;

for k = 1: K
    [cx, cy] = c_vector(copula_discrete, k);
    ax = repmat([cx; 0], 1, size(a1, 2)).*a1;
    ay = repmat([cy; 0], 1, size(a2, 2)).*a2;
    [a_tmp, t_tmp] = conv_poly(ax, t1, ay, t2, delta_t);
    try
        a_convcorr(k, :) = a_tmp(t_tmp==t_convcorr(k), :);
    catch
        disp('error\n');
    end
%     a_convcorr(k, :) = a_tmp(abs(t_tmp-t_convcorr(k))<tol*deltat, :);
    t_convcorr(k+1) = tc_1+k*deltat;
end

end