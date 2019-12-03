function [a_convcorr, t_convcorr] = conv_poly_corr(a1, t1, a2, t2, copula_discrete, delta_t)
% Piecewise polynomial convolution with correlation
% USE: [a_convoluted, t_convoluted] = CONV_POLY_CORR(a1, t1, a2, t2, copula_discrete, delta_t)
%
% INPUT:
% a1, a2: 2-dim array of coefficients of piecewise polynomial functions. 
% Each row includes a segment, and the coefficients are sorted by degree 
% in an ascending order, starting from 0, i.e., constant terms, to the 
% highest degree term (M - 1).
% Size: I by M, I is number of segments + 1 (last segment is always 0), 
% M is the highest order + 1 (lowest order is always 0).
%
% t1, t2: column vectors of starting points of all segments.
% Size: I by 1, the first element represents the starting point of the
% first non-zero segment (the first segment is always 0), and the last
% element represents the starting point of the last segment (always 0), or
% the ending point of the last non-zero segment.
%
% copula_discrete: 2-dim array of average copula PDF over the discretized
% grid by t1 and t2.
% Size: (I2 - 1) by (I1 - 1). Because a1, t1 is the x direction function,
% which is row-wise and the 2nd dimension, and a2, t2 is the y direction
% function, which is column-wise and the 1st dimension. I1 and I2 are
% number of all points, so the number of segments are (I1 - 1) and (I2 -
% 1), respectively.
% 
% delta_t: segment size of both piecewise polynomial functions, i.e.,
% difference of any two adjacent elements in t1 and t2. 
% Size: Scaler.
% 
% OUTPUT:
% a_convoluted: 2-dim array of coefficients of the convoluted piecewise
% polynomial function in the same format as a1 and a2.
%
% t_convoluted: column vector of starting points of all segments of the
% convoluted piecewise polynomial function.

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
    [cx, cy] = copula_decompose(copula_discrete, k);
    a_sumtmp = zeros(1, size(a_convcorr, 2));
    for j = 1: size(cx, 2)
        ax = repmat([cx(:, j); 0], 1, size(a1, 2)).*a1;
        ay = repmat([cy(:, j); 0], 1, size(a2, 2)).*a2;
        [a_tmp, t_tmp] = conv_poly(ax, t1, ay, t2, delta_t);
        try
            a_sumtmp = a_sumtmp + a_tmp(t_tmp==t_convcorr(k), :);
        catch
            disp('error\n');
        end
    end
%     a_convcorr(k, :) = a_tmp(abs(t_tmp-t_convcorr(k))<tol*deltat, :);
    a_convcorr(k, :) = a_sumtmp;
    t_convcorr(k+1) = tc_1+k*deltat;
end

end