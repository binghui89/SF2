function [x_conv, p_conv] = conv_ez_corr(cell_bincenter, cell_p, bin_width, copula_pdf_nd, cell_lu)
% Copula enhanced multidimensional discrete convolution up to 5-dim
% USE: [x_conv, p_conv] = CONV_EZ_CORR(cell_bincenter, cell_p, bin_width, copula_pdf_nd)
%
% INPUT:
% cell_bincenter: N-dim cell array of bin centers of all convolved
% variables. The ith cell is a column vector of M_i by 1.
%
% cell_p: N-dim cell array of probabilitis of each bin. The ith cell is a 
% column vector of M_i by 1.
%
% bin_width: Bin width, should be equal across all variables. 
%
% copula_pdf_nd: N-dim array of discretized mean copula PDF. Size is M1 x
% M2 x M3 x ... x M_N.
% 
% OUTPUT:
% x_conv: column vector of bin centers of the convolved variable. Bin size
% is the same as bin_width.
%
% p_conv: column vector of probabilities of all bins of the convolved
% variable.

N = ndims(copula_pdf_nd);
ind = allsub_nd(copula_pdf_nd);
sumind = sum(ind, 2);

flat_p_alldim = zeros(numel(copula_pdf_nd), N);
for i = 1: N
    flat_p_alldim(:, i) = cell_p{i}(ind(:, i));
end
flat_copulapdf_and_p = [copula_pdf_nd(:), flat_p_alldim];

%%%%%%%%%%%%%%%%%%%%
% Find the index of the cap and floor
if nargin == 4
    cell_lu = cell(N, 1);
end

for i = 1: N
    if isempty(cell_lu{i})
        continue;
    else
        lb = cell_lu{i}(1);
        ub = cell_lu{i}(2);
        i_lb = max(find(cell_bincenter{i}<lb))+1;
        i_ub = min(find(cell_bincenter{i}>ub))-1;
        if ~isempty(i_lb)
            ind(ind(:, i)<i_lb, i) = i_lb;
        end
        if ~isempty(i_ub)
            ind(ind(:, i)>i_ub, i) = i_ub;
        end
    end
end
% ind(ind>2) = 2;
sumind = sum(ind, 2);
%%%%%%%%%%%%%%%%%%%%
clear ind copula_pdf_nd flat_p_alldim;

ar_delta = [N: max(sumind)]'; % Array of number of deltas in the convolved variable
x_conv0 = 0;
for i = 1:N
    x_conv0 = x_conv0 + cell_bincenter{i}(1) - bin_width;
end
x_conv = x_conv0 + ar_delta.*bin_width; % Bin center of convolved results
p_conv = zeros(numel(ar_delta), 1); % Probability of each bin in the results
for k = 1: numel(ar_delta)
    thissum = ar_delta(k);
    disp(thissum);
    p_conv(k) = sum(prod(flat_copulapdf_and_p(sumind==thissum, :), 2), 1);
end
end

function sub = allsub_nd(A)
% Now only supports up to 5-dim array
N = ndims(A);

sub = zeros(numel(A), N);

switch N
    case 2
        [sub(:, 1), sub(:, 2)] = ind2sub(size(A), [1: numel(A)]');
    case 3
        [sub(:, 1), sub(:, 2), sub(:, 3)] = ind2sub(size(A), [1: numel(A)]');
    case 4
        [sub(:, 1), sub(:, 2), sub(:, 3), sub(:, 4)] = ind2sub(size(A), [1: numel(A)]');
    case 5
        [sub(:, 1), sub(:, 2), sub(:, 3), sub(:, 4), sub(:, 5)] = ind2sub(size(A), [1: numel(A)]');
end

end