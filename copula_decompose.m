function [c1, c2] = copula_decompose(copula_discrete, k)
% Construct vector multiplier for piecewise convolution. The returned c
% vectors are used in the convolution over interval 
% (x^c_1 + (k-1) \delta x, x^c_1 + k \delta x)
% Input:
% copula_discrete is the discretized copula, size: n x m (m >= n). 1 <= k <= (m + n).
% Note m is length of the first dimension (horizontal axis), and n is 
% length of the second dimension (vertical axis). 
% Output:
% c1 is a (m x 1) vector, c2 is a (n x 1) vector.

% If m < n, transpose it
if size(copula_discrete, 2) < size(copula_discrete, 1)
    copula_discrete = copula_discrete';
    flag_t = true; % Flag of transpose
else
    flag_t = false;
end

[n, m] = size(copula_discrete);

copula_left_sub  = [min(k, n): -1: max(k-m, 1); max(k-n, 0): 1: min(m, k-1)]';
copula_right_sub = copula_left_sub;
copula_right_sub(:, 2) = copula_right_sub(:, 2) + 1;
copula_left_sub(any(copula_left_sub==0, 2), :) = nan; % For column subscripts less than 1
copula_right_sub(any(copula_right_sub>m, 2), :) = nan; % For column subscripts greater than m
copula_left_ind = sub2ind([n, m], copula_left_sub(:, 1), copula_left_sub(:, 2));
copula_right_ind = sub2ind([n, m], copula_right_sub(:, 1), copula_right_sub(:, 2));
copula_ind = reshape([copula_left_ind copula_right_ind]', numel(copula_left_ind)+numel(copula_right_ind), 1);
copula_ind(isnan(copula_ind))=[];

% Detect separating point: zeros with both non-zeros on the left and right
copula_to_decompose = copula_discrete(copula_ind);
sep_logical = [false;copula_to_decompose(1:end-1)~=0] & [copula_to_decompose(2:end)~=0;false] & (copula_to_decompose==0);
sep_logical_cumsum = cumsum(sep_logical);
copula_ind_even = copula_ind(mod(sep_logical_cumsum, 2) == 0);
copula_ind_odd  = copula_ind(mod(sep_logical_cumsum, 2) == 1);

% Separate the original discrete copula matrix into two matrices that will
% not lead to confliction due to 0
% copula_discrete_even = zeros(n, m);
% copula_discrete_even(copula_ind_even) = copula_discrete(copula_ind_even);
% copula_discrete_odd = zeros(n, m);
% copula_discrete_odd(copula_ind_odd) = copula_discrete(copula_ind_odd);

if isempty(copula_ind_odd)
    cell_copula_select = cell(1, 1);
    copula_discrete_even = zeros(n, m);
    copula_discrete_even(copula_ind_even) = copula_discrete(copula_ind_even);
    cell_copula_select{1} = copula_discrete_even;
else
    cell_copula_select = cell(2, 1);
    copula_discrete_even = zeros(n, m);
    copula_discrete_even(copula_ind_even) = copula_discrete(copula_ind_even);
    cell_copula_select{1} = copula_discrete_even;
    copula_discrete_odd = zeros(n, m);
    copula_discrete_odd(copula_ind_odd) = copula_discrete(copula_ind_odd);
    cell_copula_select{2} = copula_discrete_odd;
end
c1 = zeros(m, length(cell_copula_select)); % First dimension, horizontal direction
c2 = zeros(n, length(cell_copula_select)); % Second dimension, vertical direction

for r = 1: length(cell_copula_select)
    % Find non-zero blocks
    copula_select = cell_copula_select{r};
    copula_select_vector = copula_select(copula_ind);
    start_logical = [true; copula_select_vector(1:end-1)==0] & (copula_select_vector~=0);
    end_logical = [copula_select_vector(2:end)==0; true;] & (copula_select_vector~=0);
    istart = find(start_logical);
    iend   = find(end_logical);
    assert(length(istart)==length(iend));

    for i = 1: length(istart)
        nonzero_block_ind = copula_ind(istart(i):iend(i));
        for j = 1:length(nonzero_block_ind)
            ind = nonzero_block_ind(j);
            [i2, i1] = ind2sub([n, m], ind);
            if i2+i1-k==0 % left
                if c1(i1, r) == 0
                    c1(i1, r) = 1;
                    c2(i2, r) = copula_select(ind);
                else
                    c2(i2, r) = copula_select(ind)/c1(i1, r);
                end
            else % right
                if c2(i2, r) == 0
                    c2(i2, r) = 1;
                    c1(i1, r) = copula_select(ind);
                else
                    c1(i1, r) = copula_select(ind)/c2(i2, r);
                end
            end
        end
    end

end
% if k <= n
%     c1(1) = 1;
%     c2(k) = copula_discrete(k, 1);
%     for i = 2:k
%         c1(i) = c1(i-1)*copula_discrete(k-i+1, i)/copula_discrete(k-i+1, i-1);
%         c2(k-i+1) = c2(k-i+2)*copula_discrete(k-i+1, i-1)/copula_discrete(k-i+2, i-1);
%     end
% elseif (k > n) && (k <= m)
%     c1(k-n) = 1;
%     c2(n) = copula_discrete(n, k-n);
%     for i = 1: n-1
%         c1(k-n+i) = c1(k-n+i-1)*copula_discrete(n-i+1, k-n+i)/copula_discrete(n-i+1, k-n+i-1);
%         c2(n-i) = c2(n-i+1)*copula_discrete(n-i, k-n+i)/copula_discrete(n-i+1, k-n+i);
%     end
%     i = n;
%     c1(k-n+i) = c1(k-n+i-1)*copula_discrete(n-i+1, k-n+i)/copula_discrete(n-i+1, k-n+i-1);
% else % k > m
%     c1(k-n) = 1;
%     c2(n) = copula_discrete(n, k-n);
%     for i = 1: m+n-k
%         c1(k-n+i) = c1(k-n+i-1)*copula_discrete(n-i+1, k-n+i)/copula_discrete(n-i+1, k-n+i-1);
%         c2(n-i) = c2(n-i+1)*copula_discrete(n-i, k-n+i)/copula_discrete(n-i+1, k-n+i);
%     end
% end

if flag_t
    tmp = c1;
    c1 = c2;
    c2 = tmp;
end

end