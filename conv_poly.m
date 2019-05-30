function [a_convoluted, t_convoluted] = conv_poly(a1, t1, a2, t2)
% Based on R. J. Polge and R. D. Hays, "Numerical Technique for the 
% Convolution of Piecewise Polynomial Functions," in IEEE Transactions on 
% Computers, vol. C-22, no. 11, pp. 970-975, Nov. 1973.

% Convert both functions from piecewise polynomial functions to the delta
% function forms
c1 = c_by_a_and_t(a1, t1);
c2 = c_by_a_and_t(a2, t2);

% Convert both functions from the delta function form to the triplet form
d1 = to_triplet_form(c1, t1);
d2 = to_triplet_form(c2, t2);
r1 = size(d1, 1);
r2 = size(d2, 1);

d_convoluted = nan(r1*r2, 1);
for i = 1: 3
    tmp1 = repmat(d1(:, i), r2, 1);
    tmp2 = repmat(d2(:, i), 1, r1)'; tmp2 = tmp2(:);
    if i == 1
        d_convoluted(:, i) = tmp1.*tmp2;
    else
        d_convoluted(:, i) = tmp1 + tmp2;
    end
end

% Merge similar terms (the same offset and degree)
d_convoluted_unique = unique(d_convoluted(:, 2:3), 'row');
d_convoluted_unique = [nan(size(d_convoluted_unique, 1), 1), d_convoluted_unique];
for i = 1: size(d_convoluted_unique, 1)
    row_selected = (...
        d_convoluted(:, 2) == d_convoluted_unique(i, 2) ...
        & d_convoluted(:, 3) == d_convoluted_unique(i, 3) ...
        );
    d_convoluted_unique(i, 1) = sum(d_convoluted(row_selected), 1);
end

% Remove terms with zero coefficents
d_convoluted_reduced = d_convoluted_unique;
d_convoluted_reduced(d_convoluted_reduced(:, 1)==0, :)=[];

% Convert the resulting triplet form to the delta function forms
[c_convoluted, t_convoluted] = from_triplet_form(d_convoluted_reduced);

% Conver the resulting delta function forms into the piecewise polynomial
% function forms
a_convoluted = a_by_c_and_t(c_convoluted, t_convoluted);
end

function [c, t] = from_triplet_form(d)
M = max(-d(:, 3)); % The highest integral degree of delta function
t = sort(unique(d(:, 2)));
I = size(t, 1);
c = zeros(I, M);
for r = 1: size(d, 1)
    i = (t==d(r, 2));
    j = -d(r, 3);
    c(i, j) = d(r, 1);
end
end

function d = to_triplet_form(c, t)
% Each row of d is: coefficient, offset (t_i), degree of integral of the
% delta function
[row, col] = size(c);
t_matrix = repmat(t, 1, col);
degree_matrix = -repmat(1:col, row, 1);
d = [c(:) t_matrix(:) degree_matrix(:)];
d(d(:, 1)==0, :) =[];
end