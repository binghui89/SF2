function [a_convolved, t_convolved] = conv_poly(a1, t1, a2, t2, delta_t)
% Piecewise polynomial convolution
% USE: [a_convoluted, t_convoluted] = CONV_POLY(a1, t1, a2, t2, delta_t)
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
%
% Reference:
% Based on R. J. Polge and R. D. Hays, "Numerical Technique for the 
% Convolution of Piecewise Polynomial Functions," in IEEE Transactions on 
% Computers, vol. C-22, no. 11, pp. 970-975, Nov. 1973.

% delta_t1 = mean(diff(t1));
% delta_t2 = mean(diff(t2));
% if abs(delta_t1 - delta_t2)/delta_t1 >= 0.01
%     warning('X grid size and Y grid size differs greather than 1%!');
% end
% delta_t = delta_t1;

% Convert both functions from piecewise polynomial functions to the delta
% function forms
c1 = c_by_a_and_t(a1, t1);
c2 = c_by_a_and_t(a2, t2);
% Use integer offset since the unique function works better with integers
t1_intoffset = round((t1 - t1(1))./delta_t);
t2_intoffset = round((t2 - t2(1))./delta_t);

% Convert both functions from the delta function form to the triplet form
% Each row is [coefficient, offset, degree]
% d1 = to_triplet_form(c1, t1);
% d2 = to_triplet_form(c2, t2);
d1 = to_triplet_form(c1, t1_intoffset);
d2 = to_triplet_form(c2, t2_intoffset);
r1 = size(d1, 1);
r2 = size(d2, 1);

d_convolved = nan(r1*r2, 1);
for i = 1: 3
    tmp1 = repmat(d1(:, i), r2, 1);
    tmp2 = repmat(d2(:, i), 1, r1)'; tmp2 = tmp2(:);
    if i == 1
        d_convolved(:, i) = tmp1.*tmp2;
    else
        d_convolved(:, i) = tmp1 + tmp2;
    end
end

% Merge similar terms (the same offset and degree)
tc1 = t1(1) + t2(1); % The left end of the first interval of the convolved distribution
% int_offset = round((d_convoluted(:, 2) - tc1)/delta_t);
% int_degree = round(d_convoluted(:, 3));
int_offset_degree = d_convolved(:, 2:3);
int_offset_degree_unique = unique(int_offset_degree, 'row');
d_convolved_unique_intoffset = nan(size(int_offset_degree_unique, 1), 3);
for i = 1: size(d_convolved_unique_intoffset, 1)
    row_selected = (...
        d_convolved(:, 2) == int_offset_degree_unique(i, 2-1) ...
        & d_convolved(:, 3) == int_offset_degree_unique(i, 3-1) ...
        );
    d_convolved_unique_intoffset(i, 1) = sum(d_convolved(row_selected), 1);
    d_convolved_unique_intoffset(i, 2:3) = int_offset_degree_unique(i, 1:2);
end

% Remove terms with zero coefficents
% d_convoluted_reduced_intoffset = d_convolved_unique_intoffset;
% d_convoluted_reduced_intoffset(d_convoluted_reduced_intoffset(:, 1)==0, :)=[];

% Convert the resulting triplet form to the delta function forms
[c_convolved, t_convolved_intoffset] = from_triplet_form(d_convolved_unique_intoffset);

% Conver the resulting delta function forms into the piecewise polynomial
% function forms
t_convolved = t_convolved_intoffset.*delta_t + tc1;
a_convolved = a_by_c_and_t(c_convolved, t_convolved);
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
% d(d(:, 1)==0, :) =[];
end