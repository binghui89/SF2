%% The Irwin�Hall distribution when n = 1 to 5
bin_width = 1;
a0 = [1; 0;]; % Size: I by M, I is number of segments + 1 (last segment is always 0), M is the highest order + 1 (lowest order is always 0)
t0 = [0; 1;]; % Size: I by 1
a = a0;
t = t0;

fig = figure();
for i = 1: 7
    plot_conv_poly(fig, a, t, rand(1,3));
    [a_convoluted, t_convoluted] = conv_poly(a0, t0, a, t, bin_width);
    a = a_convoluted; t = t_convoluted;
    a(end, :) = 0;
end
title('Convolution');

fig = figure();
for i = 1: 7
    [a, t] = irwin_hall(i);
    a = [a; zeros(1, size(a, 2))];
    plot_conv_poly(fig, a, t, rand(1,3));
end
title('Analytical form');
%% The paper example
c1 = [1; -1;]; t1 = [-3;-2;];
c2 = [0 2; 0 -4; 0 2; 0 2; 0 -4; 0 2; 0 2; 0 -4; 0 2;]; 
t2 = [0;0.5;1;2;2.5;3;4;4.5;5];

a1 = a_by_c_and_t(c1, t1); a2 = a_by_c_and_t(c2, t2);
[a_convoluted, t_convoluted] = conv_poly(a1, t1, a2, t2);
% a_convoluted = a_by_c_and_t(c_convoluted, t_convoluted);

fig = figure(); hold on;
plot_conv_poly(fig, a1, t1, 'b'); plot_conv_poly(fig, a2, t2, 'b');
plot_conv_poly(fig, a_convoluted, t_convoluted, 'k');


%% "Convolution" with correlation using triangular distributions
copula_discrete = [0.3 1; 1 1.7;];
bin_width = 1;

a = [0 1; 2 -1; 0 0;];
t = [0; 1; 2;];

% No correlation
[a_convoluted, t_convoluted] = conv_poly(a, t, a, t, bin_width);
fig = figure();
plot_conv_poly(fig, a_convoluted, t_convoluted, 'r');

% With correlation
% deltat = unique(diff(t)); % Samples must be equally placed
% tc_1 = t(1)+t(1);
% tc_e = t(end)+t(end);
% K = round((tc_e - tc_1)/deltat);
% 
% t_convcorr = nan(K+1, 1);
% a_convcorr = nan(K+1, size(a, 2) + size(a, 2));
% t_convcorr(1) = tc_1;
% a_convcorr(end, :) = 0;
% 
% for k = 1: K
%     [cx, cy] = copula_decompose(copula_discrete, k);
%     ax = repmat([cx; 0], 1, size(a, 2)).*a;
%     ay = repmat([cy; 0], 1, size(a, 2)).*a;
%     [a_tmp, t_tmp] = conv_poly(ax, t, ay, t);
%     a_convcorr(k, :) = a_tmp(t_tmp==t_convcorr(k), :);
%     t_convcorr(k+1) = tc_1+k*deltat;
% end
[a_convcorr, t_convcorr] = conv_poly_corr(a, t, a, t, copula_discrete, bin_width);
% fig = figure();
plot_conv_poly(fig, a_convcorr, t_convcorr, 'k');
