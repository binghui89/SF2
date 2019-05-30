function [z, z_pdf] = my_convolution(u, v, binwidth)
% Return the convolution of sample sets u and v, u and v are sets of
% random numbers, binwidth is the bin width to get the histograms. We
% assume within each bin u and v are uniformly distributed.
u_min = floor(min(u)./binwidth)*binwidth;
u_max = ceil(max(u)./binwidth)*binwidth;
[u_counts, u_edges] = histcounts(u, u_min:binwidth:u_max);
u_bincenter = (u_edges(1: end-1) + u_edges(2: end))/2;

v_min = floor(min(v)./binwidth)*binwidth;
v_max = ceil(max(v)./binwidth)*binwidth;
[v_counts, v_edges] = histcounts(v, v_min:binwidth:v_max);
v_bincenter = (v_edges(1: end-1) + v_edges(2: end))/2;

z = [u_edges(1)+v_edges u_edges(2:end)+v_edges(end)]; % z = u + v
z_pdf = [0 conv(u_counts./sum(u_counts), v_counts./sum(v_counts)) 0]; % f(z)

z = z(:);
z_pdf = z_pdf(:);

end