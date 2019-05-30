function h = plot_conv_poly(fig, a, t, c)
figure(fig); hold on;
for i = 1:(size(a, 1)-1)
    x_l = t(i);
    x_r = t(i+1);
    x = linspace(x_l, x_r, 10);
    y = polyval(fliplr(a(i, :)), x);
    h = plot(x, y, 'color', c);
end
end