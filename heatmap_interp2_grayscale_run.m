function heatmap_interp2_grayscale_run(pixelsize,protein,x_trim,y_trim)
if mod(4,2) == 1
    error('y_trim must be odd integer\n');
end

for n = 2
[H, stats] = heatmap_maker_kinet(cd,[1200 1800],'remove',n,0);
H_integrated = sum(H,3);
H_integrated = H_integrated(1+y_trim:size(H_integrated,1)-y_trim,1:size(H_integrated,2)-x_trim,:);
mean_abs_y = mean(abs(stats.nm_2D(:,2)));
std_abs_y = std(abs(stats.nm_2D(:,2)),1);
mean_x = mean(abs(stats.nm_2D(:,1)));
std_x = std(abs(stats.nm_2D(:,1)),1);
counts = size(stats.nm_2D,1);
plottitle = sprintf('%s',protein);

conv = 64/pixelsize;

heatmap_interp2_grayscale(sum(H_integrated,3),pixelsize,conv);
end

end