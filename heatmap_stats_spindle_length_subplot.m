function heatmap_stats_spindle_length_subplot(pixelsize,protein,x_trim,y_trim)
if mod(4,2) == 1
    error('y_trim must be odd integer\n');
end

for n = 0:5
[H, stats] = heatmap_maker_kinet(cd,[(1200+(100*n)) (1300+(100*n))],'remove',2,0);
H_integrated = sum(H,3);
H_integrated = H_integrated(1+y_trim:size(H_integrated,1)-y_trim,1:size(H_integrated,2)-x_trim,:);
mean_abs_y = mean(abs(stats.nm_2D(:,2)));
std_abs_y = std(abs(stats.nm_2D(:,2)),1);
mean_x = mean(abs(stats.nm_2D(:,1)));
std_x = std(abs(stats.nm_2D(:,1)),1);
counts = size(stats.nm_2D,1);
spmin = 1200+(100*n);
spmax = 1300+(100*n);
plottitle = sprintf('%s, Tilt = 2, n = %d\nSpindle Length = %d to %d\nX = %d +/- %d nm\n AbsY = %d +/- %d nm',...
    protein,counts,spmin,spmax,round(mean_x),round(std_x),...
    round(mean_abs_y), round(std_abs_y));

heatmap_interp2_auto_subplot(sum(H_integrated,3),pixelsize,plottitle,n+1);
end

end