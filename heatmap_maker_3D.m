function heatmap_maker_3D(pixelsize, protein)

[H, stats] = heatmap_maker_kinet_abs(cd,[1200 1800],'remove',2,0);
for n = 1:((size(H,1)-1))/2
    new_row = H(n,:) + H(size(H,1)-(n-1),:);
    H(n,:) = new_row;
    H(size(H,1)-(n-1),:) = new_row;
end
%double mid-row
mid_row = H(ceil(size(H,1)/2),:);
H(ceil(size(H,1)/2),:) = mid_row*2;

H_int = sum(H,3);
pop_counts = size(stats.nm_2D,1);
figure;
h = surf(H_int);
h.FaceColor = 'interp';
h.EdgeColor = '[0,0,0]';
h.BackFaceLighting = 'unlit';
set(gca,'XTick',[320/pixelsize 640/pixelsize 960/pixelsize ...
    1280/pixelsize 1600/pixelsize 1920/pixelsize] );
set(gca,'XTickLabel',[320 640 960 1280 1600 1920] )
set(gca,'YTick',[(-640/pixelsize)+11 (-512/pixelsize)+11 (-384/pixelsize)+11 ...
    (-256/pixelsize)+11 (-128/pixelsize)+11 (0/pixelsize)+11 (128/pixelsize)+11 ...
    (256/pixelsize)+11 (384/pixelsize)+11 (512/pixelsize)+11 (640/pixelsize)+11] );
set(gca,'YTickLabel',[-640 -512 -384 -256 -128 ...
    0 128 256 384 512 640] )
set(gca,'ZTick',[(0/100)*pop_counts (1/100)*pop_counts (2/100)*pop_counts ...
    (3/100)*pop_counts (4/100)*pop_counts (5/100)*pop_counts (6/100)*pop_counts ...
    (7/100)*pop_counts (8/100)*pop_counts (9/100)*pop_counts (10/100)*pop_counts ...
    (11/100)*pop_counts (12/100)*pop_counts (13/100)*pop_counts (14/100)*pop_counts ...
    (15/100)*pop_counts (16/100)*pop_counts (17/100)*pop_counts (18/100)*pop_counts ...
    (19/100)*pop_counts (20/100)*pop_counts]);
set(gca,'ZTickLabel',0:20);
axis tight;
xlabel('Distance (nm)');
ylabel('Distance (nm)');
zlabel('Frequency (%)');
colormap jet;
mean_abs_y = mean(abs(stats.nm_2D(:,2)));
std_abs_y = std(abs(stats.nm_2D(:,2)),1);
mean_x = mean(abs(stats.nm_2D(:,1)));
std_x = std(abs(stats.nm_2D(:,1)),1);
plottitle = sprintf('%s, Tilt = 2, n = %d\nX = %d +/- %d nm\n AbsY = %d +/- %d nm',...
    protein,pop_counts,round(mean_x),round(std_x),...
    round(mean_abs_y), round(std_abs_y));
title(plottitle)
end