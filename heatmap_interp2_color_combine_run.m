function heatmap_interp2_color_combine_run(pixelsize,proteinr,proteing,x_trim,y_trim)
% ensure that y-trim is an odd integer
if mod(4,2) == 1
    error('y_trim must be odd integer\n');
end

%% Making the GH and RH matrices for the heatmap (BH will be the sum of the two heatmaps
% choosing the two directories
for n = 2
rdir = sprintf('%s%sBioArk.bio.unc.edu%sBloomLab%sBrandon%sInner Kinetochore Project%sAll Matlab Files%smatlab_files_%s',filesep,filesep,filesep,filesep,filesep,filesep,filesep,proteinr);
gdir = sprintf('%s%sBioArk.bio.unc.edu%sBloomLab%sBrandon%sInner Kinetochore Project%sAll Matlab Files%smatlab_files_%s',filesep,filesep,filesep,filesep,filesep,filesep,filesep,proteing);

% making the matricies
[RH, Rstats] = heatmap_maker_kinet(rdir,[1200 1800],'remove',n,0);
RH_integrated = sum(RH,3);
RH_integrated = RH_integrated(1+y_trim:size(RH_integrated,1)-y_trim,1:size(RH_integrated,2)-x_trim,:);
[GH, Gstats] = heatmap_maker_kinet(gdir,[1200 1800],'remove',n,0);
GH_integrated = sum(GH,3);
GH_integrated = GH_integrated(1+y_trim:size(GH_integrated,1)-y_trim,1:size(GH_integrated,2)-x_trim,:);
%% Mean and std dev information, not used in this program
% mean_abs_y = mean(abs(stats.nm_2D(:,2)));
% std_abs_y = std(abs(stats.nm_2D(:,2)),1);
% mean_x = mean(abs(stats.nm_2D(:,1)));
% std_x = std(abs(stats.nm_2D(:,1)),1);
% counts = size(stats.nm_2D,1);

% titling the plot
plottitle = sprintf('Magenta Channel: %s vs Green Channel: %s',proteinr,proteing);

% conversion factor to scale all heatmaps to the same size
conv = 64/pixelsize;

% running the program to mirror and interpret the matricies
[activeinterpr, activeinterpg, activedatar, interpnumber] = heatmap_interp2_color_combine(sum(RH_integrated,3),sum(GH_integrated,3));

% assigning RGB channels to the matricies to get magenta and cyan heatmaps
RC = activeinterpr;
GC = activeinterpg;
BC = activeinterpr + activeinterpg;
a = zeros(size(activeinterpr, 1), size(activeinterpg, 2));

% If you want a magenta-cyan heatmap, use "cat(3,RC,GC,BC)"
% If you want a red-green heatmap, use "cat(3,RC,GC,a)"
% If you decide to do a red-green heatmap, remember to change the title of your plot
activeinterp = cat(3,RC,a,GC);

figure(6);
%create new plot
imagesc(activeinterp)
% title(plottitle)
xlabel('Distance (nm)')
ylabel('Distance (nm)')
%adjust axis labels
xlabels=round([0:2:(size(activedatar,2)-1)]*pixelsize*conv);
xlabels=(mat2cell(xlabels,1,ones(1,length(xlabels))));
ylabels=round(fliplr(([0:2:(size(activedatar,1)-1)]-floor(size(activedatar,1)/2))*pixelsize*conv));
ylabels=(mat2cell(ylabels,1,ones(1,length(ylabels))));
set(gca,'XTick',1:(2*2^interpnumber):size(activeinterp,2))
set(gca,'XTickLabel',xlabels)
set(gca,'YTick',1:(2*2^interpnumber):size(activeinterp,1))
set(gca,'YTickLabel',ylabels)
axis image
title(plottitle);
end

end