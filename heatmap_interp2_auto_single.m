function heatmap_interp2_auto_single(H,pixelsize,plottitle,i)

%% Mirror H matrix in Y direction
%sum the top and bottom rows, moving toward the middle
for n = 1:((size(H,1)-1))/2
    new_row = H(n,:) + H(size(H,1)-(n-1),:);
    H(n,:) = new_row;
    H(size(H,1)-(n-1),:) = new_row;
end
%double mid-row
mid_row = H(ceil(size(H,1)/2),:);
H(ceil(size(H,1)/2),:) = mid_row*2;

%% Create the heatmap
% code taken from Matthew Larson's Heatmap make program 
data_name = 'H';
activedata=eval(data_name);
interpnumber = 2;
activeinterp=interp2(activedata,interpnumber);%linear interpolate data
activeinterp=activeinterp/max(max(activeinterp));%standardize to max=100%
% pixelsize=input('What is your initial bin size in nanometers?\n'); %input pixel size
% % plottitle=input('What is the title of you plot?\n','s'); %input title
eval([data_name 'interp' num2str(interpnumber) '= activeinterp;'])%create name for interpolated data
eval(['interpvarname =''' data_name 'interp' num2str(interpnumber) ''';'])
figure
%create new plot
imagesc(activeinterp)
title(plottitle)
xlabel('Distance (nm)')
ylabel('Distance (nm)')
%adjust axis labels
xlabels=round([0:2:(size(activedata,2)-1)]*pixelsize);
xlabels=(mat2cell(xlabels,1,ones(1,length(xlabels))));
ylabels=round(fliplr(([0:2:(size(activedata,1)-1)]-floor(size(activedata,1)/2))*pixelsize));
ylabels=(mat2cell(ylabels,1,ones(1,length(ylabels))));
set(gca,'XTick',1:(2*2^interpnumber):size(activeinterp,2))
set(gca,'XTickLabel',xlabels)
set(gca,'YTick',1:(2*2^interpnumber):size(activeinterp,1))
set(gca,'YTickLabel',ylabels)
colorbar
colormap jet
%colormap jet %remove % at the beginning of this line for jet (rainbow) colormap
axis image