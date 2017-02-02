function [activeinterpr, activeinterpg, activedatar, interpnumber] = heatmap_interp2_color_combine(RH,GH)

%% Mirror RH and GH matrix in Y direction
% code taken from Josh Lawrimore
% sum the top and bottom rows, moving toward the middle
for n = 1:((size(RH,1)-1))/2
    new_row = RH(n,:) + RH(size(RH,1)-(n-1),:);
    RH(n,:) = new_row;
    RH(size(RH,1)-(n-1),:) = new_row;
end
% double mid-row
mid_row = RH(ceil(size(RH,1)/2),:);
RH(ceil(size(RH,1)/2),:) = mid_row*2;

% sum the top and bottom rows, moving toward the middle
for n = 1:((size(GH,1)-1))/2
    new_row = GH(n,:) + GH(size(GH,1)-(n-1),:);
    GH(n,:) = new_row;
    GH(size(GH,1)-(n-1),:) = new_row;
end
% double mid-row
mid_row = GH(ceil(size(GH,1)/2),:);
GH(ceil(size(GH,1)/2),:) = mid_row*2;

%% Create the heatmap matricies
% code taken from Matthew Larson's Heatmap make program 
interpnumber = 2;

% creating the red channel matix
data_name_r = 'RH';
activedatar=eval(data_name_r);
activeinterpr=interp2(activedatar,interpnumber);%linear interpolate data
activeinterpr=activeinterpr/max(max(activeinterpr));%standardize to max=100%
eval([data_name_r 'interpr' num2str(interpnumber) '= activeinterpr;'])%create name for interpolated data
eval(['interpvarname =''' data_name_r 'interpr' num2str(interpnumber) ''';'])

% creating the green channel matrix
data_name_g = 'GH';
activedatag=eval(data_name_g);
activeinterpg=interp2(activedatag,interpnumber);%linear interpolate data
activeinterpg=activeinterpg/max(max(activeinterpg));%standardize to max=100%
eval([data_name_g 'interpg' num2str(interpnumber) '= activeinterpg;'])%create name for interpolated data
eval(['interpvarname =''' data_name_g 'interpg' num2str(interpnumber) ''';'])

% grayscale the two matricies
colormap gray
