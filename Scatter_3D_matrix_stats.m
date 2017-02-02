function [data_cell] = Scatter_3D_matrix_stats
%% By Brandon Friedman
% spindleLengths = [];
% threeDSpindleLengths = [];
%

%% Parse all .mat files in the directory
mat_files = dir('*.mat');
data_cell = load(mat_files(1).name,'data_cell');
data_cell = data_cell.data_cell;
coords = data_cell;
step = coords{2,7};
pixelsize = coords{2,8};
plane_separation = 2;
spindle_limits = [1200 1800];
for a=2:length(mat_files)
    data_cell = load(mat_files(a).name,'data_cell');
    data_cell = data_cell.data_cell;
    coords = cat(1,coords,data_cell(2:end,:));
end

%% Remove any entries where the two spb entries are the same or tilted
same_spb_bin = cellfun(@eq,coords(2:end,5),coords(2:end,6),'Un',0);
same_spb_array = cellfun(@sum,same_spb_bin);
spindle_bin = same_spb_array == 4;
%remove entries in which spindles are more than some Z plane distance
%apart
sbp_sub = cell2mat(cellfun(@minus,coords(2:end,5),...
    coords(2:end,6),'Un',0));
z_sep = abs(sbp_sub(:,3));
z_sep_bin = z_sep > plane_separation;
remove_bin = spindle_bin | z_sep_bin;
coords = coords(~([0;remove_bin]),:);
%% Remove single lac foci or remove all others
same_lac_bin = cellfun(@eq,coords(2:end,1),coords(2:end,3),'Un',0);
same_lac_array = cellfun(@sum,same_lac_bin);
remove_lac_bin = same_lac_array == 4;
coords = coords(~([0;remove_lac_bin]),:);

x = size(coords,1);

for n = 2:x
    % Make r1 the origin by subtracting it from everything
    coords{n,6} = coords{n,6} - coords{n,5};
    coords{n,1} = coords{n,1} - coords{n,5};
    coords{n,3} = coords{n,3} - coords{n,5};
    coords{n,5} = coords{n,5} - coords{n,5};
    
    % Convert pixels to nm for X and Y
    coords{n,6}(1:2) = coords{n,6}(1:2)*pixelsize;
    coords{n,1}(1:2) = coords{n,1}(1:2)*pixelsize;
    coords{n,3}(1:2) = coords{n,3}(1:2)*pixelsize;
    coords{n,5}(1:2) = coords{n,5}(1:2)*pixelsize;
    
    % Convert steps to nm for Z
    coords{n,6}(3) = coords{n,6}(3)*step;
    coords{n,1}(3) = coords{n,1}(3)*step;
    coords{n,3}(3) = coords{n,3}(3)*step;
    coords{n,5}(3) = coords{n,5}(3)*step;
    
    % Remove 4th data point
    coords{n,6} = coords{n,6}(1:3);
    coords{n,1} = coords{n,1}(1:3);
    coords{n,3} = coords{n,3}(1:3);
    coords{n,5} = coords{n,5}(1:3);
    
    
    %  Rotate r2 around the z and x axes
    %  Rotate around the x axis
    alpha = atan(coords{n,6}(3)/coords{n,6}(2));
    
    rotX = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
    coords{n,6} = coords{n,6}*rotX;
    coords{n,1} = coords{n,1}*rotX;
    coords{n,3} = coords{n,3}*rotX;
    
    %  Rotate around the Z axis
    theta = atan(coords{n,6}(2)/coords{n,6}(1));
    
    rotZ = -[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
    coords{n,6} = -coords{n,6}*rotZ;
    coords{n,1} = -coords{n,1}*rotZ;
    coords{n,3} = -coords{n,3}*rotZ;
end

% Select for spindle length and remove values that lie outside the SPB Range
coords_mat6 = cell2mat(coords(2:end,6));
coords_mat3 = cell2mat(coords(2:end,3));
coords_mat1 = cell2mat(coords(2:end,1));
coords_lt = coords_mat6(:,1) < spindle_limits(1);
coords_gt = coords_mat6(:,1) > spindle_limits(2);
coords_neg1 = coords_mat3(:,1) > coords_mat6(:,1);
coords_neg2 = coords_mat1(:,1) < 0;
coords_int = ~(coords_lt | coords_gt | coords_neg1 | coords_neg2);
coords_int = logical([0;coords_int]);
coords = coords(coords_int,:);

data_cell = coords;
clearvars -except data_cell
end