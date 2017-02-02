function [H, stats] = heatmap_maker_kinet_abs(directory,...
    spindle_limits,single_foci,plane_separation,pixel_round)
%This funciton will parse the .mat file from heatmap_GUI.m, rotate the
%spindle coordinates, and create a matrix to generate a heatmap using
%Matthew Larson's code.

%Options for single_foci are 'keep', 'remove', or 'only'.

%% Parse all .mat files in the directory
cd(directory);
mat_files = dir('*.mat');
for n=1:length(mat_files)
    data_cell = load(mat_files(n).name,'data_cell');
    data_cell = data_cell.data_cell;
    %parse pixel size
    pixel_size = data_cell{2,8};
    %% Remove any entries where the two spb entries are the same or tilted
    same_spb_bin = cellfun(@eq,data_cell(2:end,5),data_cell(2:end,6),'Un',0);
    same_spb_array = cellfun(@sum,same_spb_bin);
    spindle_bin = same_spb_array == 4;
    %remove entries in which spindles are more than some Z plane distance
    %apart
    sbp_sub = cell2mat(cellfun(@minus,data_cell(2:end,5),...
        data_cell(2:end,6),'Un',0));
    z_sep = abs(sbp_sub(:,3));
    z_sep_bin = z_sep > plane_separation;
    remove_bin = spindle_bin | z_sep_bin;
    data_cell = data_cell(~([0;remove_bin]),:);
    %% Remove single lac foci or remove all others
    same_lac_bin = cellfun(@eq,data_cell(2:end,1),data_cell(2:end,3),'Un',0);
    same_lac_array = cellfun(@sum,same_lac_bin);
    remove_lac_bin = same_lac_array == 4;
    switch single_foci
        case 'remove'
            display('Removing all single lacO foci')
            data_cell = data_cell(~([0;remove_lac_bin]),:);
        case 'only'
            display('Only using single lacO foci');
            data_cell = data_cell(logical([1;remove_lac_bin]),:);
        case 'keep'
            display('Keeping all data');
            data_cell = data_cell(:,:);
    end
    %% Register and Rotate Coords
    %Register all coords by SPB1
    reg_lac1 = cellfun(@minus, data_cell(2:end,1),data_cell(2:end,5),'Un',0);
    reg_lac2 = cellfun(@minus, data_cell(2:end,3),data_cell(2:end,5),'Un',0);
    reg_SPB2 = cellfun(@minus, data_cell(2:end,6),data_cell(2:end,5),'Un',0);
    %% Rotate coords
    %only start this portion of program if data_cell has entries
    if size(data_cell,1) < 2
        %create x and y axes
        xDim = (-2:30)';
        yDim = (-10:10)';
        ff.H(:,:,n) = zeros(length(yDim),length(xDim));
        fs.H(:,:,n) = zeros(length(yDim),length(xDim));
        sf.H(:,:,n) = zeros(length(yDim),length(xDim));
        ss.H(:,:,n) = zeros(length(yDim),length(xDim));
    else
        %Parse just the X and Y Coords of SPB2 for theta
        spb_XY = cell2mat(cellfun(@(x) x(1:2),reg_SPB2,'Un',0));
        lac1_XY = cell2mat(cellfun(@(x) x(1:2),reg_lac1,'Un',0));
        lac2_XY = cell2mat(cellfun(@(x) x(1:2),reg_lac2,'Un',0));
        thetas = atan(spb_XY(:,2)./spb_XY(:,1));
        for i = 1:length(thetas)
            rot_mat(:,:,i) = [cos(thetas(i)), -sin(thetas(i)); sin(thetas(i)), cos(thetas(i))];
            rot_spb_XY(i,:) = spb_XY(i,:)*squeeze(rot_mat(:,:,i));
            rot_lac1_XY(i,:) = lac1_XY(i,:)*squeeze(rot_mat(:,:,i));
            rot_lac2_XY(i,:) = lac2_XY(i,:)*squeeze(rot_mat(:,:,i));
            norm1(i,1) = norm(rot_lac1_XY(i,:));
            norm1(i,2) = norm(lac1_XY(i,:));
            norm2(i,1) = norm(rot_lac2_XY(i,:));
            norm2(i,2) = norm(lac2_XY(i,:));
            normspb(i,1) = norm(rot_spb_XY(i,:));
            normspb(i,2) = norm(spb_XY(i,:));
            normlac(i,1) = norm(lac1_XY(i,:) - lac2_XY(i,:));
            normlac(i,2) = norm(rot_lac1_XY(i,:) - rot_lac2_XY(i,:));
        end
        if sum(round(norm1(:,1),4)==round(norm1(:,2),4)) ~= size(norm1,1)
            error('Error in lac1 rotation. Distance not equal\n');
        end
        if sum(round(norm2(:,1),4)==round(norm2(:,2),4)) ~= size(norm2,1)
            error('Error in lac2 rotation. Distance not equal\n');
        end
        if sum(round(normspb(:,1),4)==round(normspb(:,2),4)) ~= size(normspb,1)
            error('Error in SPB rotation. Distance not equal\n');
        end
        if sum(round(normlac(:,1),4)==round(normlac(:,2),4)) ~= size(normlac,1)
            error('Error in lac rotation. Distance not equal\n');
        end
        %% Round to nearest pixel if pixel_round equals 1
        if pixel_round == 1
            rot_spb_XY = round(rot_spb_XY,0);
            rot_lac1_XY = round(rot_lac1_XY,0);
            rot_lac2_XY = round(rot_lac2_XY,0);
        end
        %% Calculate absolute X and Y distance from origin and pole
        %Check the sign of each of the points    
        lac1_pos = rot_lac1_XY(:,1) >= 0;
        lac2_pos = rot_lac2_XY(:,1) >= 0;
        spb2_pos = rot_spb_XY(:,1) >= 0;
        sign_check = lac1_pos + lac2_pos + spb2_pos;
        %flip the sign of all the 0s
        bin_0 = sign_check == 0;
        mult_0 = (-2*bin_0) + 1;
        rot_lac1_XY(:,1) = rot_lac1_XY(:,1) .* mult_0;
        rot_lac2_XY(:,1) = rot_lac2_XY(:,1) .* mult_0;
        rot_spb_XY(:,1) = rot_spb_XY(:,1) .* mult_0;
        %run sign check again
        lac1_pos = rot_lac1_XY(:,1) >= 0;
        lac2_pos = rot_lac2_XY(:,1) >= 0;
        spb2_pos = rot_spb_XY(:,1) >= 0;
        sign_check = lac1_pos + lac2_pos + spb2_pos;
        %keep only the 3s
        bin_3 = sign_check == 3;
        rot_lac1_XY = rot_lac1_XY(bin_3,:);
        rot_lac2_XY = rot_lac2_XY(bin_3,:);
        rot_spb_XY = rot_spb_XY(bin_3,:);
        
        %% Check that the x values are in order of lac1,lac2,spb2
        x_order_bin = rot_lac1_XY(:,1) <= rot_lac2_XY(:,1) & ...
            rot_lac2_XY(:,1) <= rot_spb_XY(:,1);
        %Keep only cells in which they are in order L1,L2,S2 in X
        rot_lac1_XY = rot_lac1_XY(x_order_bin,:);
        rot_lac2_XY = rot_lac2_XY(x_order_bin,:);
        rot_spb_XY = rot_spb_XY(x_order_bin,:);
        
        %% Caclulate the xdist1 and xdist2
        xdist1_final = rot_lac1_XY(:,1);
        xdist2_final = rot_spb_XY(:,1) - rot_lac2_XY(:,1);
        ydist1_final = rot_lac1_XY(:,2);
        ydist2_final = rot_lac2_XY(:,2);
        %% Create matrix of absolute,registered, and rotated coords and spindle length
        %convert the rot_spb2 xcoords to nm for spindle length
        spindle_nm = rot_spb_XY(:,1) * pixel_size;
        mat_2D = [xdist1_final,ydist1_final,...
            xdist2_final, ydist2_final,spindle_nm];
        %% Apply spindle_limit
        spindle_bin = mat_2D(:,end) >= spindle_limits(1) & ...
            mat_2D(:,end) <= spindle_limits(2);
        mat_2D = mat_2D(spindle_bin,:);
        return_mat{n} = mat_2D;
        H_mat = [mat_2D(:,1), mat_2D(:,2);...
            mat_2D(:,3), mat_2D(:,4)];
        %% create H matrices for each category
        %create x and y axes
        xDim = (0:30)';
        yDim = (-10:10)';
        %% foci with foci
        H(:,:,n) = zeros(length(yDim),length(xDim));
        for l = 1:size(mat_2D,1)
            countX = dsearchn(xDim, H_mat(l,1));
            countY = dsearchn(yDim, H_mat(l,2));
            H(countY,countX,n) = H(countY,countX,n) + 1 ;
        end
    end
    %% clear all variables except H's mat files and n
    clearvars -except H  return_mat n mat_files directory ...
        spindle_limits single_foci pixel_size plane_separation pixel_round
end

%% calcuate x and y means
%setup the stats structure
stats.mat_2D = return_mat;
stats.nm_2D = [];
stats.pixel_size = pixel_size;
for k = 1:length(return_mat)
    if isempty(return_mat{1,k}) == 0
        stats.nm_2D = [stats.nm_2D;[return_mat{1,k}(:,1:2)*pixel_size,return_mat{1,k}(:,5);...
            return_mat{1,k}(:,3:4)*pixel_size,return_mat{1,k}(:,5)]];
    end
end
%calc means, diameter and stds
if isempty(stats.nm_2D) == 0
    stats.xmean = mean(stats.nm_2D(:,1));
    stats.xstd = std(stats.nm_2D(:,1),1);
    stats.ymean = mean(abs(stats.nm_2D(:,2)));
    stats.ystd = std(abs(stats.nm_2D(:,2)),1);
end
