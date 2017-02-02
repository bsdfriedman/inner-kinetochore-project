function [stats] = Scatter_3D_matrix_stats_run

% Load the coords file
mat_files = dir('*.mat');
data_cell = load(mat_files(1).name,'data_cell');
data_cell = data_cell.data_cell;
coords = data_cell;

% Spltit the coords file into smaller matricies
coords_mat1 = cell2mat(coords(1:end,1));
coords_mat3 = cell2mat(coords(1:end,3));
coords_mat6 = cell2mat(coords(1:end,6));

x_1 = coords_mat1(:,1);
x_2 = coords_mat6(:,1) - coords_mat3(:,1);
x_tot = cat(1,x_1,x_2);
x_tot = abs(x_tot);

y_1 = coords_mat1(:,2);
y_2 = coords_mat3(:,2);
y_tot = cat(1,y_1,y_2);
y_tot = abs(y_tot);

z_1 = coords_mat1(:,3);
z_2 = coords_mat3(:,3);
z_tot = cat(1,z_1,z_2);
z_tot = abs(z_tot);

for n = 1:length(y_tot);
    d{n,1} = sqrt(((y_tot(n,1))^2) + ((z_tot(n,1))^2));
end

d_n = cell2mat(d(1:end,1));

stats.data_points = size(x_tot,1);
stats.pixel_size = coords{2,8};
stats.step_size = coords{2,7};
stats.xmean = mean(x_tot,'omitnan');
stats.xstd = std(x_tot,'omitnan');
stats.ymean = mean(y_tot,'omitnan');
stats.ystd = std(y_tot,'omitnan');
stats.zmean = mean(z_tot,'omitnan');
stats.zstd = std(z_tot,'omitnan');
stats.axis_disp_mean = mean(d_n,'omitnan');
stats.axis_disp_std = std(d_n,'omitnan');