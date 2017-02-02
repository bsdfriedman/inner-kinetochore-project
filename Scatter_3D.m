%% By Brandon Friedman and Derek Fulton
% spindleLengths = [];
% threeDSpindleLengths = [];
%

%% Parse all .mat files in the directory
mat_files = dir('*.mat');
for n=1:length(mat_files)
    data_cell = load(mat_files(n).name,'data_cell');
    data_cell = data_cell.data_cell;
    coords = data_cell;
    step = coords{2,7};
    pixelsize = coords{2,8};
    plane_separation = 2;
    spindle_limits = [1200 1800];

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

for i = 1:x-1;
    
    r1 = coords{i+1,5}(1:3);
    r2 = coords{i+1,6}(1:3);
    g1 = coords{i+1,1}(1:3);
    g2 = coords{i+1,3}(1:3);
    
    % Put z in nm
    r1(3) = (mod(r1(3)-1, 7) + 1)*step;
    r2(3) = (mod(r2(3)-1, 7) + 1)*step;
    g1(3) = (mod(g1(3)-1, 7) + 1)*step;
    g2(3) = (mod(g2(3)-1, 7) + 1)*step;
    
    r1(1) = r1(1)*pixelsize;
    r2(1) = r2(1)*pixelsize;
    g1(1) = g1(1)*pixelsize;
    g2(1) = g2(1)*pixelsize;
    
    r1(2) = r1(2)*pixelsize;
    r2(2) = r2(2)*pixelsize;
    g1(2) = g1(2)*pixelsize;
    g2(2) = g2(2)*pixelsize;
    
    % Make r1 the origin by subtracting it from everything
    r2 = r2 - r1;
    g1 = g1 - r1;
    g2 = g2 - r1;
    r1 = r1 - r1;
    
    %  Rotate r2 around the z and x axes
    %  Rotate around the x axis
    alpha = atan(r2(3)/r2(2));
    
    rotX = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
    r2 = r2*rotX;
    g1 = g1*rotX;
    g2 = g2*rotX;
    
    %  Rotate around the Z axis
    theta = atan(r2(2)/r2(1));
    
    rotZ = -[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
    r2 = -r2*rotZ;
    g1 = -g1*rotZ;
    g2 = -g2*rotZ;
    
    %  Normalize all spindle lengths to 1500 nm
    conv = 1500/r2(1);
    r2(1) = r2(1)*conv;
    g1(1) = g1(1)*conv;
    g2(1) = g2(1)*conv;
    
    %  Correct for rounding error for r2
    r2(2) = 0;
    r2(3) = 0;
    
    %  If r2 is greater than or less than spindle limits convert to nan
    lt_test = r2(1) < spindle_limits(1);
    gt_test = r2(1) > spindle_limits(2);
    neg_test1 = g2(1) > r2(1);
    neg_test2 = g1(1) < 0;
    if lt_test == 0 && gt_test == 0 && neg_test1 == 0 && neg_test2 == 0
              
        %  Plot the data
        xR = [r1(1); r2(1)];
        yR = [r1(2); r2(2)];
        zR = [r1(3); r2(3)];
        
        figure(1), scatter3(xR, yR, zR, 800, 'red.');
        
        hold on
        scatter3(g1(1), g1(2), g1(3), 200, 'green.');
        
        hold on
        scatter3(g2(1), g2(2), g2(3), 200, 'blue.');
    end
end
end
hold off;