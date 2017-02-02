%This program takes the X and Y coordinates from the Prox. Sub mean Coords
%and calculated a population Radius of Confinement using the variance from
%a Gaussian fit.  This program assumes that the data is 2D.  Program will
%search for the Excel_Summary.xls file in the current directory.

%%Read in the coordinates
function [Rc] = Pop_variance_rc(file,bead,chain)
load(file,'xstore','ystore');
X = squeeze(xstore(bead,chain,:));
Y = squeeze(ystore(bead,chain,:));
Xclear = [];
Yclear = [];
[rows, ~] = size(X);
for i = 1:rows
testnan = isnan(X(i));
if testnan == 0
    Xclear = [Xclear, X(i)];
end
end
for i = 1:rows
testnan = isnan(Y(i));
if testnan == 0
    Yclear = [Yclear, Y(i)];
end
end
x_mean = mean(Xclear);
y_mean = mean(Yclear);
[~, clearRows] = size(Xclear);
for i = 1:clearRows
x_sub_sq(i,1) = (Xclear(i) - x_mean)^2;
y_sub_sq(i,1) = (Yclear(i) - y_mean)^2;
x_sub(i,1) = (Xclear(i) - x_mean);
y_sub(i,1) = (Yclear(i) - y_mean);
end
[mu_x, sig_x] = normfit(x_sub);
[mu_y, sig_y] = normfit(y_sub);
sig_x_sq = sig_x^2;
sig_y_sq = sig_y^2;
sig_sq_total = (sig_x_sq + sig_y_sq)/2;
r_zero_sq = (mean(x_sub_sq) + (mean(y_sub_sq)));
Rc = (5/4)*(sqrt(2*sig_sq_total + r_zero_sq));

