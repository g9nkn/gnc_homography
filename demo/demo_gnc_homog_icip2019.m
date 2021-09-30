clear
close all
addpath(genpath('../'))

%% Load graf1 and graf3 pair and point correspondences.
% 'graf_1_3.mat' was obtained by using VLFeat and graffiti dataset:
%
% H_gt = load('graf/H1to3p');
% H_gt = H_gt / norm(H_gt,'fro');
% I1 = load('graf/img1.ppm');
% I2 = load('graf/img3.ppm');
% [frame1,desc1] = vl_sift(single(rgb2gray(I1)));
% [frame2,desc2] = vl_sift(single(rgb2gray(I2)));
% match = vl_ubcmatch(desc1,desc2);
% pt1 = frame1(1:2,match(1,:));
% pt2 = frame2(1:2,match(2,:));
%
% VLFeat  : https://www.vlfeat.org/
% graffiti: https://www.robots.ox.ac.uk/~vgg/data/affine/
load('graf_1_3.mat');

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');
nexttile; imshow(I1); title('graf1');
nexttile; imshow(I2); title('graf3');


%% Run the proposed GNC for homography estimation
% See all options by 
% >> help gnc_homog_nakano_icip2019 
[H_est, inliers_est] = gnc_homog_nakano_icip2019(pt1,pt2,'method','symtrans','weightfunc','truncL2','fastconv',true);

%% Numerical result (Frobenius error)
fprintf('H_gt = \n'); disp(H_gt)
fprintf('H_est = \n');disp(H_est)
froerr = min( norm(H_gt+H_est,'fro'), norm(H_gt-H_est,'fro') );
fprintf('Frobenius error of the homography estimation is %f\n', froerr);

%% Visualization
% plot matches
I3 = [I1, I2];
pt1_in = pt1(:,inliers_est);
pt2_in = pt2(:,inliers_est);
pt2_in(1,:) = pt2_in(1,:) + size(I1,2);

nexttile(3,[1,2]); 
imshow(I3); hold on
plot([pt1_in(1,:); pt2_in(1,:)], [pt1_in(2,:); pt2_in(2,:)], 'g-o', 'linewidth', 1, 'markerfacecolor', 'g');
title('Predicted inliers'); 


% image stitching
nexttile; imshow(stitch(I1, H_gt , I2, 1)); title('Image Stitching by H_gt' ,'interpreter','none');
nexttile; imshow(stitch(I1, H_est, I2, 1)); title('Image Stitching by H_est','interpreter','none');

scr_size = get(0,'ScreenSize');
fig_w = size(I3,2);
fig_h = min(size(I3,1)*3, scr_size(4));
set(gcf,'Position', [100,1,fig_w,fig_h])

