clc; clear all;
close all;

load absMax 

% 654 coils, 218 tetrahedron, 204 surfaces
% absACAm: ACA E-field amplitude (654*218), largest of the 360 orientations
% absADMm: ADM E-field amplitude (654*218)
% ROItri: ROI surface triangle points index (204*3)
% p: All the points location
% node5: index of the surface in the tetrahedron

% % Find the largest value
% maxv0=find(absACAm==max(max(absACAm)));
% idx=mod(maxv0,654)

idx=350; % 1~654
figure; 
subplot(1,2,1);
trisurf(ROItri,p(1,:)',p(2,:)',p(3,:)',abs(absADMm(idx,node5)),'edgealpha',0,'facealpha',.8);
title('ADM Result'); 
grid on; axis equal;
xlabel('x (m)');ylabel('y (m)'); zlabel('z (m)');
set(gca,'fontsize', 15);
colorbar;

subplot(1,2,2);
trisurf(ROItri,p(1,:)',p(2,:)',p(3,:)',absACAm(idx,node5),'edgealpha',0,'facealpha',.8);
title('ACA Result'); 
grid on; axis equal;
xlabel('x (m)');ylabel('y (m)'); zlabel('z (m)');
set(gca,'fontsize', 15);
colorbar;
set(gcf,'position',[10,300,1210,600])
