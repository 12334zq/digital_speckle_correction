%%
%close all
%特征点检测
function [tri,tri1,tri2,pts1,pts2]=delaunay_mesh(I,I1,edge_thresh)
%I:参考图
%I1:变形图
%edge_thresh:sift参数，2.5
%tri：匹配正确的三角形
%tri1：参考图的三角化序列
%tri2：变形图的三角化序列
%pts1：参考图特征点x1=pts1(:,1)';y1=pts1(:,2)';
%pts2:变形图特征点
%参考文章：Ef?cient Background Segmentation and Seed Point Generation for a Single-Shot Stereo System

addpath(genpath('matchSIFT'));
% image_i=imread('KyotoA.png');
% image_j=imread('KyotoB.png');
image_i=I;image_j=I1;
%sift特征检测
%edge_thresh = 3.5;
[SIFTloc_i,SIFTdes_i] = vl_sift(single((image_i)), 'edgethresh', edge_thresh) ;
[SIFTloc_j,SIFTdes_j] = vl_sift(single((image_j)), 'edgethresh', edge_thresh) ;
%sift match
[matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j);
pts1=SIFTloc_i(1:2,matchPointsID_i)';%x,y
pts2=SIFTloc_j(1:2,matchPointsID_j)';
scale1=SIFTloc_i(3,matchPointsID_i)';
scale2=SIFTloc_j(3,matchPointsID_j)';
alpha1=SIFTloc_i(4,matchPointsID_i)';
alpha2=SIFTloc_j(4,matchPointsID_j)';
% number = 20;
% x = 10*rand(1,number);
% y = 10*rand(1,number);
x1=pts1(:,1)';
y1=pts1(:,2)';
tri1 = delaunay(x1,y1);

x2=pts2(:,1)';
y2=pts2(:,2)';
tri2 = delaunay(x2,y2);


% image = [image_i image_j];
% figure,imshow(image)
% hold on;
% for i = 1 : length(pts1)
%     color = rand(1,3);
%     
%     plot([pts1(i, 1), size(image_i, 2) + pts2(i, 1)], [pts1(i, 2) pts2(i, 2)], 'Color', color)
%     
%     scatter(pts1(i, 1), pts1(i, 2), 40, color, 'filled');
%     scatter(size(image_i, 2) + pts2(i, 1), pts2(i, 2), 40, color,'filled');
% end
% figure
% subplot 121;imshow(image_i,[]);
% hold on
% plot(x1, y1, 'r*')
% 
% for ii = 1:size(tri1, 1)
%     plot( [x1(tri1(ii,1)) x1(tri1(ii,2))], [y1(tri1(ii,1)) y1(tri1(ii,2))], 'b' )
%     plot( [x1(tri1(ii,2)) x1(tri1(ii,3))], [y1(tri1(ii,2)) y1(tri1(ii,3))], 'b' )
%     plot( [x1(tri1(ii,1)) x1(tri1(ii,3))], [y1(tri1(ii,1)) y1(tri1(ii,3))], 'b' )
% end

% subplot 122;imshow(image_j,[]);
% hold on
% plot(x2, y2, 'r*')
% 
% for ii = 1:size(tri2, 1)
%   
%     plot( [x2(tri2(ii,1)) x2(tri2(ii,2))], [y2(tri2(ii,1)) y2(tri2(ii,2))], 'b' )
%     plot( [x2(tri2(ii,2)) x2(tri2(ii,3))], [y2(tri2(ii,2)) y2(tri2(ii,3))], 'b' )
%     plot( [x2(tri2(ii,1)) x2(tri2(ii,3))], [y2(tri2(ii,1)) y2(tri2(ii,3))], 'b' )
% end
% set(gca, 'box', 'on')
% print(gcf,'-dpng','delaunary.png')

%去除误匹配
tri=[];
for ii=1:size(tri2,1)
    if sum(tri1(ii))==sum(tri2(ii)) 
        [S1,angle1]=acquire_area(x1,y1,tri1(ii,:));
        [S2,angle2]=acquire_area(x2,y2,tri1(ii,:));
        if max(S1,S2)/min(S1,S2)<1.2 && angle1==1
            tri=[tri,ii];
        end
        
    end
end

% image = [image_i image_j];
% figure,imshow(image)
% hold on;
% for i =1:length(tri)
%     ii=tri(i);
%     plot( [x1(tri1(ii,1)) x1(tri1(ii,2))], [y1(tri1(ii,1)) y1(tri1(ii,2))], 'b' )
%     plot( [x1(tri1(ii,2)) x1(tri1(ii,3))], [y1(tri1(ii,2)) y1(tri1(ii,3))], 'b' )
%     plot( [x1(tri1(ii,1)) x1(tri1(ii,3))], [y1(tri1(ii,1)) y1(tri1(ii,3))], 'b' )
%     
%     plot( [(size(image_i, 2) +x2(tri2(ii,1))) ,(size(image_i, 2) +x2(tri2(ii,2)))], [y2(tri2(ii,1)) y2(tri2(ii,2))], 'b' )
%     plot( [(size(image_i, 2) +x2(tri2(ii,2))), (size(image_i, 2) +x2(tri2(ii,3)))], [y2(tri2(ii,2)) y2(tri2(ii,3))], 'b' )
%     plot( [(size(image_i, 2) +x2(tri2(ii,1))),( size(image_i, 2) +x2(tri2(ii,3)))], [y2(tri2(ii,1)) y2(tri2(ii,3))], 'b' )
%     
%     color = rand(1,3);
%     
%     plot([x1(tri1(ii,1)), size(image_i, 2) + x2(tri2(ii,1))], [y1(tri1(ii,1)) y2(tri2(ii,1))], 'Color', color)
%     
%     scatter(x1(tri1(ii,1)), y1(tri1(ii,1)), 40, color, 'filled');
%     scatter(size(image_i, 2) + x2(tri2(ii,1)), y2(tri2(ii,1)), 40, color,'filled');
% end



dx=x1(tri1(tri,1))-x2(tri2(tri,1));
dy=y1(tri1(tri,1))-y2(tri2(tri,1));
mean(dx);
mean(dy);
end
function [S,angle]=acquire_area(x,y,i)
angle=0;
ptx=x(i);
pty=y(i);
a=sqrt(((ptx(1)-ptx(2))^2)+((pty(1)-pty(2))^2));
b=sqrt(((ptx(1)-ptx(3))^2)+((pty(1)-pty(3))^2));
c=sqrt(((ptx(3)-ptx(2))^2)+((pty(3)-pty(2))^2));
S=1/4*sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a));
cosA=(c^2+b^2-a^2)/(2*b*c);A=acos(cosA)*180/pi;
cosB=(a^2+c^2-b^2)/(2*a*c);B=acos(cosB)*180/pi;
cosC=(a^2+b^2-c^2)/(2*a*b);C=acos(cosC)*180/pi;
if min([A,B,C])>20
    angle=1;
end
end

