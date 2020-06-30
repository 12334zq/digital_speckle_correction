function [p,window_center,is_valid]=IC_GN(I,I1,window_size,window_center,iter,crition,p0)

%I：参考图，I1：变形图
%pts1：sift计算的参考图特征点，pts2：变形图参考点
%window_size:窗口大小，奇数;
%tri:去误匹配后三角形点对，tri1：参考图，tri2：变形图
%iter：迭代次数，crition：dp停止准则，0.0001
%p：优化后的位移
%window_center:种子点位置，可以为亚像素位置
%参考文章：Comparative Analysis of Warp Function for Digital Image Correlation-Based Accurate Single-Shot 3D Shape Measurement

I1=double(I1);I=double(I);

%%求hessian矩阵
%求梯度
[fx1,fy1]=gradient(I);
% %防止超出边界无法计算，对fx，fy填充。
% pad=(window_size-1)/2;
% [m,n]=size(I);
% fx=zeros(m+window_size-1,n+window_size-1);
% fx(pad+1:pad+m,pad+1:pad+n)=fx1;
% fy=zeros(m+window_size-1,n+window_size-1);
% fy(pad+1:pad+m,pad+1:pad+n)=fy1;
%window_size=27;
%梯度插值
[win_x,win_y]=meshgrid(-(window_size-1)/2:(window_size-1)/2,-(window_size-1)/2:(window_size-1)/2);
ref_window=[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2)]+window_center;
fx=interp2(fx1,ref_window(1,:),ref_window(2,:),'spline',0);fx=reshape(fx,window_size,window_size);
fy=interp2(fy1,ref_window(1,:),ref_window(2,:),'spline',0);fy=reshape(fy,window_size,window_size);
H=zeros(6,6);
F1=[];
for i=-(window_size-1)/2:(window_size-1)/2
    for j=-(window_size-1)/2:(window_size-1)/2
        W1=[1,i,j,0,0,0;0,0,0,1,i,j];
        %F=[bicubic_spline(fx,j+window_center(2)+pad,i+window_center(1)+pad),bicubic_spline(fy,j+window_center(2)+pad,i+window_center(1)+pad)]*W1;
        %F=[fx(j+window_center(2)+pad,i+window_center(1)+pad),fy(j+window_center(2)+pad,i+window_center(1)+pad)]*W1;
        F=[fx(j+(window_size-1)/2+1,i+(window_size-1)/2+1),fy(j+(window_size-1)/2+1,i+(window_size-1)/2+1)]*W1;
        H=H+F'*F;
        F1=[F1;F];
    end
end
%求p增量
p1=[];p=p0;
p1=[p1;p];
%iter=30;crition=0.000001;
for i=1:iter
warp_fun=[1+p(2),p(3),p(1);p(5),p(6)+1,p(4);0,0,1];
[win_x,win_y]=meshgrid(-(window_size-1)/2:(window_size-1)/2,-(window_size-1)/2:(window_size-1)/2);
ref_window=[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2)]+window_center;
cur_window=warp_fun*[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2);ones(1,window_size^2)];
cur_window1(1,:)=cur_window(1,:)+window_center(1);
cur_window1(2,:)=cur_window(2,:)+window_center(2);
%dp=cur_window1-ref_window;
curgray=interp2(I1,cur_window1(1,:),cur_window1(2,:),'spline',0);
refgray=interp2(I,ref_window(1,:),ref_window(2,:),'spline',0);
%
refmean=mean(refgray);curmean=mean(curgray);
dref=sqrt(sum((refgray-refmean).^2));dcur=sqrt(sum((curgray-curmean).^2));

SSD=(refgray-refmean-(dref/dcur)*(curgray-curmean));
SSD1=SSD.*F1';
dp=-inv(H)*sum(SSD1,2);
%更新p
Wdp=[1+dp(2),dp(3),dp(1);dp(5),dp(6)+1,dp(4);0,0,1];
Wp=warp_fun*inv(Wdp);
p=[Wp(1,3),Wp(1,1)-1,Wp(1,2),Wp(2,3),Wp(2,1),Wp(2,2)-1];
p1=[p1;p];
CNSSD=sum(sum((refgray+F1*dp-refmean)/dref-(curgray-curmean)/dcur).^2);
if (dp(1)^2+dp(4)^2)<crition
    break
end
end

if ssim((refgray-refmean)/dref,(curgray-curmean)/dcur)<0.4 || isnan(p(1)) || i>iter-2
    is_valid=0;
else
    is_valid=1;
    
end
end
function out=bicubic_spline(I,x,y)
x0=floor(x);
y0=floor(y);
if x0>2 && x0<size(I,1)-2 && y0>2 && y0<size(I,2)
    y1=[y0-1,y0,y0+1,y0+2];
    h=zeros(4,1);
    for i=1:4
        h(i)=[1,x-x0,(x-x0)^2,(x-x0)^3]*[0,1,0,0;-0.5,0,0.5,0;1,-2.5,2,-0.5;-0.5,1.5,-1.5,0.5]*[I(x0-1,y1(i));I(x0,y1(i));I(x0+1,y1(i));I(x0+2,y1(i))];
    end
    out=[1,y-y0,(y-y0)^2,(y-y0)^3]*[0,1,0,0;-0.5,0,0.5,0;1,-2.5,2,-0.5;-0.5,1.5,-1.5,0.5]*h;
else
    x0=round(x);y0=round(y);
    if x0>0 && x0<size(I,1) && y0>0 && y0<size(I,2)
        out=I(x0,y0);
    else
        out=0;
    end
end
end