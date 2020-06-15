%%IC-GN
%%sift初始点
function [p,window_center]=IC_GN2(I,I1,tri,i,pts1,pts2,window_size,tri1,tri2,iter,crition)
%%计算sift配准后点的初始位移，将其作为种子点
%I：参考图，I1：变形图
%pts1：sift计算的参考图特征点，pts2：变形图参考点
%window_size:窗口大小，奇数;
%tri:去误匹配后三角形点对，tri1：参考图，tri2：变形图
%iter：迭代次数，crition：dp停止准则，0.0001
%p：优化后的位移
%window_center:种子点位置，sift配准三角形中心
%参考文章：Comparative Analysis of Warp Function for Digital Image Correlation-Based Accurate Single-Shot 3D Shape Measurement


I=double(I);I1=double(I1);
x1=pts1(:,1)';
y1=pts1(:,2)';
x2=pts2(:,1)';
y2=pts2(:,2)';
index=tri(i);
seed_point1=[x1(tri1(index,:));y1(tri1(index,:))];
seed_point2=[x2(tri2(index,:));y2(tri2(index,:))];
seed_center=mean(seed_point1,2);
seed_point11=seed_point1-seed_center;
seed_point21=seed_point2-seed_center;
%%初始解
A1=[ones(3,1),seed_point11(1,:)',seed_point11(2,:)',zeros(3,3)];
A2=[zeros(3,3),ones(3,1),seed_point11(1,:)',seed_point11(2,:)'];
A=[A1;A2];
B=[(seed_point21(1,:)-seed_point11(1,:))';(seed_point21(2,:)-seed_point11(2,:))'];
p0=A\B;%初始解
p0=[p0(1) p0(2),p0(3),0,0,0,p0(4),p0(5),p0(6),0,0,0];
%%求hessian矩阵
%求梯度
[fx1,fy1]=gradient(I);
%防止后面点计算不了梯度，填充边界
pad=(window_size-1)/2;
[m,n]=size(I);
fx=zeros(m+window_size-1,n+window_size-1);
fx(pad+1:pad+m,pad+1:pad+n)=fx1;
fy=zeros(m+window_size-1,n+window_size-1);
fy(pad+1:pad+m,pad+1:pad+n)=fy1;
%window_size=27;
window_center=round(seed_center);
H=zeros(12,12);
F1=[];
for i=-(window_size-1)/2:(window_size-1)/2
    for j=-(window_size-1)/2:(window_size-1)/2
        W1=[1,i,j, 0.5* i*i,i*j,0.5*j*j ,zeros(1,6);zeros(1,6),1,i,j,0.5* i*i,i*j,0.5*j*j];
        F=[fx(j+window_center(2)+pad,i+window_center(1)+pad),fy(j+window_center(2)+pad,i+window_center(1)+pad)]*W1;
        H=H+F'*F;
        F1=[F1;F];
    end
end
%求p增量
p1=[];p=p0;
p1=[p1;p];
%iter=30;crition=0.000001;
for i=1:iter
warp_fun=[ 0.5*p(4),p(5),0.5*p(6),1+p(2),p(3),p(1);
    0.5*p(10),p(11),0.5*p(12),p(8),p(9)+1,p(7);
    zeros(1,5),1];
[win_x,win_y]=meshgrid(-(window_size-1)/2:(window_size-1)/2,-(window_size-1)/2:(window_size-1)/2);
ref_window=[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2)]+window_center;
cur_window=warp_fun*+[reshape(win_x,1,window_size^2).^2;reshape(win_x,1,window_size^2).*reshape(win_y,1,window_size^2);reshape(win_y,1,window_size^2).^2;reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2);ones(1,window_size^2)];
cur_window1(1,:)=cur_window(1,:)+window_center(1);
cur_window1(2,:)=cur_window(2,:)+window_center(2);
%dp=cur_window1-ref_window;
curgray=interp2(I1,cur_window1(1,:),cur_window1(2,:),'spline',0);
refgray=interp2(I,ref_window(1,:),ref_window(2,:),'spline',0);
%
refmean=mean(refgray);curmean=mean(curgray);
dref=sqrt(sum((refgray-refmean).^2));dcur=sqrt(sum((curgray-curmean).^2));
d=1;
SSD=(refgray-refmean-(dref/dcur)*(curgray-curmean));
SSD1=SSD.*F1';
dp=-inv(H)*sum(SSD1,2);
%更新p
Wdp=[ 0.5*dp(4),dp(5),0.5*dp(6),1+dp(2),dp(3),dp(1);
    0.5*dp(10),dp(11),0.5*dp(12),dp(8),dp(9)+1,dp(7);
    zeros(1,5),1];
%Wp=warp_fun*pinv(Wdp);
%p=p-dp';
Wp=t_warpfun(p)*inv(t_warpfun(dp));
p=[Wp(4,6),Wp(4,4)-1,Wp(4,5),2*Wp(4,1),Wp(4,2),2*Wp(4,3),Wp(5,6),Wp(5,4),Wp(5,5)-1,2*Wp(5,1),Wp(5,2),2*Wp(5,3)];
p1=[p1;p];
%CNSSD=sum(sum((refgray+F1*dp-refmean)/dref-(curgray-curmean)/dcur).^2);
ZNSSD=sum(((refgray-refmean)./dref-(curgray-curmean)./dcur).^2);
if (dp(1)^2+dp(7)^2)<crition
    break
end
end
end
function S=t_warpfun(p)
    %High-efficiency and high-accuracy digital image correlation for three-dimensional measurement
    %p=[u,ux,uy,uxx,uxy,uyy,v,vx,vy,vxx,vxy,vyy]
    %p=[1, 2, 3, 4,  5,   6,7, 8, 9, 10, 11,12]
    
    S(1,1)=1+2*p(2)+p(2)^2+p(1)*p(4);
    S(1,2)=2*p(1)*p(5)+2*(1+p(2))*p(3);
    S(1,3)=p(3)^2+p(1)*p(6);
    S(1,4)=2*p(1)*(1+p(2));
    S(1,5)=2*p(1)*p(3);
    S(1,6)=p(1)^2;

    S(2,1)=0.5*(p(7)*p(4)+2*(1+p(2))*p(8)+p(1)*p(10));
    S(2,2)=1+p(3)*p(8)+p(2)*p(9)+p(7)*p(5)+p(1)*p(11)+p(9)+p(2);
    S(2,3)=0.5*(p(7)*p(6)+2*p(3)*(1+p(9))+p(1)*p(12));
    S(2,4)=p(7)+p(7)*p(2)+p(1)*p(8);
    S(2,5)=p(1)+p(7)*p(3)+p(1)*p(9);
    S(2,6)=p(1)*p(7);
   %p=[u,ux,uy,uxx,uxy,uyy,v,vx,vy,vxx,vxy,vyy]
    %p=[1, 2, 3, 4,  5,  6,7, 8, 9, 10, 11,12] 
    S(3,1)=p(8)^2+p(7)*p(10);
    S(3,2)=2*p(7)*p(11)+2*p(8)*(1+p(9));
    S(3,3)=1+2*p(9)+p(9)^2+p(7)*p(12);
    S(3,4)=2*p(7)*p(8);
    S(3,5)=2*p(7)*(1+p(9));
    S(3,6)=p(7)^2;
    
    S(4:6,:)=[ 0.5*p(4),p(5),0.5*p(6),1+p(2),p(3),p(1);
        0.5*p(10),p(11),0.5*p(12),p(8),p(9)+1,p(7);
        zeros(1,5),1];
end