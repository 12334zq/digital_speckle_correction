function [p,window_center]=RG_DIC_IG_GN(I,I1,window_size,window_center,iter,crition,p)
%%计算邻域的位移
%I：参考图，I1：变形图
%window_center:窗口中心，window_size:窗口大小
%iter：迭代次数，crition：停止准则：0.0001
%p窗口中心的母种子点
%mode=1：选用IC_GN,mode=0:选用IC_GN2
%参考文章：Reliability-guided digital image correlation for image deformation measurement
%High-efficiency and high-accuracy digital image correlation for three-dimensional measurement

I1=double(I1);I=double(I);
%i=round(rand(1)*length(tri));

%%IC_GN1;mode=1;else ICGN2

%求梯度
[fx1,fy1]=gradient(I);
%防止超出边界无法计算，对fx，fy填充。
pad=(window_size-1)/2;
[m,n]=size(I);
fx=zeros(m+window_size-1,n+window_size-1);
fx(pad+1:pad+m,pad+1:pad+n)=fx1;
fy=zeros(m+window_size-1,n+window_size-1);
fy(pad+1:pad+m,pad+1:pad+n)=fy1;

if length(p)==6
    
    %window_size=27;
    H=zeros(6,6);
    F1=[];
    for i=-(window_size-1)/2:(window_size-1)/2
        for j=-(window_size-1)/2:(window_size-1)/2
            W1=[1,i,j,0,0,0;0,0,0,1,i,j];
            F=[fx(j+window_center(2)+pad,i+window_center(1)+pad),fy(j+window_center(2)+pad,i+window_center(1)+pad)]*W1;
            H=H+F'*F;
            F1=[F1;F];
        end
    end
    %求p增量
    p1=[];
    p1=[p1;p];
    %iter=30;crition=0.000001;
    for i=1:iter
    warp_fun=[1+p(2),p(3),p(1);p(5),p(6)+1,p(4);0,0,1];
    [win_x,win_y]=meshgrid(-(window_size-1)/2:(window_size-1)/2,-(window_size-1)/2:(window_size-1)/2);
    ref_window=[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2)]+window_center;
    cur_window=warp_fun*+[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2);ones(1,window_size^2)];
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
    i
    end

elseif length(p)==12
    
    %window_size=27;
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
    p1=[];%p=p0;
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
    %求矩阵逆High-efficiency and high-accuracy digital image correlationfor three-dimensional measurement
    Wp=t_warpfun(p)*inv(t_warpfun(dp));
    %p=p-dp';
    p=[Wp(4,6),Wp(4,4)-1,Wp(4,5),2*Wp(4,1),Wp(4,2),2*Wp(4,3),Wp(5,6),Wp(5,4),Wp(5,5)-1,2*Wp(5,1),Wp(5,2),2*Wp(5,3)];
    p1=[p1;p];
    %CNSSD=sum(sum((refgray+F1*dp-refmean)/dref-(curgray-curmean)/dcur).^2);
    ZNSSD=sum(((refgray-refmean)./dref-(curgray-curmean)./dcur).^2);
    if (dp(1)^2+dp(7)^2)<crition
        break
    end
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