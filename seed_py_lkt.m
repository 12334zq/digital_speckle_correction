%利用金字塔计算位移
%金字塔的层数与预估最大位移有关
%参考论文[2]Recursive-iterative digital image correlation based on salient features
%先求出特征点，这里是shi——Tomas，利用这些特征点为种子点，计算它在金字塔最高层（低分辨率，假设初始位移为0）的最优位移
%将高层位移*2做为低层的初始值
[I,I1]=testsample(400,400,1,2000,3.5,[14.8,0.000,0.000,-7,0,0.0],0.01,250);
noisy_per=std(I(:))/max(I(:));layer=4;
corners=detectMinEigenFeatures(I,'MinQuality',noisy_per);
%point=corners.selectStrongest(1000);
seed_point=corners.Location;
metric=corners.Metric;
figure,imshow(I,[]);
hold on ,plot(seed_point(:,1),seed_point(:,2),'+');
%%分块
[m,n]=size(I);
block=floor(seed_point/20);
xunique=unique(block(:,1));
yunique=unique(block(:,2));
first_point=[];
for i=1:size(xunique,1)
    [indexx,~]=find(block(:,1)==xunique(i));
    if isempty(indexx)
        continue
    else
        for j =1:size(yunique,1)
            [indexy,~]=find(block(indexx,2)==yunique(j));
            if isempty(indexy)
                continue
            else
                index=indexx(indexy);
                index2=find(metric(index)==max( metric(index)));
                seed_point(index,3)=seed_point(index(min(index2)),1);
                seed_point(index,4)=seed_point(index(min(index2)),2);  
                first_point=[first_point;[seed_point(index(min(index2)),1),seed_point(index(min(index2)),2)]];%需要先计算的种子点
                seed_point(index,5)=size(first_point,1);
            end
        end
    end
end
%%高斯金字塔
window_size=27;iter=30;crition=0.01;
ref=cell(4,1);ref{1}=I;
cur=cell(4,1);cur{1}=I1;
[m,n]=size(I);b=0.01;T=250;
[fx,fy]=gradient(I);
noise=max(std(fy(:)),std(fx(:)))*(window_size^2)*0.1;
for i=2:layer
    ref{i}=impyramid(ref{i-1},'reduce');
    cur{i}=impyramid(cur{i-1},'reduce');
end
p=zeros(size(first_point,1),2);
res=zeros(size(first_point,1),6);
window_size=27;iter=30;crition=0.0001;
flag=zeros(size(first_point,1),1);
for j=1:size(first_point,1)
    is_valid=1;
    for i=layer:-1:1
    refimage=ref{i};
    curimage=cur{i};
    window_center=first_point(j,:)/(2^(i-1)); 
        if i>1 
        %intial_p=2*[p(j,1),0,0,p(j,4),0,0];
            intial_p=2*[p(j,1),p(j,2)];
            [p(j,:),~,is_valid]=IC_GN_t(refimage,curimage,window_size,window_center',iter,noise,crition,intial_p); 
        else
            intial_p=2*[p(j,1),0,0,p(j,2),0,0];
            [res(j,:),~,is_valid]=IC_GN(refimage,curimage,window_size,window_center',iter,crition,intial_p); 
        end
        flag(j)=is_valid;
    end  
end
 res1=zeros(size(seed_point,1),6);
for i=1:size(seed_point,1)
    window_center=seed_point(i,1:2);
    intial_p=res(seed_point(i,5),:);
    intial_p=[intial_p(1),0,0,intial_p(4),0,0];
    if flag(seed_point(i,5))==1
        [res1(i,:),~,is_valid]=IC_GN(I,I1,window_size,window_center',30,crition,intial_p);
        valid(i)=is_valid;
    else
        continue
    end
end
        
[~,index]=find(valid==1);
error=[res1(index,1),res1(index,4)]+[14.8,-7];
error=error(:,1).^2+error(:,2).^2;
[index1,~]=find(error<2);
per=size(index1,1)/size(index,2);
num=size(index,2)/m/n
imshow(I,[])
hold on
plot(seed_point(index,1),seed_point(index,2),'+')
function [p,window_center,is_valid]=IC_GN_t(I,I1,window_size,window_center,iter,crition1,crition2,p0)

I1=double(I1);I=double(I);

%%求hessian矩阵
%求梯度
[fx1,fy1]=gradient(I);

%梯度插值
[win_x,win_y]=meshgrid(-(window_size-1)/2:(window_size-1)/2,-(window_size-1)/2:(window_size-1)/2);
ref_window=[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2)]+window_center;
fx=interp2(fx1,ref_window(1,:),ref_window(2,:),'spline',0);fx=reshape(fx,window_size,window_size);
fy=interp2(fy1,ref_window(1,:),ref_window(2,:),'spline',0);fy=reshape(fy,window_size,window_size);
H=zeros(2,2);
F1=[];
for i=-(window_size-1)/2:(window_size-1)/2
    for j=-(window_size-1)/2:(window_size-1)/2
        F=[fx(j+(window_size-1)/2+1,i+(window_size-1)/2+1),fy(j+(window_size-1)/2+1,i+(window_size-1)/2+1)];
        H=H+F'*F;
        F1=[F1;F];
    end
end
%求矩阵特征值
[~,lam]=eig(H);
if min(diag(lam))<crition1
    is_valid=0;
    p=p0;
else
    %求p增量
    p=p0;
    %iter=30;crition=0.000001;
    for i=1:iter
    warp_fun=[1,0,p(1);0,1,p(2);0,0,1];
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
    %Wdp=[1,0,dp(1);0,1,dp(2);0,0,1];
    %Wp=warp_fun*inv(Wdp);
    p=p-dp';
   
    CNSSD=sum(sum((refgray+F1*dp-refmean)/dref-(curgray-curmean)/dcur).^2);
    if (dp(1)^2+dp(2)^2)<crition2
        break
    end
    end

    if ssim((refgray-refmean)/dref,(curgray-curmean)/dcur)<0.4 || isnan(p(1))
        is_valid=0;
    else
        is_valid=1;

    end
end
end

