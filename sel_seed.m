function [queue,computed]=sel_seed(I,I1,window_center,p,window_size,grid_step,queue,computed)
%grid_step=5;
prior_root=window_center;
I=double(I);I1=double(I1);
neibour1=[0,0,-grid_step,grid_step;grid_step,-grid_step,0,0]+window_center;
%%判断领域是否超出边界
[m,n]=size(I);
neibour=[];
for i=1:4
    if neibour1(1,i)<=n && neibour1(1,i)>=1
        if neibour1(2,i)<=m && neibour1(2,i)>=1
            neibour=[neibour,neibour1(:,i)];
        end
    end
end
newqueue=[];
%%计算四个邻域的ZNCC,ZNSSD
% if size(neibour,2)==0
%     queue=queue;
%     computed=computed;
% else
ZNSSD=zeros(1, size(neibour,2));
    for i=1:size(neibour,2)
        sel_pixel=neibour(:,i);
        if computed(sel_pixel(2),sel_pixel(1))==0    %判断该邻域点是否已经计算过
            warp_fun=[1+p(2),p(3),p(1);p(5),p(6)+1,p(4);0,0,1];
            [win_x,win_y]=meshgrid(-(window_size-1)/2:(window_size-1)/2,-(window_size-1)/2:(window_size-1)/2);
            ref_window=[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2)]+sel_pixel;
            cur_window=warp_fun*+[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2);ones(1,window_size^2)];
            cur_window1(1,:)=cur_window(1,:)+sel_pixel(1);
            cur_window1(2,:)=cur_window(2,:)+sel_pixel(2);
            curgray=interp2(I1,cur_window1(1,:),cur_window1(2,:),'spline',0);
            refgray=interp2(I,ref_window(1,:),ref_window(2,:),'spline',0);
            %
            refmean=mean(refgray);curmean=mean(curgray);
            dref=sqrt(sum((refgray-refmean).^2));dcur=sqrt(sum((curgray-curmean).^2));
            %ZNCC(i)=sum(sum((refgray-refmean).*(curgray-curmean)))./dref./dcur;
            ZNSSD(i)=sum(((refgray-refmean)./dref-(curgray-curmean)./dcur).^2);
            newqueue=[sel_pixel;ZNSSD(i);prior_root];
            queue=[queue,newqueue];
            computed(sel_pixel(2),sel_pixel(1))=1;
        end
    end

%     %%根据ZNSSD排序
    
%     %%弹出ZNSSD最小点
%     [~,index]=find(queue(3,:)==min(queue(3,:)));
%     window_center=newqueue(1:2,min(index));
end