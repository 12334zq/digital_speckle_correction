[I,I1]=testsample(400,400,1,2000,3.5,[14.9,0.004,0.050,-7,0,0.03],0.01,250);
%粒子群优化算法计算初始位移
%参考论文：Real-Time Digital Image Correlation for DynamicStrain Measurement
window_center=[200;200];window_size=21;
c1 = 2;
c2 = 2;
 
Gmax = 5;   %更新代数
sizepop = 25+10; %粒子数 根据情况确定

pxmax=15;pxmin=-25;%x,y方向最大位移，根据实际情况更改参数值
pymax=15;pymin=-15;
Vxmax = 0.5*(pxmax-pxmin);Vxmin=-Vxmax;
Vymax = 0.5*(pymax-pymin);Vymin=-Vymax;
[x,y]=meshgrid(-6:3:6,-6:3:6);
intial(1,:)=reshape(x,1,25);
intial(2,:)=reshape(y,1,25); 
p_intial=zeros(2,sizepop);
V=zeros(2,sizepop);
fitness=zeros(1,sizepop);
tic
%初始化粒子位置计算适应度
for i = 1:sizepop
    if i>25
        p_intial(1,i) =( (pxmax+pxmin)+(pxmax-pxmin)*rands(1,1))/2; 
        p_intial(2,i) =( (pymax+pymin)+(pymax-pymin)*rands(1,1))/2; 
    else
        p_intial(1,i)=intial(1,i);
        p_intial(2,i)=intial(2,i);
    end
    V(:,i) = rands(2,1);    
    fitness(i) = compute_zncc(I,I1,window_center,window_size,p_intial(:,i))  ; 
end
 %找到最优的粒子索引
[bestfitness, bestindex] = max(fitness);%代中所有粒子全局最优
gbest = p_intial(:,bestindex);   
pbest = p_intial;    %每个粒子的历史最优
fitnesspbest = fitness;   
fitnessgbest = bestfitness;   
Ctrust=0.75; 
for i = 1:Gmax
    %判断是否收敛
    if fitnessgbest>Ctrust
        break;
    end
    w=0.9-(i/(2*Gmax));
    for j = 1:sizepop   
        V(:,j) = w*V(:,j) + c1*rand*(pbest(:,j) - p_intial(:,j)) + c2*rand*(gbest - V(:,j));
        %V(:,j) = w*V(:,j) + c1*rand*(pbest(:,j) - p_intial(:,j)) + c2*rand*(gbest - p_intial(:,j));
        %判断粒子速度是否超出边界
        if V(1,j)>Vxmax
            V(1,j)=Vxmax;
        elseif V(1,j)<Vxmin
            V(1,j)=Vxmin;
        end
        
        if V(2,j)>Vymax
            V(2,j)=Vymax;
        elseif V(2,j)<Vymin
            V(2,j)=Vymin;
        end
        temp=p_intial(:,j);
        p_intial(:,j) = p_intial(:,j) + V(:,j);
        %判断粒子是否超出边界
       if p_intial(1,j)>pxmax
            p_intial(1,j)=pxmax;
        elseif p_intial(1,j)<pxmin
            p_intial(1,j)=pxmin;
        end
        
        if p_intial(2,j)>pymax
            p_intial(2,j)=pymax;
        elseif p_intial(2,j)<pymin
            p_intial(2,j)=pymin;
        end
        if temp(1)~=p_intial(1,j) || temp(2)~=p_intial(2,j)
            fitness(j) = compute_zncc(I,I1,window_center,window_size,p_intial(:,j)); 
            %更新粒子的历史最优解
            if fitness(j) > fitnesspbest(j)
                pbest(:,j) = p_intial(:,j);
                fitnesspbest(j) = fitness(j);
            end
            %更新全局最优解
            if fitnesspbest(j) > fitnessgbest
                gbest = pbest(:,j);
                fitnessgbest =fitnesspbest(j);
            end
        end
    end           
end
toc
gbest
fitnessgbest

%%BBGDS
Cd=0.9;
while fitnessgbest<Cd
    %八领域
    x=gbest(1);
    y=gbest(2);
    Candidate_p=[x-1,y-1;x-1,y;x-1,y+1;
        x,y-1;x,y;x,y+1;x+1,y-1;x+1,y;x+1,y+1];
    for i=1:9
        Candidate_zncc(i)=compute_zncc(I,I1,window_center,window_size,Candidate_p(i,:));
    end
    [fitnessgbest, bestindex] = max(Candidate_zncc);
    gbest=Candidate_p(bestindex,:);
    %如果中心点是最大的停止迭代
    if bestindex==5
        break;
    end
end
%IC-GN
p0=[gbest(1),0,0,gbest(2),0,0];
iter=30;crition=0.001;
[p,window_center,is_valid]=IC_GN(I,I1,window_size,window_center,iter,crition,p0)
function ZNCC=compute_zncc(I,I1,window_center,window_size,p_intial)
I1=double(I1);I=double(I);
[win_x,win_y]=meshgrid(-(window_size-1)/2:(window_size-1)/2,-(window_size-1)/2:(window_size-1)/2);
ref_window=[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2)]+window_center;
warp_fun=[1,0,p_intial(1);0,1,p_intial(2);0,0,1];
cur_window=warp_fun*[reshape(win_x,1,window_size^2);reshape(win_y,1,window_size^2);ones(1,window_size^2)];
cur_window1(1,:)=cur_window(1,:)+window_center(1);
cur_window1(2,:)=cur_window(2,:)+window_center(2);
curgray=interp2(I1,cur_window1(1,:),cur_window1(2,:),'spline',0);
refgray=interp2(I,ref_window(1,:),ref_window(2,:),'spline',0);
refmean=mean(refgray);curmean=mean(curgray);
refsum=sum(sum((refgray-refmean).^2));cursum=sum(sum((curgray-curmean).^2));
ZNCC=sum(sum((refgray-refmean).*(curgray-curmean)))./(sqrt(refsum)*sqrt(cursum));
end
