tic
window_size=27;iter=30;crition=0.0001;mode=1;
quick=1;d=10;
[m,n]=size(I);
seed_queue=[];deform=cell(m,n);
valid=zeros(m,n);%��Ǹ�λ�õ�p�Ƿ��Ż�
computed=zeros(m,n);%��Ǹ�λ���Ƿ��Ż����Ƿ��ڶ�����
queue=[];
grid_step=20;
[tri,tri1,tri2,pts1,pts2]=delaunay_mesh(I,I1,2.5);
for i=1:length(tri) 
    if mode==1
        [p,window_center]=IC_GN1(I,I1,tri,i,pts1,pts2,window_size,tri1,tri2,iter,crition);
    else 
        [p,window_center]=IC_GN2(I,I1,tri,i,pts1,pts2,window_size,tri1,tri2,iter,crition);
    end
    seed_queue=[seed_queue,window_center];
    deform{window_center(2),window_center(1)}=p;
    valid(window_center(2),window_center(1))=1;
    if quick
        %%�ӿ�����ٶȣ�����������Ѿ��ڶ��У�����dΪ���ڵĵ㲻�ټ���
        lu=max(window_center(2)-d,1);ld=min(window_center(2)+d,m);
        ru=max(window_center(1)-d,1);rd=min(window_center(2)+d,n); 
        computed(lu:ld,ru:rd)=1;
    else
        computed(window_center(2),window_center(1))=1;
    end
    %�����ĸ�����,���������������
    [queue,computed]=sel_seed(I,I1,window_center,p,window_size,grid_step,queue,computed);
%         if valid(window_center(2),window_center(1))==0
%             [p,window_center]=RG_DIC_IG_GN(I,I1,window_size,window_center,iter,crition,p);
%             deform{window_center(2),window_center(1)}=p;
%             valid(window_center(2),window_center(1))=1;
%         end
end

%%�����в�Ϊ�㣬����valid==0���������ӵ�����
while isempty(queue(3,:))==0
    [~,index]=find(queue(3,:)==min(queue(3,:)));
    if isempty(index)==1
        index=1;
    end
    select_center=queue(1:2,min(index));
    prior_root=queue(4:5,min(index));
    queue(:,index)=[];
    if valid(select_center(2),select_center(1))==0
        
        p=deform{prior_root(2),prior_root(1)};
        [p,select_center]=RG_DIC_IG_GN(I,I1,window_size,select_center,iter,crition,p);%������е����ĵ��λ��
        deform{select_center(2),select_center(1)}=p;
        valid(select_center(2),select_center(1))=1;
            if quick
            %%�ӿ�����ٶȣ�����������Ѿ��ڶ��У�����dΪ���ڵĵ㲻�ټ���
                lu=max(select_center(2)-d,1);ld=min(select_center(2)+d,m);
                ru=max(select_center(1)-d,1);rd=min(select_center(1)+d,n);
                computed(lu:ld,ru:rd)=1;
            else
                computed(select_center(2),select_center(1))=1;
            end
        %computed(select_center(2),select_center(1))=1;
        [queue,computed]=sel_seed(I,I1,select_center,p,window_size,grid_step,queue,computed);%���¶���
    end
end
toc
[indexx,indexy]=find(valid==1);
p=[];
for i=1:length(indexx)
    p=[p,deform(indexx(i),indexy(i))];
end
    