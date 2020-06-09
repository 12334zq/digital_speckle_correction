I=zeros(512,512);
[m,n]=size(I);
%amount=fix(m*n*0.6);
amount=2000;
sampledata=zeros(amount,2);
r=3^2;
sampledata=[0.9*m*rand(1,amount);0.9*n*rand(1,amount)]';
% for i=1:amount
%     x=randi(m,1);
%     y=randi(n,1);
%     sampledata(i,:)=[x,y];
%     %I(x,y)=255;
% end
%imshow(I)
%[x,y]=meshgrid(1:m,1:n);
for i=1:m
    for j=1:n
        x=i-sampledata(:,1);
        y=j-sampledata(:,2);
        a=exp(-(x.^2+y.^2)/r);
        I(j,i)=1*sum(a);
    end
end
%I=I/max(I(:));
%imshow(I,[])
s=0;%0表示刚性平移形变，1一阶形变，2表示2阶形变
%%形变图
if s==0
    %sampledata1=sampledata+[0.2*m*(rand(1,amount)-0.5)',0.2*n*(rand(1,amount)-0.5)'];
    sampledata1=sampledata+[2.5,7.8];
elseif s==1
    ux=10*rand(1,amount)';
    uy=10*rand(1,amount)';
    sampledata1=sampledata+[ux,uy];
elseif s==2
end
for i=1:m
    for j=1:n
        x=i-sampledata1(:,1);
        y=j-sampledata1(:,2);
        a=exp(-(x.^2+y.^2)/r);
        I1(j,i)=1*sum(a);
    end
end
figure,imshow(I1,[])
hold on,plot(sampledata1(1,1),sampledata1(1,2),'+')
figure,imshow(I,[])
hold on,plot(sampledata(1,1),sampledata(1,2),'+')