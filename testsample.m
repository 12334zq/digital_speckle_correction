function [noisy_img,wimg]=testsample(mode,amount,r,p,b,T)
%%产生测试用例
%%产生散斑
%
%mode=0;表示非仿射变形，mode=1,仿射变形
%amount：散斑数量2000
%r:散斑尺寸3
%p:变形参数
%b：变形参数0.001，T周期250
I=zeros(216,216);
[m,n]=size(I);
%amount=fix(m*n*0.6);
%amount=2000;
%sampledata=zeros(amount,2);
r=r^2;
sampledata=[0.9*m*rand(1,amount);0.9*n*rand(1,amount)]';
for i=1:m
    for j=1:n
        x=i-sampledata(:,1);
        y=j-sampledata(:,2);
        a=exp(-(x.^2+y.^2)/r);
        I(j,i)=1*sum(a);
    end
end
%%产生变形图
img=double(I);

[h,w]=size(img);
%contol_point=[-n/2,n/2,-n/2,n/2;-m/2,-m/2,m/2,m/2];
%warp_point=warp_fun*[contol_point;ones(4,1)];

% Get all points in destination to sample
[xg, yg] = meshgrid(-w/2+1:w/2, -h/2+1:h/2);
xy = [reshape(xg, prod(size(xg)), 1)'; reshape(yg, prod(size(yg)), 1)'];
xy = [xy; ones(1,size(xy,2))];
if mode==1
    % Transform into source
   % p=[2.3,0.0101,0,3.5,-0.013,-0.0034];%
    M=[1+p(2),p(3),p(1);p(5),p(6)+1,p(4);0,0,1];
    uv = M * xy;
    % Remove homogeneous
    uv = uv(1:2,:)';
else
    %b=0.001;T=250;
    uv(1,:)=xy(1,:)+b*T*sin(2*pi*xy(1,:)/T);
    uv(2,:)=xy(2,:)+b*T*sin(2*pi*xy(2,:)/T);
    uv=uv';
end

% Sample
xi = reshape(uv(:,1),h,w)+w/2;
yi = reshape(uv(:,2),h,w)+h/2;
wimg = interp2(img, xi, yi, 'spline');

% Check for NaN background pixels - replace them with a background of 0
idx = find(isnan(wimg));
if ~isempty(idx)
	wimg(idx) = 0;
end
I1=wimg;
tmplt_pixel_sigma=0.4/max(wimg(:));%8;6;12;
image_pixel_sigma=0.4/max(I(:));
% Add noise to template
	if tmplt_pixel_sigma > 0
		wimg= wimg + (randn(size(wimg)) * tmplt_pixel_sigma);
	end
	
	% Add noise to image
	if image_pixel_sigma > 0
		noisy_img = I + (randn(size(I)) * image_pixel_sigma);
	else
		noisy_img = I;
	end
	
%[p,select_center]=RG_DIC_IG_GN(I,I1,61,[w/2;h/2],iter,crition,ans)