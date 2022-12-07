
clearvars
close all

% load('/scratch/spxdavis/15mhz_3d/PSdata_10172022/center_20MHz_remap.mat');%00_00_10.mat')
path = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-10-11\';
filename = 'lessOffCenter_20MHz_remap';
% filename = 'Data 2022-12-01 19-55-30_avg-remap';

load([path, filename, '.mat']);
data3 = single(mean(VOLTAGE_remap,3)');
%data3 = data3(:,:,:)/size(data3,2)/size(data3,2);

load('calibration_TAT.mat')

Nk = 512;

shift1 = (2^13-1)/80*20;


angles = 0;





data3 = data3(:,:,1:1:end);
xk = single(xk(1:1:end));
yk = single(yk(1:1:end));
zk = single(zk(1:1:end));



%t = sqrt((xt-xp).^2 + (xt-xp).^2 + (xt-xp).^2)/vp;


X = (-64:63)*0.25+xp(1); % FOV
Y = (-64:63)*0.25+yp(1);
Z = 0;%(-64:63)*0.5+zp(2);
%Z = Z(149);
gpuDevice(1)

% for v = linspace(1.46,1.6,10)
% fun = @(x) shift_v_sharpness_3d(v,data3, X, Y,Z,x,xk,yk,zk);
% x = fminsearch(fun,[shift1*1.49/v])  
% reIMG3 = subfunc_3d_cuda_reduce_weights(data3,X*v/1.49,Y*v/1.49,Z*v/1.49,v,x,xk,yk,zk,40);
% figure; imshow(squeeze(max(-reIMG3,[],3)),[])
% end
modifier = v(1);
%modifier = 1.4822;
%modifier = 1.46*linspace(0.99,1.01,31);
reIMG3 = zeros(numel(X),numel(Y),numel(Z));
ii=1;
 for i = 1:length(modifier)     

tic
reIMG3 = subfunc_3d_cuda_reduce(data3(:,:,1:1:end),X,Y,Z,modifier(i),shift1,xk(:,1:1:end),yk(:,1:1:end),zk(:,1:1:end),20);
%toc
%reIMG4 = subfunc_3d_cuda_reduce_weights(data3(:,:,2:2:end),X,Y,Z,modifier(i),shift1,xk(:,2:2:end),yk(:,2:2:end),zk(:,2:2:end),40);



%figure; imshow(squeeze(max(-reIMG3,[],1)),[])
%figure; imshow(squeeze(max(-reIMG3,[],2)),[])
figure; imshow(squeeze(max(reIMG3,[],3)),[])
%figure; imshow(squeeze(max(-reIMG4,[],1)),[])
%figure; imshow(squeeze(max(-reIMG4,[],2)),[])
%figure; imshow(squeeze(max(-reIMG4,[],3)),[])
%variance(i,ii) = -gather(sum(reIMG3(:).^2));

[mx(i,ii),ind] = max(reIMG3(:));

 end



    %figure; plot(modifier*1.49,mx);
figure; plot(mx);
%figure; plot(variance);

% reIMG3 = -reIMG3.*(reIMG3<0);    
% reIMG3 = padarray(rescale(reIMG3,0,255),[100 100 0],0,'both');
% 


% for i = 0:359
% red4 = squeeze(3*max(imrotate3(reIMG3,i,[0 0 1],'cubic','crop'),[],2));
% rgb = uint8(red4);
% 
% %size(rgb)
% imwrite(rgb,strcat(num2str(i,'%03.f'),'.tif'))
% 
% end    
    

%figure; plot(modifier*1.49,vr);
