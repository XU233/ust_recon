% clearvars
% close all

path = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-10-11\';
filename = 'lessOffCenter_20MHz_remap';
% filename = 'Data 2022-12-01 19-55-30_avg-remap';

Data_all = load([path, filename, '.mat']).VOLTAGE_remap;

% Data_all = Data_all(:,:,1);

%%
path = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-12-01\';
load([path, 'calibrated_coordinates\calibrated_coords_2.mat']);
valid_id = 1:512;
% valid_id(fault_id) = [];
element_num = length(valid_id);

xk = tx_x(valid_id,1);
yk = tx_x(valid_id,2);

fs = 20;    % sample rate [MHz]
delay = (2^13-1)/80*20; %750nm

%%
DAQ_all = single(Data_all(valid_id,:));
% DAQ_all = single(DAQ_all - mean(DAQ_all,1));
DAQ_all = permute(DAQ_all,[2,1,3]);

DAQ_all = reshape(DAQ_all,size(DAQ_all,1),512,(numel(DAQ_all)/512/size(DAQ_all,1)));

%% 
x_range = 128;   % [mm]
y_range = 128;   % [mm]
res_factor = 2; % resolution factor: no. of pixels per mm
center_x = 0;   % [mm]
center_y = 0;   % [mm]

X = (-(x_range*res_factor/2-1):(x_range*res_factor/2))/res_factor - center_x; % FOV
Y = (-(y_range*res_factor/2-1):(y_range*res_factor/2))/res_factor - center_y;

% X = (-255:256)*0.2-0; % FOV
% Y = (-255:256)*0.2+0;

% v = linspace(1.486,1.502,6);
v = 1.483; % water speed of sound [km/s]

for i = 1:length(delay)
% for i = 1:length(v)

tic
% reIMG3 = subfunc_3d_cuda_reduce(DAQ_all,X,Y,Z,v(i),delay,xk,yk,zk,fs);
reIMG2 = subfunc_2d_cuda_reduce(DAQ_all,X,Y,v,delay(i),xk,yk,fs);
toc

figure('Name','XY'); 
XY_view = squeeze(max((reIMG2(:,:)),[],3));
imshow(XY_view,[]) % Top view, XY
xlabel ('--y++')
ylabel ('--x++')

% [mx(i),ind] = max(reIMG3(:));
end

%% display
dynamic_range = [-10,0];
imagedata = abs(reIMG2) ./ max(abs(reIMG2(:)));
imagelog = 20 * log10(imagedata);
figure();
% imagesc(Display.xm(x_ind) * 1e2,Display.ym(y_ind) * 1e2, imagelog(y_ind,x_ind), dynamic_range);colormap(gray); colorbar;
imagesc(imagelog, dynamic_range);colormap(gray); colorbar;
% xlabel('X (cm)');ylabel('Y (cm)');axis('image');
