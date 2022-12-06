% clearvars
% close all

path = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-12-01\';
filename = 'Data 2022-12-01 19-55-30_avg-remap';

Data_all = load([path, filename, '.mat']).VOLTAGE_remap;


%%
load([path, '\calibrated_coordinates\calibrated_coords_2.mat']);
valid_id = 1:512;
valid_id(fault_id) = [];
element_num = length(valid_id);

xk = tx_x(valid_id,1);
yk = tx_x(valid_id,2);

fs = 5;    % sample rate [MHz]
delay = 0; %750nm

%%
DAQ_all = Data_all;
DAQ_all = single(DAQ_all - mean(DAQ_all,1));
DAQ_all = permute(DAQ_all,[2,1,3]);

DAQ_all = reshape(DAQ_all,size(DAQ_all,1),512,numel(DAQ_all)/512/size(DAQ_all,1));

%% 
x_range = 80;   % [mm]
y_range = 80;   % [mm]
z_range = 0;    % [mm]
res_factor = 3; % resolution factor: no. of pixels per mm
center_x = 0;   % [mm]
center_y = 0;   % [mm]
center_z = 0;   % [mm]

X = (-(x_range*res_factor/2-1):(x_range*res_factor/2))/res_factor - center_x; % FOV
Y = (-(y_range*res_factor/2-1):(y_range*res_factor/2))/res_factor - center_x;

% v = linspace(1.486,1.502,6);
v = 1.483; % water speed of sound [km/s]

for i = 1:length(delay)
% for i = 1:length(v)

tic
% reIMG3 = subfunc_3d_cuda_reduce(DAQ_all,X,Y,Z,v(i),delay,xk,yk,zk,fs);
reIMG3 = subfunc_2d_cuda_reduce(DAQ_all,X,Y,v,delay(i),xk,yk,fs);
toc

figure('Name','XY'); 
XY_view = squeeze(max((reIMG3(:,:,:)),[],3));
% XY_view = squeeze(max((reIMG3(:,:,size(X,2)/2:end)),[],3));
imshow(XY_view,[]) % Top view, XY
xlabel ('--y++')
ylabel ('--x++')

% [mx(i),ind] = max(reIMG3(:));
end

