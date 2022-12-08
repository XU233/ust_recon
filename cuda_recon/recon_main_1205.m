clear all; clc

%% load data
Para.datapath = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-10-11\';
Para.dataname = 'lessOffCenter_20MHz_remap'; %

Para.display = 0; % Display raw signal or not
Para.hilbert = 0; % single to complex single

[Data_all, Acq]= f_load_data(Para);

% Data_all = raw_data(:,:,1:10);
% Data_all = sum(raw_data,3);

DAQ_all = single(Data_all(:,:));
DAQ_all = single(DAQ_all - mean(DAQ_all,1));

DAQ_all = reshape( DAQ_all, 512, size(DAQ_all,2), numel(DAQ_all)/512/size(DAQ_all,2) );

% %% transducer setting
% Trans.t_foc_bias = 72.7e-3; %72 * 1e-3; % distance between surface and arc, % for air-air
% 
% Trans.t_init_angle = 90; %180; % degrees, % -87
% 
% Trans.t_Nsteps = 403; %433; %595; %Acq.Nf;% no. of positions = no. of frames
% Trans.t_foclens = 1.00 * 25.4 * 1e-3;% focus length
% Trans.t_clockwise = 1; % rotation direction, Tbc.
% Trans.r_Aapo = 180; %  (degrees) Angle range of detection array elements for recon
% Trans.x_offset = 0e-3; %-5e-3; %0; %12e-3;
% Trans.y_offset = 0e-3; %-5e-3; %17e-3;
% Trans.t_scan_angle = 2*pi; % 2*pi
% Trans = f_trans_array(Trans); % Trans.rapo_list;
% 
% % xk = tx_x(valid_id,1)*1e3;
% % yk = tx_x(valid_id,2)*1e3;
% 
% % DAQ_all = DAQ_all(Trans.valid_id,:,:);
% 
% %% display field of view
% Display.x_range = 80 * 1e-3; % m
% Display.y_range = 80 * 1e-3; % m
% Display.z_range = 0 * 1e-3; % m
% Display.res_factor = 3; % resolution factor: no. of pixels per mm
% Display.center_x = 0*1e-3;
% Display.center_y = 0*1e-3;
% Display.center_z = 0 * 1e-3; % m
% 
% % acquisition para
% Acq.c = 1.483e3; %1.495 * 1e3;% speed of sound (m/s)
% Acq.delay = 0; %-(Trans.t_foclens)/Acq.c; % trigger delay of DAQ between US transmission
% Acq.delay_revise_m =(0)*1e-6 * Acq.c + Acq.delay; % DAQ delay , 25 sampling error of DAQ = 0.625us

%%
path = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-12-01\';
load([path, 'calibrated_coordinates\calibrated_coords_2.mat']);

load('calibration_TAT.mat');

xk = xk * 1e-3;
yk = yk * 1e-3;

fs = 20e6;    % sample rate
delay = (2^13-1)/80*20; % index

%% field of view
x_range = 128 * 1e-3;
y_range = 128 * 1e-3;
res_factor = 2; % resolution factor: no. of pixels per mm
center_x = 0 * 1e-3;
center_y = 0 * 1e-3;

% X = (-(x_range*res_factor/2-1):(x_range*res_factor/2))/res_factor - center_x; 
% Y = (-(y_range*res_factor/2-1):(y_range*res_factor/2))/res_factor - center_y;

X = (-255:256)*0.2-0; % FOV
Y = (-255:256)*0.2+0;

X = X * 1e-3; % FOV
Y = Y * 1e-3;

% v = linspace(1.486,1.502,6);
% v = 1.483; % water speed of sound [km/s]
v = 1.4899e3; % water speed of sound [m/s]

%%
% v = 1.483; % water speed of sound [km/s]

for i = 1:length(delay)
% for i = 1:length(v)

tic
% reIMG3 = subfunc_3d_cuda_reduce(DAQ_all,X,Y,Z,v(i),delay,xk,yk,zk,fs);
reIMG2 = subfunc_2d_cuda_reduce(DAQ_all,X*1e3,Y*1e3,v/1e3,delay(i),xk*1e3,yk*1e3,fs/1e6);
toc

figure('Name',string(delay)); 
XY_view = squeeze(max((reIMG2(:,:)),[],3));
imshow(XY_view,[]) % Top view, XY
xlabel ('--y++')
ylabel ('--x++')

% [mx(i),ind] = max(reIMG3(:));
end

%% display
% dynamic_range = [-10,0];
% imagedata = abs(reIMG2) ./ max(abs(reIMG2(:)));
% imagelog = 20 * log10(imagedata);
% figure();
% % imagesc(Display.xm(x_ind) * 1e2,Display.ym(y_ind) * 1e2, imagelog(y_ind,x_ind), dynamic_range);colormap(gray); colorbar;
% imagesc(imagelog, dynamic_range);colormap(gray); colorbar;
% % xlabel('X (cm)');ylabel('Y (cm)');axis('image');
