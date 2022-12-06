
%% load data
clc
clear 
% close all
Para.datapath = './';
% Para.dataname = 'glycerin_air';
% Para.dataname = 'water_only';

%Para.dataname = '../phantom_2022-12-01 11-22-26_avg-remap';
Para.dataname = '../Data 2022-12-01 19-55-30 avg-remap'; %

% Para.dataname = 'glycerin_air_shifted';
% Para.dataname = 'glycerin_alcohol_shifted';
Para.display = 0; % Display raw signal or not
Para.hilbert = 0; % single to complex single

[raw_data, Acq]= f_load_data(Para);% 512 x 4000 x 320


%%
% figure; 
% for i=1:435
%     clim = 2048;
%     imagesc(squeeze(abs(raw_data(:,:,i))),[-clim,clim])
%     title(['Frame: ', num2str(i)])
%     drawnow;
%     pause(0.05)
% end

%% Para setting

% transducers
%Trans.tfoc_bias = 85 * 1e-3; % distance between surface and arc, for g,airs
Trans.tfoc_bias = 72.7e-3; %72 * 1e-3; % distance between surface and arc, % for air-air

Trans.t_init_angle =90; %180; % degrees, % -87

Trans.t_Nsteps = 403; %433; %595; %Acq.Nf;% no. of positions = no. of frames
Trans.tfoclens = 1.00 * 25.4 * 1e-3;% focus length
Trans.tclockwise = 1; % rotation direction, Tbc.
Trans.r_Aapo = 180; %  (degrees) Angle range of detection array elements for recon
Trans.x_offset = 0e-3; %-5e-3; %0; %12e-3;
Trans.y_offset = 0e-3; %-5e-3; %17e-3;
Trans.scan_angle = 2*pi; % 2*pi
Trans = f_trans_array(Trans); % Trans.rapo_list;

%raw_data = raw_data(Trans.valid_id,:,:);

% display roi
Display.x_range = 80 * 1e-3; % m
Display.y_range = 80 * 1e-3; % m
Display.z_range = 0 * 1e-3; % m
Display.res_factor = 3; % resolution factor: no. of pixels per mm
Display.center_x = 0*1e-3;
Display.center_y = 0*1e-3;
Display.center_z = 0 * 1e-3; % m

% acquisition para
Acq.c = 1.483e3; %1.495 * 1e3;% speed of sound (m/s)
Acq.delay = 0; %-(Trans.tfoclens)/Acq.c; % trigger delay of DAQ between US transmission
Acq.delay_revise_m =(0)*1e-6 * Acq.c + Acq.delay; % DAQ delay , 25 sampling error of DAQ = 0.625us

%% make an image
R = 610.8 * 1e-3/2;    %mm, For Breast transducer
ini_angle = 0; % the angle of the first step (in degree)
Nstep = 512; %Trans.relements;
angle_per_step = 2*pi/Nstep;   % angle per step, 1.5 degree/step (in radian)
angle_step1 = ini_angle/180*pi;
detector_angle = (0:Nstep-1)*angle_per_step+angle_step1;

x_receive = Trans.x_receive; %cos(detector_angle)*R; %xk./1000; %
y_receive = Trans.y_receive; %sin(detector_angle)*R; %yk./1000; %
N_receive = length(x_receive);

Rt = R - Trans.tfoclens - Trans.tfoc_bias; % radius of the transmission focus
Nstep = Trans.t_Nsteps;
angle_per_step = Trans.scan_angle/Nstep;   % angle per step, 1.5 degree/step (in radian)
angle_step1 = Trans.t_init_angle/180*pi;
if Trans.tclockwise == 1
transmit_angle = (0:Nstep-1)*angle_per_step + angle_step1;
else
transmit_angle = (Nstep-1:-1:0)*angle_per_step - angle_step1; % minus? Tbc.
end

% transmit_angle (0,2*pi)
transmit_angle(transmit_angle<0) = transmit_angle(transmit_angle<0) + 2 * pi; % negative 
transmit_angle(transmit_angle>2*pi) = transmit_angle(transmit_angle>2*pi) - 2 * pi;

x_transmit = -cos(transmit_angle)*Rt + Trans.x_offset;
y_transmit = sin(transmit_angle)*Rt + Trans.y_offset;


% 5 MHz sampling
dt = 1/5e6;
t_vec = 0:dt:dt*3999;

% first pos
toa = zeros(N_receive,1);

delay_add = -(Trans.tfoclens)/Acq.c;

%%
figure;
for j = 1:403
    for i=1:N_receive
        toa(i) = (sqrt( (x_receive(i) - x_transmit(j)).^2 + (y_receive(i) - y_transmit(j)).^2) ./ Acq.c) - delay_add;
    end

    toa_ind = round(toa./dt);
    toa_img = zeros(N_receive,4000);

    for i=1:N_receive
        toa_img(i,toa_ind(i):toa_ind(i)+5) = 1;
    end

%     figure;
%     imagesc(toa_img)

    img_ex = squeeze(abs(raw_data(:,:,j)));
    img_ex = img_ex./max(max(img_ex));

    img_over = img_ex*4 - toa_img;


    clim = 2048;
    imagesc(img_over,[-1,1])
    drawnow
    pause(0.02)
end

%% Recon

tic
frame_list = 1 : Trans.t_Nsteps;
n = 1;

plotOption = 0;
figure; hold on;

for i = 10:10:400 %i = frame_list

    Trans.tfocus = Trans.txyz(:,i);% position of transmission focus
    Trans.apo_receive = Trans.rapo_list(:,i);

    disp(['recon:',num2str(frame_list(i)),'/',num2str(Trans.t_Nsteps)]);

    % disp(' start obtaining recon parameters..');
    [Reconpara, Display]= f_recon_rusct_para(Acq,Trans,Display);
    Reconpara.map_epim_xoy(Reconpara.map_epim_xoy < 1) = 1;% remove zeros points

    % disp(' start reconstruction.. ');
    recondata = f_recon_rusct(Trans,Display,Reconpara,squeeze(raw_data(:,:,frame_list(i))));
    % disp('reconstructed.');

    if n == 1
        ReconData = zeros(size(recondata,1),size(recondata,2),length(frame_list));
        recon_compound = zeros(size(recondata,1),size(recondata,2));
    end
    ReconData(:,:,n) = recondata;
    recon_compound = recon_compound + recondata;
    n = n + 1;
    
    
    if plotOption && not(mod(i,5))
        dynamic_range = [-35,0];
        %imagedata = abs(recon_compound) ./ max(abs(recon_compound(:)));
        imagedata = abs(recondata) ./ max(abs(recondata(:)));
        imagelog = 20 * log10(imagedata);
        imagesc(Display.xm * 1e2,Display.ym * 1e2,imagelog,dynamic_range);colormap(gray); colorbar;
        title(num2str(i));
        pause();
    end
    
    
end
toc

%% display
% x_ind = 76:115; y_ind = 201:240;
x_ind = 1:240; y_ind = 1:240;

dynamic_range = [-10,0];
imagedata = abs(recon_compound) ./ max(abs(recon_compound(:)));
imagelog = 20 * log10(imagedata);
figure();
imagesc(Display.xm(x_ind) * 1e2,Display.ym(y_ind) * 1e2, imagelog(y_ind,x_ind), dynamic_range);colormap(gray); colorbar;
% imagesc(imagelog, dynamic_range);colormap(gray); colorbar;
xlabel('X (cm)');ylabel('Y (cm)');axis('image');
title(Para.dataname)


figure;
imagesc(imagelog(y_ind,x_ind), dynamic_range);colormap(gray); colorbar;

% imagesc(Display.xm * 1e3,Display.ym * 1e3,imagelog,dynamic_range);colormap(gray); colorbar;
% xlabel('X (mm)');ylabel('Y (mm)');axis('image');

%%
ind = 93;
 figure; plot(Display.xm .* 1000, abs(recon_compound(:,ind))./max(abs(recon_compound(:,ind))));
 xlabel('X (mm)')
%
ind = 225;
 figure; plot(Display.xm .* 1000, abs(recon_compound(ind,:))./max(abs(recon_compound(ind,:))));
 xlabel('Y (mm)')
