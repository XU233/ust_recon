% default SI units

clear all; clc

%% load data
%Para.datapath = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-12-01\';
% Para.dataname = 'Data 2022-12-01 19-55-30 avg-remap short'; 

Para.datapath = 'E:\OneDrive - California Institute of Technology\PhD\experiments\ust\2022-12-01\';
%Para.dataname = 'Data 2022-12-01 19-55-30 avg-remap'; 
Para.dataname = 'phantom_2022-12-01 11-22-26_avg-remap'; 

Para.display = 0; % Display raw signal or not
Para.hilbert = 0; % single to complex single

[raw_data, Acq]= f_load_data(Para);

% Data_all = raw_data(:,:,1:10);
% Data_all = sum(raw_data,3);

raw_data = single(raw_data);
% raw_data = single(raw_data - mean(raw_data,1));
raw_data = reshape( raw_data, 512, size(raw_data,2), numel(raw_data)/512/size(raw_data,2) );


%% transducer setting
Trans.t_foc_bias = 71.2e-3; %72 * 1e-3; % distance between surface and arc, % for air-air

Trans.t_init_angle = 90; %180; % degrees, % -87

Trans.t_Nsteps = 403; %433; %595; %Acq.Nf;% no. of positions = no. of frames
Trans.t_foclens = 1.00 * 25.4 * 1e-3;% focus length
Trans.t_clockwise = 0; % rotation direction, Tbc.
Trans.r_Aapo = 180; %  (degrees) Angle range of detection array elements for recon
Trans.x_offset = 0e-3; %-5e-3; %0; %12e-3;
Trans.y_offset = 0e-3; %-5e-3; %17e-3;
Trans.t_scan_angle = 2*pi; % 2*pi
Trans = f_trans_array(Trans); % Trans.rapo_list;

%% display field of view
Display.x_range = 128 * 2 * 1e-3; % m
Display.y_range = 128 * 2 * 1e-3; % m
Display.z_range = 0 * 1e-3; % m
Display.res_factor = 8; % resolution factor: no. of pixels per mm
Display.center_x = 0 * 1e-3;
Display.center_y = 0 * 1e-3;
Display.center_z = 0 * 1e-3; % m

%% acquisition para
Acq.c = 1.483e3; %1.495 * 1e3;% speed of sound (m/s)
Acq.delay = -Trans.t_foclens/Acq.c; % trigger delay of DAQ between US transmission
Acq.delay_revise_m =(0)*1e-6 * Acq.c + Acq.delay; % DAQ delay , 25 sampling error of DAQ = 0.625us

%% show sinogram
% f_show_sinogram(Trans, Acq, raw_data)

%% set display
Display = f_set_display(Display);

%% scale data by apodization

raw_data_used = single(zeros(size(raw_data)));
for i=1:Trans.t_Nsteps
    raw_data_used(:,:,i) = raw_data(:,:,i) .* Trans.r_apo_list(:,i);
end


%% reconstruction
tic
% for i = 1:10
%     
%     Trans.t_focus = Trans.t_xyz(:,i);
%     Trans.apo_receive = Trans.r_apo_list(:,i);
%     
%     stat_disp = fprintf('recon: %i / %i \n', i, Trans.t_Nsteps);
%     
%     [Reconpara, Display]= f_recon_rusct_para(Acq,Trans,Display);
%     Reconpara.map_epim_xoy(Reconpara.map_epim_xoy < 1) = 1; % remove zeros points
%     
%     reIMG2 = subfunc_2d_cuda_reduce(raw_data(:,:,i), Acq.c/1e3, Acq.delay*Acq.fs, Acq.fs/1e6, ...
%         Trans.x_receive*1e3, Trans.y_receive*1e3, Display.xm*1e3, Display.ym*1e3, Reconpara.map_epim_xoy);
% 
%     ReconData = zeros(size(reIMG2,1),size(reIMG2,2),Trans.t_Nsteps);
%     recon_compound = zeros(size(reIMG2,1),size(reIMG2,2));
%     ReconData(:,:,i) = reIMG2;
%     recon_compound = recon_compound + reIMG2;
%     
% end

reIMG_final = zeros(length(Display.xm), length(Display.ym));

plotOpt = 0; 
for i=1:403
    raw_data_used_i = raw_data_used(:,:,i);

%     reIMG2 = subfunc_2d_cuda_reduce(raw_data_used, Acq.c/1e3, Acq.delay*Acq.fs, Acq.fs/1e6, ...
%         Trans.x_receive*1e3, Trans.y_receive*1e3, Trans.x_transmit*1e3, Trans.y_transmit*1e3,...
%         Display.xm*1e3, Display.ym*1e3);

    reIMG2 = subfunc_2d_cuda_reduce(raw_data_used_i, Acq.c/1e3, Acq.delay*Acq.fs, Acq.fs/1e6, ...
        Trans.x_receive*1e3, Trans.y_receive*1e3, Trans.x_transmit(i)*1e3, Trans.y_transmit(i)*1e3,...
        Display.xm*1e3, Display.ym*1e3);

    if plotOpt
        imshow(reIMG2,[]) % Top view, XY
        xlabel ('--y++')
        ylabel ('--x++')
        drawnow;
        pause(0.01)
    end
    
    reIMG_final = reIMG_final + reIMG2;
end

    
toc



%% display
dynamic_range = [-80,0];
imagedata = abs(reIMG_final) ./ max(abs(reIMG_final(:)));
imagelog = 20 * log10(imagedata);
%imagelog = 20 * log(imagedata);

figure();
% imagesc(Display.xm(x_ind) * 1e2,Display.ym(y_ind) * 1e2, imagelog(y_ind,x_ind), dynamic_range);colormap(gray); colorbar;
imagesc(imagelog.', dynamic_range);colormap(gray); colorbar;
axis equal;
% xlabel('X (cm)');ylabel('Y (cm)');axis('image');
