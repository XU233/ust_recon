
%{ 
Aim: Reconstruction parameters for the 4 X 256 curved array

Display.x_range = 25; % mm
Display.y_range = 25; % mm
Display.res_factor = 15; % resolution factor
Display.center_x = 0;
Display.center_y = 0;

Trans.tfocus = [0,0,25.4] * 1e-3;% coordination of focus
Trans.tfoclens = 25.4 * 1e-3;% focus length
Trans.rxyz = zeros(3,512);% x, y, z for each elements
Trans.relements = 512;

Init.c = 1540; % m/s
Init.fs = 40 * 1e6;% DAQ sampling
Init.delay_revise = 50*1e-6; % DAQ delay

%}

function [Reconpara, Display] = f_recon_rusct_para(Acq,Trans,Display)


% % coordinate of each pixel
Display.Nx = round(Display.x_range * 1e3 * Display.res_factor) ;
Display.Ny = round(Display.y_range * 1e3 * Display.res_factor);
Display.Nz = round(Display.z_range * 1e3 * Display.res_factor);

Display.xm = ((1:Display.Nx)-Display.Nx/2)*Display.x_range/(Display.Nx)+Display.center_x;% x axis coordinates
Display.ym = ((1:Display.Ny)-Display.Ny/2)*Display.y_range/(Display.Ny)+Display.center_y; % y axis coordinates
Display.zm = ((1:Display.Nz)-Display.Nz/2)*Display.z_range/(Display.Nz)+Display.center_z; % y axis coordinates
% Display.zm = ((Display.Nz:-1:1)-Display.Nz/2)*Display.z_range/(Display.Nz)+Display.center_z; % y axis coordinates

%% xoy reconstruction parameters
[image_xm,image_ym] = meshgrid(Display.xm,Display.ym); % coordinates of all pixels in m

% transmit delay for each pixel
%% CHANGED INDEX
d_emit = sqrt((image_xm - Trans.tfocus(1)).^2 + (image_ym - Trans.tfocus(2)).^2 + Trans.tfocus(3)^2);% distance from source
%%
d_emit = d_emit + Trans.tfoclens; % consider transducer focus

% plus reception delay for each pixel
d_receive = zeros(Display.Ny, Display.Nx,Trans.relements);
d_delay = zeros(Display.Ny, Display.Nx,Trans.relements);
for i = 1 : Trans.relements
    d_receive(:,:,i) = sqrt((image_xm - Trans.rxyz(1,i)).^2 + (image_ym - Trans.rxyz(2,i)).^2 + Trans.rxyz(3,i).^2);
    d_delay(:,:,i) = d_receive(:,:,i) + d_emit; 
%     d_delay(:,:,i) = d_receive(:,:,i)+Trans.tfoclens;
end
d_delay = d_delay - Acq.delay_revise_m;% consider the delay bias;

% positioins or profiles in # of channel data
map_epim = ceil(d_delay/Acq.c*Acq.fs);

Reconpara.map_epim_xoy = single(map_epim);


% solid angle weight
%{
dr = 0;
for i = 1:Trans.relements
    r0 = sqrt(Trans.rxyz(1,i)^2+Trans.rxyz(2,i)^2);
    dx = image_xm-Trans.rxyz(1,i);
    dy = image_ym-Trans.rxyz(2,i);
    rr0 = sqrt(dx.^2+dy.^2+dr^2);         % distance to the detector
    dsx = image_xm-Trans.tfocus(1);
    dsy = image_ym-Trans.tfocus(2);
    dsz = Trans.tfocus(3)+dr;
    rs0 = sqrt(dsx.^2+dsy.^2+dsz^2);      % distance to the ps
    % normalize the reconstruction at each location by a total solid angle
    % alpha = acos(min(abs((-x_receive(iStep)*dx-y_receive(iStep)*dy)/r0./rr0),0.999))+1.0e-5;
    alpha = acos(abs((-Trans.rxyz(1,i)*dx-Trans.rxyz(2,i)*dy)/r0./rr0));
    angle_weight(:,:,i) = real(cos(alpha))./rr0.^2;
end

% total_angle_weight = sum(angle_weight,3);
total_angle_weight = max(angle_weight,[],3);
temp_weight = repmat(total_angle_weight, [1,1,Trans.relements]);

Reconpara.sg_weight = single(angle_weight ./ temp_weight);% final solid angle
%}
% pause

%% yoz reconstruction paramters
% [image_ym,image_zm] = meshgrid(Display.ym,Display.zm); % coordinates of all pixels in m
% 
% % transmit delay for each pixel
% d_emit = sqrt((image_ym - Trans.tfocus(2)).^2 + (image_zm - Trans.tfocus(3)).^2 + Trans.tfocus(1)^2);% distance from source
% d_emit = d_emit + Trans.tfoclens; % consider transducer focus
% 
% % plus reception delay for each pixel
% d_receive = zeros(Display.Nz, Display.Ny,Trans.relements);
% d_delay = zeros(Display.Nz, Display.Ny,Trans.relements);
% for i = 1 : Trans.relements
%     d_receive(:,:,i) = sqrt((image_ym - Trans.rxyz(2,i)).^2 + (image_zm - Trans.rxyz(3,i)).^2 + Trans.rxyz(1,i).^2);
%     d_delay(:,:,i) = d_receive(:,:,i) + d_emit; % only reception if com
% %     d_delay(:,:,i) = d_receive(:,:,i)+Trans.tfoclens;
% end
% d_delay = d_delay - Init.delay_revise_m;% consider the delay bias;
% 
% % positioins or profiles in # of channel data
% map_epim = ceil(d_delay/Init.c*Init.fs);
% 
% Reconpara.map_epim_yoz = single(map_epim);

%% xyz 3D reconstruction parameters

% [image_xm,image_ym,image_zm] = meshgrid(Display.xm,Display.ym,Display.zm); % coordinates of all pixels in m
% 
% % transmit delay for each pixel
% d_emit = sqrt((image_xm - Trans.tfocus(1)).^2 + (image_ym - Trans.tfocus(2)).^2 + (image_zm - Trans.tfocus(3)).^2);% distance from source
% d_emit = d_emit + Trans.tfoclens; % consider transducer focus
% 
% % plus reception delay for each pixel
% d_receive = zeros(Display.Ny, Display.Nx,Display.Nz,Trans.relements);
% d_delay = zeros(Display.Ny, Display.Nx,Display.Nz,Trans.relements);
% for i = 1 : Trans.relements
%     d_receive(:,:,:,i) = sqrt((image_xm - Trans.rxyz(1,i)).^2 + (image_ym - Trans.rxyz(2,i)).^2 + (image_zm - Trans.rxyz(3,i)).^2);
%     d_delay(:,:,:,i) = d_receive(:,:,:,i) + d_emit; 
% %     d_delay(:,:,:,i) = d_receive(:,:,:,i)+Trans.tfoclens; 
%     
% end
% d_delay = d_delay - Acq.delay_revise_m;% consider the delay bias;
% 
% % positioins or profiles in # of channel data
% map_epim = ceil(d_delay/Acq.c*Acq.fs);
% 
% Reconpara.map_epim_xyz = single(map_epim);

end

