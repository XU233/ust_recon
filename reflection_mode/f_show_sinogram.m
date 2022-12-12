function f_show_sinogram(Trans, Acq, raw_data)
R = 610.8 * 1e-3/2;    %mm, For Breast transducer
ini_angle = 0; % the angle of the first step (in degree)
Nstep = 512; %Trans.relements;
angle_per_step = 2*pi/Nstep;   % angle per step, 1.5 degree/step (in radian)
angle_step1 = ini_angle/180*pi;
detector_angle = (0:Nstep-1)*angle_per_step+angle_step1;

x_receive = Trans.x_receive; %cos(detector_angle)*R; %xk./1000; %
y_receive = Trans.y_receive; %sin(detector_angle)*R; %yk./1000; %
N_receive = length(x_receive);

Rt = R - Trans.t_foclens - Trans.t_foc_bias; % radius of the transmission focus
Nstep = Trans.t_Nsteps;
angle_per_step = Trans.t_scan_angle/Nstep;   % angle per step, 1.5 degree/step (in radian)
angle_step1 = Trans.t_init_angle/180*pi;
if Trans.t_clockwise == 1
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

delay_add = -(Trans.t_foclens)/Acq.c;

%%
figure;
for j = 1:size(raw_data,3)
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
end