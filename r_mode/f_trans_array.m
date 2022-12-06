

function Trans = f_trans_array(Trans)

calibFile = 1;
% Detection array

if calibFile
    load('../calibrated_coordinates/calibrated_coords_2.mat')
    R = R_calibrated; 
    
    n_receive = 1:512;
    n_receive(fault_id) = [];
    x_receive = tx_x(n_receive,2).';
    y_receive = tx_x(n_receive,1).';
else
    n_receive = 1:512;
    R = 610.8 * 1e-3/2;    %mm, For Breast transducer
    ini_angle = 0; % the angle of the first step (in degree)
    Nstep = 512; %Trans.relements;
    angle_per_step = 2*pi/Nstep;   % angle per step, 1.5 degree/step (in radian)
    angle_step1 = ini_angle/180*pi;
    detector_angle = (0:Nstep-1)*angle_per_step+angle_step1;

    x_receive = cos(detector_angle)*R; %xk./1000; %

    %% CHANGED SIGN
    y_receive = sin(detector_angle)*R; %yk./1000; %
end

Trans.x_receive = x_receive;
Trans.y_receive = y_receive;
Trans.valid_id = n_receive; 
Trans.relements = length(n_receive);

%%

z_receive = zeros(1,length(x_receive));
Trans.rxyz = single(cat(1,x_receive, y_receive, z_receive));

%% Transmission array
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

%% CHANGED SIGN
x_transmit = -cos(transmit_angle)*Rt + Trans.x_offset;

%% 
y_transmit = sin(transmit_angle)*Rt + Trans.y_offset;


%%

z_transmit= zeros(1,Nstep);
Trans.txyz = single(cat(1,x_transmit, y_transmit, z_transmit));

% detection apodization
Trans.r_Napo = round((Trans.r_Aapo/360) * Trans.relements); % no. of elements
% temp_center = linspace(1,Trans.relements,Trans.t_Nsteps);% center of apo = center of transmission position

r_list = linspace(1,Trans.relements,Trans.relements);
temp_center = round(transmit_angle./(2*pi) * Trans.relements); % center for 
Trans.rapo_center = temp_center;

temp_alist = round(Trans.relements/2 - Trans.r_Napo/2) : round(Trans.relements/2+Trans.r_Napo/2);
temp_alist = temp_alist(1 : Trans.r_Napo);% for center apo range

Trans.rapo_list = zeros(Trans.relements,Trans.t_Nsteps);
for i = 1 : Trans.t_Nsteps
    temp0_list = circshift(r_list, round(Trans.relements/2));% center if no bias of initial angle
    temp_list = circshift(temp0_list , -temp_center(i)); % 133-387-512-120 (initial = 387)
    temp_apo_list = temp_list(temp_alist);
    temp_apo = tukeywin(length(temp_apo_list),0.25); % add tukey win for suppression of side lobes
    Trans.rapo_list(temp_apo_list,i) = temp_apo; % uniform apodization
end

end

% Note, check the Trans.rapo_list and Trans.rapo_center for debug.

