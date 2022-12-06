
%% load data
clc
clear 
% close all
Para.datapath = './';
% Para.dataname = 'glycerin_air';
% Para.dataname = 'water_only';

%Para.dataname = 'neck_2_remap';
%Para.dataname = 'abdomen_3_remap'; %
%Para.dataname = '../water_only_30dB_remap'; %;
%Para.dataname = '../line_damped_30dB_remap'; %'../dsiplaced&damped_30dB_remap';
Para.dataname = '../dsiplaced&damped_30dB_remap';
%Para.dataname = '../neck_1_remap';

% Para.dataname = 'glycerin_air_shifted';
% Para.dataname = 'glycerin_alcohol_shifted';
Para.display = 1; % Display raw signal or not
Para.hilbert = 1; % single to complex single

[raw_data, Acq]= f_load_data(Para);% 512 x 4000 x 320

%% Para setting

x_off_vec = 0; %-20e-3:10e-3:20e-3;
y_off_vec = -20e-3; %-20e-3:10e-3:20e-3;
angle_vec = 190; %150:10:210;

for xi = x_off_vec
    for yi = y_off_vec
        for ai = angle_vec

            Trans.x_offset = xi; %0; %12e-3;
            Trans.y_offset = yi; %17e-3;
            Trans.t_init_angle = ai; % degrees, % -87

            % transducers
            % Trans.tfoc_bias = 85 * 1e-3; % distance between surface and arc, for g,airs
            Trans.tfoc_bias = 65 * 1e-3; % distance between surface and arc, % for air-air



            Trans.t_Nsteps = 595; %Acq.Nf;% no. of positions = no. of frames
            Trans.tfoclens = 0.75 * 25.4 * 1e-3;% focus length
            Trans.tclockwise = 1; % rotation direction, Tbc.
            Trans.r_Aapo = 90; %  (degrees) Angle range of detection array elements for recon

            Trans = f_trans_array(Trans); % Trans.rapo_list;

            % display roi
            Display.x_range = 300 * 1e-3; % m
            Display.y_range = 300 * 1e-3; % m
            Display.z_range = 0 * 1e-3; % m
            Display.res_factor = 0.5; % resolution factor: no. of pixels per mm
            Display.center_x = 0*1e-3;
            Display.center_y = 0*1e-3;
            Display.center_z = 0 * 1e-3; % m

            % acquisition para
            Acq.c = 1.505e3; %1.495 * 1e3;% speed of sound (m/s)
            Acq.delay = 0; %-(Trans.tfoc_bias)/Acq.c; % trigger delay of DAQ between US transmission
            Acq.delay_revise_m =(0)*1e-6 * Acq.c + Acq.delay * Acq.c; % DAQ delay , 25 sampling error of DAQ = 0.625us

            %% Recon

            tic
            frame_list = 1 : Trans.t_Nsteps;
            n = 1;

            plotOption = 0;
            figure; hold on;

            for i = frame_list

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


                if plotOption && not(mod(i,10))
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

            % display

            dynamic_range = [-25,0];
            imagedata = abs(recon_compound) ./ max(abs(recon_compound(:)));
            imagelog = 20 * log10(imagedata);
            figure();
            imagesc(Display.xm * 1e2,Display.ym * 1e2,imagelog,dynamic_range);colormap(gray); colorbar;
            xlabel('X (cm)');ylabel('Y (cm)');axis('image');
            title([num2str(xi), num2str(yi), num2str(ai)])
            savefig(['./', num2str(xi), num2str(yi), num2str(ai),'.fig'])
            
            save(['./', num2str(xi), num2str(yi), num2str(ai),'.mat'],'recon_compound')
            % imagesc(Display.xm * 1e3,Display.ym * 1e3,imagelog,dynamic_range);colormap(gray); colorbar;
            % xlabel('X (mm)');ylabel('Y (mm)');axis('image');


            %figure; plot(Display.xm .* 1000, abs(recon_compound(:,75))./max(abs(recon_compound(:,75))));
            %xlabel('X (mm)')

            %figure; plot(Display.xm .* 1000, abs(recon_compound(116,:))./max(abs(recon_compound(116,:))));
            %xlabel('Y (mm)')
        end
    end
end
