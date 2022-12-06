


% raw_data 512 x 4000 x 320


function [raw_data, Acq]= f_load_data(Para)

load([Para.datapath,Para.dataname]);
Acq.fs = 5e6; % sampling rate
Acq.Nf = 512;% number of frames

% load channel maps
% CH_1_128 = load('.\ch_map\array2DAQ_1').array2DAQ_1(:,1:2);
% CH_129_256 = load('.\ch_map\array2DAQ_129').array2DAQ_129(:,1:2);
% CH_257_384 = load('.\ch_map\array2DAQ_257').array2DAQ_257(:,1:2);
% CH_385_512 = load('.\ch_map\array2DAQ_385').array2DAQ_385(:,1:2);
% 
% [N_ch,N_samp,N_fr] = size(VOLTAGE);
% VOLTAGE_remap = zeros(N_ch,N_samp,N_fr);
% 
% VOLTAGE_remap(CH_257_384(:,1),:,:) = VOLTAGE(1:128,:,:);
% VOLTAGE_remap(CH_129_256(:,1),:,:) = VOLTAGE(129:256,:,:);
% 
% VOLTAGE_remap(CH_1_128(:,1),:,:) = VOLTAGE((1:128)+256,:,:);
% VOLTAGE_remap(CH_385_512(:,1),:,:) = VOLTAGE((129:256)+256,:,:);

raw_data = double(VOLTAGE_remap);
clear VOLTAGE_remap

if Para.display == 1
    fs = 5 * 1e6;
    sos = 1500; % m/s
    d = 1 : size(raw_data,2)./fs * sos * 1e2;% mm
    ch = 1 : size(raw_data,1);
    figure(1001)
    for i = 1 : 5: size(raw_data,3)
%     for i = 1
        %    plot(squeeze(VOLTAGE(1,:,i)));
        %     imshow(squeeze(raw_data(:,:,i)),[-1,1]*1000)
        %     imagesc(squeeze(raw_data(:,:,i)),[-1,1]*3000);
        %     xlabel('Samples'),ylabel('Amplitude');
        imagesc(d,ch,abs(squeeze(raw_data(:,:,i))),[-1,1]*2000);
        xlabel('Distance (cm)'),ylabel('Channel number');
        %     axis('image');
        title(num2str(i))
        pause(0.05)
    end
end

if Para.hilbert == 1
    %cdata = single(raw_data);% data: 8192 * 512 * 7
    clear data
%     temp = single(zeros(200,size(cdata,2),size(cdata,3)));% add zeros
%     cdata = cat(1,cdata,temp);
    for i = 1:size(raw_data,3)
        raw_data(:,:,i) = hilbert(bandpass(raw_data(:,:,i), [0.3e6,1.3e6],Acq.fs));
        %raw_data(:,:,i) = abs(hilbert(raw_data(:,:,i)));
    end
    %raw_data = cdata;
end

end




