function [raw_data, Acq]= f_load_data(Para)

load([Para.datapath,Para.dataname]);
Acq.fs = 5e6; % sampling rate
Acq.Nf = 512;% number of frames

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
        %raw_data(:,:,i) = hilbert(bandpass(raw_data(:,:,i), [0.3e6,1.3e6],Acq.fs));
        raw_data(:,:,i) = abs(hilbert(raw_data(:,:,i)));
    end
    %raw_data = cdata;
end

end




