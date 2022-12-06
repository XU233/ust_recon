% clearvars
% close all

path = 'C:\Users\Legion12\Downloads\sdk1.1 20220429\SDK 1.1\examples\matlab\Data 2022-12-01\';
filename = ['Data 2022-12-01 19-55-30_avg-remap'];

Data_all = load([path, filename, '.mat']).VOLTAGE_remap;

Data_all = Data_all(:,:,1:10);

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
% DAQ_all = Data_all_f;
DAQ_all = Data_all;
DAQ_all = single(DAQ_all - mean(DAQ_all,1));
DAQ_all = permute(DAQ_all,[2,1,3]);

%%
DAQ_all = reshape(DAQ_all,size(DAQ_all,1),512,numel(DAQ_all)/512/size(DAQ_all,1));


%% 
X = (-255:256)*0.2-0; % FOV
Y = (-255:256)*0.2+0;

%
% v = linspace(1.486,1.502,6);
v = 1.483; % water speed of sound [km/s]
% v = 1.525;
% v = 1.471;

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

