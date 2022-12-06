

%rawdata = zeros(1000,512);
%Trans.apo_receive = ones(512,1);

function Recondata = f_recon_rusct(Trans,Display,Reconpara,rawdata)

% xoy recon
recon_data_xoy = zeros(Display.Ny,Display.Nx);
recon_data_xoy = single(recon_data_xoy);
for k = 1 : Trans.relements
% temp_data = rawdata(:,k);
temp_data = rawdata(k,:); % correct
temp_map = Reconpara.map_epim_xoy(:,:,k);
temp_map(temp_map > 4000) = 4000;
temp_map(temp_map < 1) = 1;
map_data_each = temp_data(temp_map);
recon_data_xoy = recon_data_xoy + Trans.apo_receive(k) .* map_data_each;% * or /
% temp_sgweight = Reconpara.sg_weight(:,:,k);
% recon_data = recon_data + Trans.apo_receive(k) .* map_data_each .* temp_sgweight;% * or /
end
% Recondata.recon_xoy = recon_data_xoy;
Recondata= recon_data_xoy; % only recon xoy

% yoz recon
% recon_data_yoz = zeros(Display.Nz,Display.Ny);
% recon_data_yoz = single(recon_data_yoz);
% for k = 1 : Trans.relements
% temp_data = rawdata(:,k);
% temp_map = Reconpara.map_epim_yoz(:,:,k);
% map_data_each = temp_data(temp_map);
% recon_data_yoz = recon_data_yoz + Trans.apo_receive(k) .* map_data_each;% * or /
% % temp_sgweight = Reconpara.sg_weight(:,:,k);
% % recon_data = recon_data + Trans.apo_receive(k) .* map_data_each .* temp_sgweight;% * or /
% end
% Recondata.recon_yoz = recon_data_yoz;

% % % xyz recon
% recon_data_xyz = zeros(Display.Ny,Display.Nx,Display.Nz);
% recon_data_xyz = single(recon_data_xyz);
% for k = 1 : Trans.relements
% temp_data = rawdata(:,k);
% temp_map = Reconpara.map_epim_xyz(:,:,:,k);
% map_data_each = temp_data(temp_map);
% recon_data_xyz = recon_data_xyz + Trans.apo_receive(k) .* map_data_each;% * or /
% % temp_sgweight = Reconpara.sg_weight(:,:,k);
% % recon_data = recon_data + Trans.apo_receive(k) .* map_data_each .* temp_sgweight;% * or /
% end
% % Recondata.recon_xyz = recon_data_xyz;
% Recondata= recon_data_xyz;

end