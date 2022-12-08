function reIMG = subfunc_2d_cuda_reduce(data,v1,fixedDelay,fs,xk,yk,X,Y,map_epim)
% data: raw data axes transducer, time, frame
% v1: speed of sound (km/s)
% fixed_delay: index of surface signal relative to start of data (index #)
% fs: sampling frequency (MHz)
% xk: x coordinates of receiving transducers(mm)
% yk: y coordinates of receiving transducers(mm)
% X: x coordinates of reconstructed frame (mm units)
% Y: y coordinates of reconstructed frame (mm units)
% map_epim: grid of distances from the transmitter to the pixels in reconstructed frame (index #)



tic
% [Nt, Nk, Nphi] = size(data);   % Nstep : number of steps, Nsample : number of samples per step   
[Nk, Nt, Nphi] = size(data);   % Nstep : number of steps, Nsample : number of samples per step   
%fs = 40;                            % sampling frequency (MHz)
%fs = 1/kgrid.dt/10^6;
tv = fs/v1;



%parameters of image
Nx = numel(X);
Ny = numel(Y);
Nz = numel(0);

% receiver position
xk = gpuArray(single(xk));
yk = gpuArray(single(yk));
zk = gpuArray(single(0));

% image position
X = gpuArray(single(X));
Y = gpuArray(single(Y));
Z = gpuArray(single(0));

map_epim = gpuArray(single(map_epim));


reIMG = gpuArray.zeros(Nx,Ny,Nz,'single');
data = padarray(data,[0,1,0],'both');
data = gpuArray(data);

    

k = parallel.gpu.CUDAKernel('pct_ubp_reduce.ptx','pct_ubp_reduce.cu');
k.ThreadBlockSize = [Nk 1 1];
k.GridSize = [Nx*Nphi,Ny,Nz];
reIMG = feval(k,reIMG,X,Y,Z,xk(:),yk(:),zk(:),Nx,uint16(log2(Nx)),Nphi,Nt+2,fixedDelay,map_epim,tv,data);

reIMG = gather(reIMG);
%reIMG = reIMG./weights;
%     k = parallel.gpu.CUDAKernel('pct_ubp.ptx','pct_ubp.cu');
%     k.ThreadBlockSize = [Nphi 1 1];
%     k.GridSize = [Nx,Ny,Nz];
% 
%     
%     reIMG = feval(k,reIMG,X,Y,Z,xk(:),yk(:),zk(:),Nphi,Nk,Nt+2,fixedDelay,tv,data);

end
