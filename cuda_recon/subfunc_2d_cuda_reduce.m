function reIMG = subfunc_3d_cuda_reduce(data,X,Y,v1,fixedDelay,xk,yk,fs)
% data: raw data axes time, transducer, frame
% X: x coordinates of reconstructed frame, size Ny by Nx by Nz (mm units)
% Y: y coordinates of reconstructed frame, size Ny by Nx by Nz (mm units)
% Z: z coordinates of reconstructed frame, size Ny by Nx by Nz (mm units)
% v1: speed of sound (km/s)
% fixed_delay: index of surface signal relative to start of data
% dphi: angle between acquisitions (radians)


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


X = gpuArray(single(X));
Y = gpuArray(single(Y));
Z = gpuArray(single(0));
xk = gpuArray(single(xk));
yk = gpuArray(single(yk));
zk = gpuArray(single(0));


reIMG = gpuArray.zeros(Nx,Ny,Nz,'single');
data = padarray(data,[0,1,0],'both');
data = gpuArray(data);

    

k = parallel.gpu.CUDAKernel('pct_ubp_reduce.ptx','pct_ubp_reduce.cu');
k.ThreadBlockSize = [Nk 1 1];
k.GridSize = [Nx*Nphi,Ny,Nz];
reIMG = feval(k,reIMG,X,Y,Z,xk(:),yk(:),zk(:),Nx,uint16(log2(Nx)),Nphi,Nt+2,fixedDelay,tv,data);

reIMG = gather(reIMG);
%reIMG = reIMG./weights;
%     k = parallel.gpu.CUDAKernel('pct_ubp.ptx','pct_ubp.cu');
%     k.ThreadBlockSize = [Nphi 1 1];
%     k.GridSize = [Nx,Ny,Nz];
% 
%     
%     reIMG = feval(k,reIMG,X,Y,Z,xk(:),yk(:),zk(:),Nphi,Nk,Nt+2,fixedDelay,tv,data);

end
