function Display = f_set_display(Display)
% coordinate of each pixel
Display.Nx = round(Display.x_range * 1e3 * Display.res_factor) ;
Display.Ny = round(Display.y_range * 1e3 * Display.res_factor);
Display.Nz = round(Display.z_range * 1e3 * Display.res_factor);

Display.xm = ((1:Display.Nx)-Display.Nx/2)*Display.x_range/(Display.Nx)+Display.center_x;   % x axis coordinates
Display.ym = ((1:Display.Ny)-Display.Ny/2)*Display.y_range/(Display.Ny)+Display.center_y;   % y axis coordinates
Display.zm = ((1:Display.Nz)-Display.Nz/2)*Display.z_range/(Display.Nz)+Display.center_z;   % z axis coordinates

end