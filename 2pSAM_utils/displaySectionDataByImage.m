function displaySectionDataByImage(sectionData,PSFParameters,x_label,y_label,title_name, lower_log, upper_log)
%% A Visualization Function to draw specific figures, yields a Surf figure 
%% x_label: choice of x-coordinate (1) time domain: 'x'  or  'y'; 
%%                                 (2) frequency domain :’kx' or 'ky'
%% y_label: choice of y-coordinate (1) time domain: 'y'  or  'z'; 
%%                                 (2) frequency domain :’ky' or 'kz'
%% title_name: picture title, a name string;
%%

%% time domian parameters
xysizePixels = PSFParameters.xysize; %xy point num
delta_xy = PSFParameters.dxy;   %nm, x,y direction FFT sampling period T 
T_xy = xysizePixels * delta_xy; %nm, side length of PSF slice, x,y direction FFT basic period T (periodical function )

zsizePixels = length(PSFParameters.defocus);  %z point num
delta_z  = PSFParameters.defocusInterval*1000;   %nm, z direction FFT sampling period T
T_z = zsizePixels*delta_z; %nm, height og PSF,z direction FFT basic function period T
%%

%% frequency domain parameters
delta_kxy = 2*3.14159/T_xy; % wave number(x,y direction FFT baseic angular frequency) 2*3.14/T_xy=w0
delta_kz = 2*3.14159/T_z; % z direction FFT baseic angular frequency w0
W_xy = xysizePixels * delta_kxy; % equal to the sampling angular frequency in x y direction 2*3.14/delta_xy
W_z = zsizePixels * delta_kz; % equal to the sampling angular frequency in z direction 2*3.14/delta_z
%%

%% light parameters
kMax = PSFParameters.duringGeneration.kMax;   %light wave number or angular frequency k= 2*3.14/wave length
kMaxPixels =  kMax/delta_kxy; % light max frequency point 
%%

  
%% find label
mum_ticks = 11;
if strcmp(x_label,'x')  ||  strcmp(x_label,'X') 
  x_label_name = 'X(nm)';
  x_ticks = linspace(-T_xy/2+delta_xy/2,+T_xy/2-delta_xy/2,mum_ticks);
  x_data = (linspace(1,xysizePixels,xysizePixels)-(xysizePixels+1)/2.0)*delta_xy;
end

if strcmp(x_label,'y') ||  strcmp(x_label,'Y') 
  x_label_name = 'Y(nm)';
  x_ticks = linspace(-T_xy/2+delta_xy/2,+T_xy/2-delta_xy/2,mum_ticks);
  x_data = (linspace(1,xysizePixels,xysizePixels)-(xysizePixels+1)/2.0)*delta_xy;
end

if strcmp(x_label,'kx') ||  strcmp(x_label,'Kx') 
  x_label_name = 'Kx(1/nm)';
  x_ticks = linspace(-W_xy/2+delta_kxy,+W_xy/2,mum_ticks);
  x_data = (linspace(1,xysizePixels,xysizePixels)-xysizePixels/2.0)*delta_kxy;
end

if strcmp(x_label,'ky') ||  strcmp(x_label,'Ky') 
  x_label_name = 'Ky(1/nm)';
  x_ticks = linspace(-W_xy/2+delta_kxy,+W_xy/2,mum_ticks);
  x_data = (linspace(1,xysizePixels,xysizePixels)-xysizePixels/2.0)*delta_kxy;
end


%%%%%%%%%%%%%%%%%%%
if strcmp(y_label,'y')  ||  strcmp(y_label,'Y')
  y_label_name = 'Y(nm)';
  y_ticks = linspace(-T_xy/2+delta_xy/2,+T_xy/2-delta_xy/2,mum_ticks);
  y_data = (linspace(1,xysizePixels,xysizePixels)-(xysizePixels+1)/2.0)*delta_xy;
end

if strcmp(y_label,'z') ||  strcmp(y_label,'Z')
  y_label_name = 'Z(nm)';
  y_ticks = linspace(-T_z/2+delta_z/2,+T_z/2-delta_z/2,mum_ticks);
  y_data = (linspace(1,zsizePixels,zsizePixels)-(zsizePixels+1)/2.0)*delta_z;
end

if  strcmp(y_label,'ky') ||  strcmp(y_label,'Ky')
  y_label_name = 'Ky(1/nm)';
  y_ticks = linspace(-W_xy/2+delta_kxy,+W_xy/2,mum_ticks);
  y_data = (linspace(1,xysizePixels,xysizePixels)-xysizePixels/2.0)*delta_kxy;
end

if  strcmp(y_label,'kz') ||  strcmp(y_label,'Kz')
  y_label_name = 'Kz(1/nm)';
  y_ticks = linspace(-W_z/2+delta_kz,+W_z/2,mum_ticks);
  y_data = (linspace(1,zsizePixels,zsizePixels)-zsizePixels/2.0)*delta_kz;
end

%%

%% data scale
  dataMin = min(min(sectionData));
  dataMax = max(max(sectionData));
  fprintf(num2str(dataMin));
  fprintf('\n');
  fprintf(num2str(dataMax));
  fprintf('\n');
  
  %% use original data to map color, the mapping domain(color_limits)can be adjusted
% scaleSectionData = sectionData;
% color_limits = [dataMin dataMax];

  %% use scale data to map color, the mapping domain(color_limits) can be adjusted,
   %  try [0.5 1] 0r [0,0.5]
  data_Scale = 1.0/(dataMax-dataMin);
  %scaleSectionData = log(sectionData + eps(0));
  scaleSectionData = log(data_Scale*(sectionData-dataMin));
  color_limits = [-5, 0];
  %color_limits = [lower_log, upper_log];
%%

 %% draw image now
 [X,Y] = meshgrid(x_data,y_data);
 surf(X,Y,scaleSectionData,'EdgeColor','none');
 colormap(hot); 
 colorbar
 caxis(color_limits)
 title(title_name);
 grid on;

 ax = gca;
 ax.GridColor = [1, 1, 1];
 ax.Layer = 'top';
 ax.LineWidth = 3;
 ax.GridLineWidth = 2;


 view([0 0 1]);
 xlabel(x_label_name);
 xtickformat('%.3g');
 xticks(x_ticks);
 ylabel(y_label_name);
 ytickformat('%.3g');
 yticks(y_ticks);

 %%
  
end