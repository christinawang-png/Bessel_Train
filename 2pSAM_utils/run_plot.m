%% A plot drawing Code that takes in a group of saved PSFs from run_PSFs.m
%% Can output the following to a single folder:
%% 1. TIF files of PSFs from all angles
%% 2. FIG and JPG files of middle sections of XY, YZ plots of all angles.
%% 3. FIG and JPG files of direct fourier transforms of middle sections of XY, YZ plots of all angles.
%% 4. FIG and JPG files of summed FFTXY, FFTYZ of all angles.
%% 5. Halfmax Plot of Beam, needs to be manually saved.

clc; clear; close all; addpath(genpath('utils'));
%% Change this according to PSF saved position
PSFs_Path = 'data/PSFs/default'; %Check if PSF exists before generating new PSF
%%%%%%%
PSFs_savePath = ['results//plots_',datestr(now, 'YYYYmmDD_HHMMSS')];%Save results here
mkdir(PSFs_savePath)

%% Change the load path according to PSF saved position
disp('Loading PSFs...');
load([PSFs_Path,'//2p_1.05NA//PSF_mid_zoom6ca128up1psf128_-25-1-25_scanMode2_2-4-8_maxX8.5000Y8.5000biasX0.0Y0.0_pupNO_en0_K0bias0//PSFParameters.mat']);
for i = 1:PSFParameters.angleNum
    load([PSFs_Path,'//2p_1.05NA//PSF_mid_zoom6ca128up1psf128_-25-1-25_scanMode2_2-4-8_maxX8.5000Y8.5000biasX0.0Y0.0_pupNO_en0_K0bias0//psf_all_',num2str(i),'.mat'],'psf_thisAngle');
    psfs(:,:,:,i) = psf_thisAngle;
end
disp('PSFs loaded');
%%%%%%%%

[psf_r, psf_c, psf_s, angleNum]=size(psfs);

xysizePixels = PSFParameters.xysize; %xy point num
delta_xy = PSFParameters.dxy;   %nm, x,y direction FFT sampling period T 
T_xy = xysizePixels * delta_xy; %nm, side length of PSF slice, x,y direction FFT basic period T (periodical function )
zsizePixels = length(PSFParameters.defocus);  %z point num
delta_z  = PSFParameters.defocusInterval*1000;   %nm, z direction FFT sampling period T
T_z = zsizePixels*delta_z; %nm, height og PSF,z direction FFT basic function period T


FFTyz_total = zeros(zsizePixels, xysizePixels); % The preset FFTYZ summed figure
FFTxy_total = zeros(xysizePixels, xysizePixels); % The preset FFTXY summed figure

for N_beam=1:PSFParameters.angleNum
%     %%%%%%%%%%%
% 
%   %Save tif files, can comment out if don't want this
%   stackTempForSave = gather(squeeze(psfs(:,:,:, N_beam)));
%   saveastiff_overwrite(stackTempForSave, [PSFs_savePath, '/Beam',num2str(N_beam),'.tif'],0,1);
%   %%%%%%%%%
% 
%   %compute yz section
  yz = figure('name',['YZ-', num2str(N_beam)]);
% 
%   %% Adjust this base on where the middle section of PSF is, or can customize. EX: a 128*128 size is 64
  X_section = 64;%middle, 1---128
  %%%%%%%%

  imageSctionYZ = (squeeze(psfs(:,X_section,:,N_beam)) + squeeze(psfs(:,X_section+1,:,N_beam)))/2;
  imageSctionYZshow = flip(imageSctionYZ',2);
  displaySectionDataByImage(imageSctionYZshow,PSFParameters,'y','z','YZ Slice X=0', -7, 0);
% 
%   %Save yz section, can comment out if don't want this
%   saveas(yz,[PSFs_savePath,'//yz', num2str(N_beam), '.fig']);
%   saveas(yz,[PSFs_savePath,'//yz', num2str(N_beam), '.jpg']);
%   %%%%%%%%%
  close(yz);
% 
%   %Compute FFTYZ section
  fftyz = figure('name',['FFTYZ-', num2str(N_beam)]);
  FFTimageSctionYZ = fftForImage(imageSctionYZ);
  FFTimageSctionYZshow = flip(FFTimageSctionYZ',2);
  FFTyz_total = FFTyz_total + abs(FFTimageSctionYZshow);
  displaySectionDataByImage(abs(FFTimageSctionYZshow),PSFParameters,'ky','kz','FFT YZ Slice Z=0', -4, 2);
% 
%   %Save fftyz section, can comment out if don't want this
  % saveas(fftyz,[PSFs_savePath,'//fftyz', num2str(N_beam), '.fig']);
  % saveas(fftyz,[PSFs_savePath,'//fftyz', num2str(N_beam), '.jpg']);
%   %%%%%%%%%
  close(fftyz);
% 
% 
%   %Compute xy section
  xy = figure('name',['XY-', num2str(N_beam)]);
% 
%   %% Adjust this base on where the middle Z Layer of PSF is, or can customize. EX: a 51-depth Z would be 25
  z_section = 25;% middle,1--129
  %%%%%%%%

  imageSctionXY = squeeze(psfs(:,:,z_section,N_beam));
  displaySectionDataByImage(imageSctionXY,PSFParameters,'x','y','XY Slice Z=0', -11, -6.5);
% 
%   %Save xy section, can comment out if don't want this
%   saveas(xy,[PSFs_savePath,'//xy', num2str(N_beam), '.fig']);
%   saveas(xy,[PSFs_savePath,'//xy', num2str(N_beam), '.jpg']);
%   %%%%%%%%
  close(xy)
% 
  %Compute FFTXY
  fftxy = figure('name',['FFTXY-', num2str(N_beam)]);
  FFTimageSctionXY = fftForImage(imageSctionXY);
  FFTxy_total = FFTxy_total + abs(FFTimageSctionXY);
  displaySectionDataByImage(abs(FFTimageSctionXY),PSFParameters,'kx','ky','FFT XY Slice Z=0', -4, 0);
% 
%   %Save FFTXY section, can comment out if don't want this
  % saveas(fftxy,[PSFs_savePath,'//fftxy', num2str(N_beam), '.fig']);
  % saveas(fftxy,[PSFs_savePath,'//fftxy', num2str(N_beam), '.jpg']);
%   %%%%%%%%%%
  close(fftxy)
% 
% 
end
fftxy_total = figure('name','FFTXY-total');
displaySectionDataByImage(abs(FFTxy_total),PSFParameters,'kx','ky','FFT XY Total', -4, 0);

%Save Summed FFTXY of all angles
% saveas(fftxy_total,[PSFs_savePath,'//fftxy-total.fig']);
% saveas(fftxy_total,[PSFs_savePath,'//fftxy-total.jpg']);

fftyz_total = figure('name','FFTYZ-total');
displaySectionDataByImage(abs(FFTyz_total),PSFParameters,'ky','kz','FFT YZ Total', -4, 0);

%Save Summed FFTYZ of all angles
% saveas(fftyz_total,[PSFs_savePath,'//fftyz-total.fig']);
% saveas(fftyz_total,[PSFs_savePath,'//fftyz-total.jpg']);


%% Computing Halfmax plot, can comment out below code if don't want it

%% Adjust this according to PSF size to take center of a XY slice. EX: 128*128 xy sizes yields 64 and 64.
X_section = 64;%middle, 1---128
Y_section = 64;
%%%%%%%%%

z_label_name = 'Z(nm)';
z_data_raw = (linspace(1,zsizePixels,zsizePixels)-(zsizePixels+1)/2.0)*delta_z;
z_data = (linspace(1,zsizePixels,zsizePixels*2)-(zsizePixels+1)/2.0)*delta_z;

imageSctionYZ = (squeeze(psfs(Y_section,X_section,:,1)) + squeeze(psfs(Y_section,X_section+1,:,1)))/2;
imageSctionYZshow_raw = flip(imageSctionYZ',2);
imageSctionYZshow = interp1(z_data_raw,imageSctionYZshow_raw,z_data);
pks = max(findpeaks(imageSctionYZshow));
half_pks = pks/2;
const = half_pks*ones(zsizePixels*2,1);

%% Adjust this value manually. Takes all the sampled points that are in the range of:
%% Half_of_Max-smallNumber < points< Half_of_Max + smallNumber
%% Then compute and draw Figure based on Min of Z and Max of Z.
smallNumber = 0.03;

[~,index]=find(abs(imageSctionYZshow-const)<smallNumber);
index_r = unique(index(:).');
len = length(index_r);

inter = z_data(index_r);
bottom = inter(1);
top = inter(len);
half_max = top - bottom;
bottom_intens = imageSctionYZshow(index_r(1));
top_intens = imageSctionYZshow(index_r(len));

p = figure('name','Halfmax_plot');
plot(z_data,imageSctionYZshow, 'o-','MarkerIndices', [index_r(1), index_r(len)]);
hold on;
plot(z_data, const);
xlabel(z_label_name);
title(['HalfMax =', num2str(half_max, 3), 'nm'])
text(bottom-10000, imageSctionYZshow(index_r(1))+0.01,['(', num2str(bottom, 3), ', ' num2str(bottom_intens, 3), ')']);
text(top, imageSctionYZshow(index_r(len))+0.01,['(', num2str(top, 3), ', ' num2str(top_intens, 3), ')']);




