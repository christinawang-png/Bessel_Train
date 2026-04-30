%% Reconstruct captured projections
% PSFLoading - imageLoading - reconstruction
% ELi, 20220726
% ELi, 20230210, add multi-site DAO options

%% clear and paths
clc; clear; close all; addpath(genpath('utils'));
PSFs_Path = 'data//PSFs//default'; %Check if PSF exists before generating new PSF
EI_path = 'data//sim_data//default//Bessel_0.2NA_0.7_annulus_1';%file path
EI_name = 'testing.tif';%file name
resultSave_path = ['results//recon_',datestr(now, 'YYYYmmDD_HHMMSS')];%Save results here

% loading PSFs 
disp('Loading PSFs...');
load([PSFs_Path,'//Bessel_0.2NA_0.7_annulus//PSF_mid_zoom6ca128up1psf128_-25-1-25_scanMode2_2-4-8_maxX8.5000Y8.5000biasX0.0Y0.0_pupNO_en0_K0bias0//PSFParameters.mat']);
for i = 1:PSFParameters.angleNum
    load([PSFs_Path,'//Bessel_0.2NA_0.7_annulus//PSF_mid_zoom6ca128up1psf128_-25-1-25_scanMode2_2-4-8_maxX8.5000Y8.5000biasX0.0Y0.0_pupNO_en0_K0bias0//psf_all_',num2str(i),'.mat'],'psf_thisAngle');
    psfs(:,:,:,i) = psf_thisAngle;
end
disp('PSFs loaded');

% disp('Loading PSFs...');
% load([PSFs_Path,'//25X//PSF_cx3cr1//PSFParameters.mat']);
% for i = 1:PSFParameters.angleNum
%     load([PSFs_Path,'//25X//PSF_cx3cr1//psf_all_',num2str(i),'.mat'],'psf_thisAngle');
%     psfs(:,:,:,i) = psf_thisAngle;
% end
% disp('PSFs loaded');

%% loading images
disp('Loading raw data...');
projs = single(loadtiff([EI_path,'//',EI_name]));
disp('Raw data loaded');

%% reconstruct
%%%see reconPipeline.m for more options
reconOpts.DAO = 1; % DAO on
[~,shiftMap] = reconPipeline(psfs,projs,PSFParameters,resultSave_path,reconOpts);
phaseRecon(PSFParameters, shiftMap, resultSave_path); % phase reconstruction