%% Code to Generate Projections given the 3D olume and PSFs.

clc; clear; close all; addpath(genpath('utils'));
PSFs_Path = 'data//PSFs//400_xy_pixels_dz_1'; %Check if PSF exists before generating new PSF
EI_path = 'data//images';%file path
%data\sim_results\2pSAM_original_results\recon_fRL_zoom4ca512up1_-30-0.5-30_pupNO_en0_K0bias0_DAO1_EW0.5_mI30_mC0_frames1-13-13
EI_name = '3d_400um_to_500um_cropped_xp4um.tif';
%EI_name = '3d_400um_to_500um_cropped_xp4um.tif';%file name
%EI_name = 'balls.tif';%file name
resultSave_path = ['results//generate_',datestr(now, 'YYYYmmDD_HHMMSS')];%Save results here

%data\sim_results\2pSAM_original_results\recon_fRL_zoom4ca512up1_-30-0.5-30_pupNO_en0_K0bias0_DAO1_EW0.5_mI30_mC0_frames1-13-13

disp('Loading PSFs...');
load([PSFs_Path,'//2pSAM//PSF_mid_zoom6ca128up1psf128_-25-1-25_scanMode2_2-4-8_maxX4.2426Y6.0000biasX0.0Y0.0_pupNO_en0_K0bias0//PSFParameters.mat']);
for i = 1:PSFParameters.angleNum
    load([PSFs_Path,'//2pSAM//PSF_mid_zoom6ca128up1psf128_-25-1-25_scanMode2_2-4-8_maxX4.2426Y6.0000biasX0.0Y0.0_pupNO_en0_K0bias0////psf_all_',num2str(i),'.mat'],'psf_thisAngle');
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

[psf_r, psf_c, psf_s, angleNum]=size(psfs);

%% loading images
disp('Loading raw data...');
projs = single(loadtiff([EI_path,'//',EI_name]));
[projr_r,projr_c,projr_num]=size(projs);


%% Clip or Pad to adjust image slice numbers
if mod(projr_num, 2) == 0
    projs = projs(:,:,1:projr_num - 1);
    projr_num = projr_num - 1;
end
if projr_num > psf_s
    difference = projr_num - psf_s;
    clip = difference/2;
    projs = projs(:,:,clip+1:projr_num-clip);
end
if projr_num < psf_s
    difference = psf_s - projr_num;
    pad = difference/2;
    pad_zeros = zeros(projr_r, projr_c, pad);
    projs = cat(3, pad_zeros, cat(3, projs, pad_zeros));
end
disp('Raw data loaded');


%% Normalize PSF and data before calculation if norm flag set 1
norm = 0;
if norm
    disp('Normalizing psfs and data...')
    for psfCount=1:angleNum
        cur_psf = psfs(:, :, :, psfCount);
        psf_maximum = max(cur_psf,[],"all");
        disp(['psf maximum ', num2str(psf_maximum)])
        psfs(:, :, :, psfCount) = psfs(:, :, :, psfCount) ./ psf_maximum;
    end
    projs_maximum = max(projs,[],"all");
    disp(['data maximum ', num2str(projs_maximum)])
    projs = projs ./ projs_maximum;
end


[proj_r,proj_c,proj_num]=size(projs);

Proj_init = zeros(proj_r,proj_c,angleNum);

for angleNow = 1:angleNum
    %largepsf((proj_r-psf_r)/2+1:(proj_r+psf_r)/2,(proj_c-psf_c)/2+1:(proj_c+psf_c)/2,:)=squeeze(psfs(:,:,:,angleNow));
    ProjNow = forwardProj_RL_GPU(squeeze(psfs(:,:,:,angleNow)),projs);
    Proj_init(:,:, angleNow) = ProjNow;
end
stackTempForSave = gather(Proj_init);
saveastiff_overwrite(stackTempForSave,[resultSave_path '//testing.tif'],0,1);

saveastiff_overwrite(projs,[resultSave_path '//data.tif'],0,1);