%% This script shows how to generate point spread functions (PSFs) in 2pSAM
% PSFGeneration
% ELi, 20220726

%% clear and paths
clc; clear; close all; addpath(genpath('utils'));
PSFs_savePath = ['results//PSFs_',datestr(now, 'YYYYmmDD_HHMMSS')];%Check if PSF exists before generating new PSF
systemParameters_path = 'data//systemParameters_simu//25X';

%% generate PSFs
%%%see PSFsGenerator.m for more PSF parameters
PSFParameters.imagingMode = 'mid';% use 'mid' for mid-NA 2pSAM configuration / use 'min' for min-NA 2pSAM configuration 
[psfs,PSFParameters] = PSFsGenerator(PSFParameters,systemParameters_path,PSFs_savePath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
N_beam = 1;
%%%%%%%%%%%%
  figure('name','YZ-1');
  X_section = 64;%middle, 1---128
  imageSctionYZ = (squeeze(psfs(:,X_section,:,N_beam)) + squeeze(psfs(:,X_section+1,:,N_beam)))/2;
  imageSctionYZshow = flip(imageSctionYZ',2);
  displaySectionDataByImage(imageSctionYZshow,PSFParameters,'y','z','YZ Slice X=0');
  
  figure('name','FFTYZ-1');
  FFTimageSctionYZ = fftForImage(imageSctionYZ);
  FFTimageSctionYZshow = flip(FFTimageSctionYZ',2);
  displaySectionDataByImage(abs(FFTimageSctionYZshow),PSFParameters,'ky','kz','FFT YZ Slice Z=0');
 
  figure('name','XY-1');
  z_section = 65;% middle,1--129
  imageSctionXY = squeeze(psfs(:,:,z_section,N_beam));
  displaySectionDataByImage(imageSctionXY,PSFParameters,'x','y','XY Slice Z=0');
  
  figure('name','FFTXY-1');
  FFTimageSctionXY = fftForImage(imageSctionXY);
  displaySectionDataByImage(abs(FFTimageSctionXY),PSFParameters,'kx','ky','FFT XY Slice Z=0');
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 