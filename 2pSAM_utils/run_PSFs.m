%% This script shows how to generate point spread functions (PSFs) in 2pSAM
% PSFGeneration
% ELi, 20220726

%% clear and paths
%pinholes = [1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0 6.4 6.8 7.2 7.6 8.0 8.4 8.8]

clc; clear; close all; addpath(genpath('utils'));
%pinholes = [1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0 6.4 6.8 7.2 7.6 8.0 8.4 8.8]
pinholes = [1]; %pinholes: how many times bigger the pinhole is set relative to the original pinhole size.
for p = 1 : length(pinholes)
    pinhole_size = pinholes(p);
    PSFs_savePath = ['results//PSFs_',datestr(now, 'YYYYmmDD_HHMMSS')];%Check if PSF exists before generating new PSF
    systemParameters_path = 'data//systemParameters_simu//25X';

    %% generate PSFs

    %%%see PSFsGenerator.m for more PSF parameters
    PSFParameters.imagingMode = 'mid';% use 'mid' for mid-NA 2pSAM configuration / use 'min' for min-NA 2pSAM configuration
    [psfs,PSFParameters] = PSFsGenerator(pinhole_size, PSFParameters,systemParameters_path,PSFs_savePath);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for N_beam=1:PSFParameters.angleNum
    N_beam = 1;
    %%%%%%%%%%%%
  yz = figure('name',['YZ-', num2str(N_beam)]);
  X_section = 64;%middle, 1---128
  imageSctionYZ = (squeeze(psfs(:,X_section,:,N_beam)) + squeeze(psfs(:,X_section+1,:,N_beam)))/2;
  imageSctionYZshow = flip(imageSctionYZ',2);
  displaySectionDataByImage(imageSctionYZshow,PSFParameters,'y','z','YZ Slice X=0', -7, 0);
  % saveas(yz,[PSFs_savePath,'//yz', num2str(N_beam), '.fig']);
  % saveas(yz,[PSFs_savePath,'//yz', num2str(N_beam), '.jpg']);
  % close(yz);
% 
%   fftyz = figure('name',['FFTYZ-', num2str(N_beam)]);
%   FFTimageSctionYZ = fftForImage(imageSctionYZ);
%   FFTimageSctionYZshow = flip(FFTimageSctionYZ',2);
%   displaySectionDataByImage(abs(FFTimageSctionYZshow),PSFParameters,'ky','kz','FFT YZ Slice Z=0', -4, 2);
%   saveas(fftyz,[PSFs_savePath,'//fftyz', num2str(N_beam), '.fig']);
%   saveas(fftyz,[PSFs_savePath,'//fftyz', num2str(N_beam), '.jpg']);
%   close(fftyz);
% 
%   xy = figure('name',['XY-', num2str(N_beam)]);
%   z_section = 61;% middle,1--129
%   imageSctionXY = squeeze(psfs(:,:,z_section,N_beam));
%   displaySectionDataByImage(imageSctionXY,PSFParameters,'x','y','XY Slice Z=0', -11, -6.5);
%   saveas(xy,[PSFs_savePath,'//xy', num2str(N_beam), '.fig']);
%   saveas(xy,[PSFs_savePath,'//xy', num2str(N_beam), '.jpg']);
%   close(xy)
% 
%   fftxy = figure('name',['FFTXY-', num2str(N_beam)]);
%   FFTimageSctionXY = fftForImage(imageSctionXY);
%   displaySectionDataByImage(abs(FFTimageSctionXY),PSFParameters,'kx','ky','FFT XY Slice Z=0', -4, 0);
%   saveas(fftxy,[PSFs_savePath,'//fftxy', num2str(N_beam), '.fig']);
%   saveas(fftxy,[PSFs_savePath,'//fftxy', num2str(N_beam), '.jpg']);
%   close(fftxy)
% end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%
    %N_beam = 3;
    %%%%%%%%%%%%
    % figure('name','YZ-3');
    % X_section = 64;%middle, 1---128
    % imageSctionYZ = (squeeze(psfs(:,X_section,:,N_beam)) + squeeze(psfs(:,X_section+1,:,N_beam)))/2;
    % imageSctionYZshow = flip(imageSctionYZ',2);
    % imageShow(imageSctionYZshow);
    %
    % figure('name','FFTYZ-3');
    % FFTimageSctionYZ = fftForImage(imageSctionYZ);
    % FFTimageSctionYZshow = flip(FFTimageSctionYZ',2);
    % imageShow(abs(FFTimageSctionYZshow));
    %
    % figure('name','XY-3');
    % z_section = 64;% middle,1--129
    % imageSctionXY = squeeze(psfs(:,:,z_section,N_beam));
    % imageSctionXYshow = flip(imageSctionXY',2);
    % imageShow(imageSctionXYshow);
    %
    % figure('name','FFTXY-3');
    % FFTimageSctionXY = fftForImage(imageSctionXY);
    % FFTimageSctionXYshow = flip(imageSctionXY',2);
    % imageShow(abs(FFTimageSctionXY));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end