%% The main function to Generate Random Balls as a testing 3D dataset and save them as TIF files.

clc; clear; close all; addpath(genpath('utils'));
resultSave_path = ['results//generate_balls_',datestr(now, 'YYYYmmDD_HHMMSS')];%Save results here


num_ball = 30;% number of balls in body
max_radius = 10;% unit(um),reference to the size of body: 172 X 172 X 60;
[X_body] = generatingCubicBodyWithBalls(num_ball,max_radius);

%% draw  image of given slice 
[x_mum,y_num,z_num] = size(X_body);

disp('Saving circles as a .mat and .tif file...');
stackTempForSave = gather(X_body);
saveastiff_overwrite(stackTempForSave,[resultSave_path '//circles.tif'],0,1);
save([resultSave_path,'/','circles.mat'],'stackTempForSave','-v7.3','-nocompression');
disp('Done!');

% for slice = 1:5:z_num
%   figure;
%   Ximage = X_body(:,:,slice);
%   imshow(Ximage);
%   title(num2str(slice));
% end
%%