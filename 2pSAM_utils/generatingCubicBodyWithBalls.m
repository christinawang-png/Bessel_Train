function [X_body] = generatingCubicBodyWithBalls(num_ball,max_radius)
%% A function to Generate Random Balls given parameters
%  num_ball: number of balls in body 
%  max_radius (um): max radius of balls, every ball radius is generated
%  randomly within (0, max_radius)
%% original 2pSAM testing data example parameter, can be changed
FOV_zoom8_um = 86;% x-y view when zoomFactor = 8: unit (um)
zoomFactor = 4; % zoom ratio
XY_example_dimension = FOV_zoom8_um * 8 / zoomFactor; % 172 (um)
XYsize = 512; % x or y direction pixel number
xyInterval = XY_example_dimension/XYsize; % 172/512 = 0.336 (um)
defocusUpper = 30; % z max size : unit (um)
defocusLower = -30; % z min size : unit (um)
defocusInterval = 0.5; % z interval : unit (um)
Zsize = (defocusUpper - defocusLower) / defocusInterval + 1;% (30-(-30))/0.5+1 = 121 (slice number)
%%

%% generating new 3-d body parameter
xy_pixel_number = XYsize; % 512, may be adjusted
xy_pixel_size = xyInterval * 1000; % 0.336 *1000 = 336(nm),do not modify if comparing with original example
z_slice_number = Zsize; % 121, may be adjusted
z_slice_interval = defocusInterval*1000; % 0.5*1000 = 500(nm),do not modify if comparing with original example
xy_length = xy_pixel_number * xy_pixel_size; 
z_height = (z_slice_number-1) * z_slice_interval;
%%

%% body initial
X_body = zeros(xy_pixel_number,xy_pixel_number,z_slice_number); %initialized as zero
%%

%%
%% generating ball centers ,radius and grey value 
ball_numbers = num_ball; % number of balls ; may be adjusted
ball_radius_max = max_radius*1000;% 5(um) = 5000(nm) max radius; may be adjusted
x_centers = random('Uniform',0,1,ball_numbers)*xy_length; % x coordinates of ball centers
y_centers = random('Uniform',0,1,ball_numbers)*xy_length; % y coordinates of ball centers
z_centers = random('Uniform',0,1,ball_numbers)*z_height;  % z coordinates of ball centers
r_radius = random('Uniform',0,1,ball_numbers)*ball_radius_max; % radius of balls
grey_balls = random('Uniform',0,1,ball_numbers); % grey value of balls
%%
%% computing gray values of every grid point of body
for ii = 1:1:xy_pixel_number
    for jj = 1:1:xy_pixel_number
        for kk = 1:1:z_slice_number
            % coordinates of grid points
            x = (0.5 + (ii-1))*xy_pixel_size;
            y = (0.5 + (jj-1))*xy_pixel_size;
            z = (kk-1)*z_slice_interval;
            [greyValue] = findAveragedGreyValue(x_centers,y_centers,z_centers,r_radius,grey_balls,...
                             ball_numbers,x,y,z);
            X_body(ii,jj,kk) =  greyValue;
        end
    end
end
%%

end

function [greyValue] = findAveragedGreyValue(x_centers,y_centers,z_centers,r_radius,grey_balls,...
                             ball_numbers,x,y,z)
%
%
%
  greyValue = 0; num = 0; % intial grey value
for ii = 1:1:ball_numbers
    % each ball
    xc = x_centers(ii);
    yc = y_centers(ii);
    zc = z_centers(ii);
    r = r_radius(ii);
    greyBall = grey_balls(ii);
    % if or not in ball
    [flag_in_ball] = isInBall(xc,yc,zc,r,x,y,z);
    % find the averaged grey value if ball overlaps
    if flag_in_ball 
      greyValue = greyValue + greyBall; 
      num = num +1;
    end
end
 % average
if num ~= 0
    greyValue = greyValue/num;
end

end

function [flag_in_ball] = isInBall(xc,yc,zc,r,x,y,z)
%
% judge if point in ball
%
  f = (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)-r*r;
if f > 0
   flag_in_ball = false;
else
  flag_in_ball = true;
end

end

