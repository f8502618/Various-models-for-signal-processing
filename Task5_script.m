% ******** Introduction: ********
%
% This script performs the computation of the pressure generated by an
% acoustic point source across an extended xy-plane rather than in a
% single point.The comp_press_field_point_source function was reused to
% compute the pressure for all points in space and time.And the
% comp_Gaussian_tone_burst function was used to implement a Gaussian tone
% burst on the pressure.


% ******** Clear Space: ********

% Starts timing
tic;

% Clear all contents from command window
clc;

% Clear all variables from workspace
clear

% Close all figures
close all;


% ******** Parameters: ********

% Pressure field at acoustic pressure source point(x0,y0,z0) = (0,0,0),
% hence define (x0, y0, z0) position in unit of m
x0 = 0;

y0 = 0;

z0 = 0;

x_s = x0;

y_s = y0;

z_s = z0;

% The pressure generated varies across the x-axis and the y-axis, and
% maintains at z = 0 at the z-direction.
% In the spatial grid, -3 mm <= x <= +3 mm,-3 mm <= y <= +3 mm, hence we
% define the x and y position in unit of m.

xcoors = -0.003: 1.0000e-05 : 0.003;

ycoors = -0.003: 1.0000e-05 : 0.003;

zcoors = 0;

x = xcoors;

y = ycoors;

z = zcoors;

% We use sound speed c = 1500 m/s to simulate water or soft biological
% tissue
c = 1500;

% We use the initial pressure amplitude p_0 = 1 Pa.m
p_0 = 1;

% For Gaussian tone burst, sigma is the standard deviation (in seconds) of
%the Gaussian envelope given by exp(-(t^2)/(2*(sigma^2))), where it
% determines the width in time of the tone burst.Unit: [s]

% We define the standard deviation sigma as:
sigma = 5.0000e-08;

% f_0 is the centre frequency which is a measure of a central frequency
% between the upper and lower cutoff frequencies, with a unit, [Hz].

% We define the centre frequency, f_0, as:
f_0 = 1.0000e07;

% Pressure varys as a function of time, and the time range is
% 0<= t <= 3 miliseconds with a temporal step size ?t = 10 ns.
% So define time t in unit of seconds

dt = 1.0000e-08;

upper_t = 3.0000e-06;

lower_t = 0;

t = lower_t:dt:upper_t;

% Gaussian tone burst has a finite duration(-4*sigma<= t <= +4*sigma) in
% time, so we defien the lower time limit of the Gaussian tone burst as
% lower_tG and the uppere time limit of the Gaussian tone burst as
% upper_tG:

lower_tG = - (4 .* sigma);

upper_tG = (4 .* sigma);


% ******** Methods: ********

% Define the spatial steps size Numx,Numy,Numz, along x,y and z-dimension,
% the temporal step size, Numt, is also defined.
% Compute the pressure as a function of time over a three-dimensional (3D)
% grid of sample points, with a unit of [Pa]. t_end is the time the user
% wishes to terminate and inx finds the index of the time.

Numx = length(xcoors);

Numy = length(ycoors);

Numz = length(zcoors);

Numt = length(t);

t_end = 1.0000e-06;

inx = find(t == t_end);

[pressure,t]  = comp_press_field_point_source(x,x_s,y,y_s,z,z_s,p_0,c,upper_t,lower_t,dt,t_end,inx);


% ******** Plot the impulse response in xy plane at t = 1 us: ********
% Divide the figure into an 2-by-2 grid and creates axes in the position
% specified by 1.
subplot(2,2,1);

% Plot the presseure with a time up to inx, and x is (-3:10e-3:3)mm and
% y is (-3:10e-3:3)mm
imagesc(-3:10e-3:3,-3:10e-3:3,pressure(:,:,1,inx));

% Adjust the xyz ratio
daspect([Numx,Numy,Numz]);

% Choose the colormap mode
colormap hot;

% Set up the colorbar and tells you the pressure value
colorbar;

% Constrain the colorbar to a certain range
caxis([-60,60]);

% label the x-axis
xlabel('x [mm]');

% label the y-axis
ylabel('y [mm]');

% giving a title on the figure
title('Impulse response in xy plane, as t = 1 miliseconds.');


% ******** Plot the field generated by Gaussian tone burst at t = 1 us: ***

% Compute the amplitude of the Gaussian tone burst
[ G, time ] = comp_Gaussian_tone_burst(upper_tG,lower_tG,dt,sigma,f_0 );

% Transform the G into a 4D array
NumG = length(G);

G4d = reshape(G,[1,1,1,NumG]);

% Perform the convolution of pressure and G4d over time, central part of
% the convolution will be same size as pressure.
s_conver = convn(pressure,G4d,'same');

% Divide the figure into an 2-by-2 grid and creates axes in the position
% specified by 2.
subplot(2,2,2);

% Plot the the response to a Gaussian tone burst with a time up to inx,
% and x is (-3:10e-3:3)mm and y is (-3:10e-3:3)mm
imagesc(-3:10e-3:3,-3:10e-3:3,s_conver(:,:,1,inx));

% Adjust the xyz ratio
daspect([Numx,Numy,Numz]);

% Set up the colormap
colormap hot;

% colour bar that explains which colour corresponds to which Gaussian
% response value.
colorbar;

% Limit the colour bar to a certain range
caxis([-60,60]);

% label the x-axis
xlabel('x [mm]');

% label the y-axis
ylabel('y [mm]');

% giving a title on the figure
title('Field generated by Gaussian tone burst.');


% ******** Error checking and finish timing: ********

% Stop the code and indicate the error
dbstop if error;

% Finish timing and return the total running time
fprintf('The total running time of Task5_script is: %.6f seconds. \n',toc');
