% ******** Background: ********
%
% The acoustic pressure given by equation (1) in Task 1  describes the
% impulse response of an idealised acoustic system where an acoustic point
% source emits a pulse of ultrasound. This pulse, which is described by
% the ��-distribution, is of infinitesimal duration in time and hence
% exhibits infinite bandwidth (i.e., all frequency components are equally
% represented within the signal). the Gaussian tone burst ,G, of
% Task 2, exhibits a finite bandwidth and only contains frequency
% components across a finite range of frequencies.
%
% ******** Principle: ********
%
% This script performs an over time convolution of the impulse response,
% pressure, computed in Task 1 with the Gaussian tone burst, G, computed
% using comp_Gaussian_tone_burst function.This operation is performed based
% on the equation below:
% s(x,y,z,t) =p(x,y,z,t) *t G(t;f0,sigma),
% where s or s_excited is the the pressure generated by an excitation
% signal of finite temporal duration and *t denotes the convolution
% operator performed over time. This whole operation is equivalent to apply
% a linear filter that suppresses those frequency components not present
% in the pulse shape of finite duration in time.


% ******** Clear Space: ********

% Starts timing
tic;

% Clear all contents from command window
clc;

% Clear all variables from workspace
clear

% Close all figures
close all;


% ******** Computation of p and s: ********

% Import functions and results from Task1_script.
Task1_script;

% sigma is the standard deviation (in seconds) of the Gaussian envelope
% given by exp((-(t.^2)./(exp((-2*t^2)/(2*(sigma^2)))))), where it
% determines the width in time of the tone burst.Unit: [s]

% We define the standard deviation sigma as:
sigma = 5.0000e-08;

% Time temporal step size,dt, is the time increment in every step over a
% finite duration that the pulse can travel, with a unit, [s].

%We define the time temporal step size,dt, as:
dt = 1.0000e-08;

% Define lower time limit of the Gaussian tone burst
lower_t1 = - (4 * sigma);

% Define upper time limit
upper_t1 = (4 * sigma);

% f_0 is the centre frequency which is a measure of a central frequency
% between the upper and lower cutoff frequencies, with a unit, [Hz].

% We define the centre frequency, f_0, as:
f_0 = 1.0000e07;

% Compute the amplitude of the amplitude of the gaussian tone burst
% waveform and time variation range
[ G, time ] = comp_Gaussian_tone_burst(upper_t1,lower_t1,dt,sigma,f_0 );

% Tranform the 4D pressure array into a 2D array,with dimensions
% ((Numx*Numy*Numz, Numt)
pressure_n = reshape(pressure, [], length(t));

% Perform the convolution of the transformed 2D pressure with the Gaussian
% tone burst, returns only the central part of the convolution, the same
% size as pressure_n.s_excited is the pressure generated by an
% excitation signal of finite temporal duration
s_excited = conv( pressure_n, G, 'same');

% Define the time duration of the Gaussian tone burst
t1 = lower_t1:dt:upper_t1;

% Combine the reshaped pressure with the pressure generated by an
% excitation signal of finite temporal duration.
combined_data = [ pressure_n ; s_excited];


% ******** Plot of p and s: ********

% Plot the resulting pressures against time in miliseconds
P = plot(t*1e06, combined_data);

% Change the color of the line of p(x,y,z,t) to blue
set(P(1),'Color','blue');

% Change the color of the line of s(x,y,z,t) to green
set(P(2),'Color','green');

% Label each line with the corresponding name
legend([P(1) P(2)], 'p(x,y,z,t)','s(x,y,z,t)');

% label the x-axis
xlabel('Time [us]');

% label the y-axis
ylabel('Pressure [Pa]');

% giving a title on the figure
title('Pressure as function of time for an acoustic point source');


% ******** Error checking and finish timing: ********

% Find the step that an error occured if there is one
dbstop if error;

% Display the running time of the code and finish timing
fprintf('The total running time of Task3_script is: %.6f seconds. \n',toc');


