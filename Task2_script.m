% ******** Introduction: ********
%
% Actual acoustic sources generate acoustic pulses that have a finite
% duration in time and contain frequency components over only a
% limited range of frequencies (i.e. the acoustic energy is mainly
% contained within a finite range of frequencies).
% Gaussian tone burst allows the computor to perform simulations more close
% to a realistic acoutic pulse.A Gaussian over the range -4sigma to +4sigma
% is a good approximation to a Gaussian function, which is used in defining
% the time range. Thea normalised Gaussian tone burst we consisering, G,
% is given by the below equation:
% G(t;f0,sigma) = exp((-2*t^2)/(2*(sigma^2)))*sin(2*pi*f0*t);
% This script uses comp_Gaussian_tone_burst function to compute a Gaussian
% tone burst for the following parameters: f_0 = 10 MHz, sigma = 50 ns, and
% dt = 10 ns.


% ******** Method: ********

% Start timing
tic;

% Clear all contents from command window
clc;

% Clear all variables from workspace
clear

% Close all figures
close all;

% sigma is the standard deviation (in seconds) of the Gaussian envelope
% given by exp((-(t.^2)./(exp((-2*t^2)/(2*(sigma^2)))))), where it determines
% the width in time of the tone burst.Unit: [s]

% We define the standard deviation sigma as:
sigma = 5.0000e-08;

% Time temporal step size,dt, is the time increment in every step over a
% finite duration that the pulse can travel, with a unit, [s].

% We define the time temporal step size,dt, as:
dt = 1.0000e-08;

% Define the lower time limit as:
lower_t = - (4 .* sigma);

% Define the upper time limit as:
upper_t = (4 .* sigma);

% f_0 is the centre frequency which is a measure of a central frequency
% between the upper and lower cutoff frequencies, with a unit, [Hz].

% We define the centre frequency, f_0, as:
f_0 = 1.0000e07;


% ******** Methods: ********

% Compute the amplitude of the gaussian tone burst waveform and time
% variation

[ Amplitude, t ] = comp_Gaussian_tone_burst( upper_t,lower_t,dt,sigma,f_0 );


% ******** Results: ********

% Plot thhe Gaussian tone burst in the time region -4*sigma<= t <= 4*sigma;
% T is the transpose of time,t. The relation betweeen the amplitude and the
% time can only be plotted when  two arrays assigning the values properly.

% Transpose t to make the data align with the Amplitude dataset. Plot the
% Amplitude of the Gaussian tone burst against the time in milisecond.
T =t.';

plot(T*1e06, Amplitude);

% label the x-axis
xlabel('Time [us]');

% label the y-axis
ylabel('Amplitude [a.u.]');

% giving a title on the figure
title('Gaussian tone burst');


% ******** Error checking and finish timing: ********

% Stop the code and indicate the error
dbstop if error;

% Display the time taken to run the code,and finish timing
fprintf('The total running time of Task2_script is: %.6f seconds. \n',toc');