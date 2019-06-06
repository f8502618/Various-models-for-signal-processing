function [ Amplitude, t ] = comp_Gaussian_tone_burst( upper_t,lower_t,dt,sigma,f_0 )

% ******** Background: ********
%
% Actual acoustic sources generate acoustic pulses that have a finite
% duration in time and contain frequency components over only a limited
% range of frequencies (i.e. the acoustic energy is mainly contained
% within a finite range of frequencies).To extend the simulation by
% incorporating a more realistic acoustic pulse shape, we consider a
% waveform called a normalised Gaussian tone burst, G, whereas its equation
% is given by:
% G(t;f_0,sigma) = exp(-(t^2)/(2*(sigma^2))) * sin(2*pi*f_0*t).
% Here t is time (in seconds), f_0 is the centre frequency(in Hz)
% ,and sigma is the standard deviation determines the width in time of the
% tone burst.
%
% ******** Function: ********
%
%[ Amplitude, t ]= comp_Gaussian_tone_burst( upper_t,lower_t,dt,sigma,f_0 )
%
% This function calculates the Gaussian tone burst for arbitrary input
% values of centre frequency f_0, standard deviation, sigma, and temporal
% step size dt.
%
% INPUTS:
%
% f_0: the centre frequency, [Hz]
%
% sigma: standard deviation,determines the width in time of the tone burst,
% [s]
%
% lower_t: the lower limit of the time range, [s]
%
% upper_t: the upper limit of the time range, [s]
%
% dt: the temporal step size,[s]
%
% OUTPUT:
%
% Amplitude: the amplitude of the gaussian tone burst waveform, [a.u.]
%
% t: the finite duration of the acoutic pulses travelled past a fixed point
% in time, or, say, the time range of the gaussian tone burst, [s]

% ******** Time range: ********

% As a Gaussian tone burst over the range -4 * sigma to +4 * sigma is
% a good approximation to a Gaussian function. Hence we define the time
% range as:

t = lower_t:dt:upper_t;

% ******** ERROR CHECKING: ********
%
% Check the inputs are real and numeric, if the condition is not satisfied
% , display the error.Check the size of the time range is a (1 row *
% length(t) coloumns) array, if the condition is not satisfied, display the
% error.

if ~isnumeric([upper_t,lower_t,dt,sigma,f_0]) || ~isreal([upper_t,lower_t,dt,sigma,f_0])
    error('Input [upper_t,lower_t,dt,sigma,f_0] is expected to be numeric and real-valued');
end

if length(size(t))~= 2
    error('Input time t is expected to be a 2D array');
end

% Only if all the above test passed successfully, execute remainder.

% ******** Computation: ********

% calculate the denominator part in the exponential formula
exp_de = 2 .*(sigma.^2);

% add eps to avoid exp_de yielding "not a number" (nan):
exp_part1 = exp_de + eps;

% Calculate the exponential part of the Gaussian tone burst,Gaussian
% Envelop
G_envelop = exp((-(t.^2)./(exp_part1)));

% Calculate the sine part of the Gaussian tone burst
sin_part = sin(2.*(pi+eps).* f_0.*t);

% Compute the Gaussian tone burst
Amplitude = G_envelop .* sin_part;

% Compute all coordinates as two arrays:
[t, Amplitude] = meshgrid(t, Amplitude);

% Exact the first coloumn dataset, which is the valid data needed to use,
% others coloumns values in Amplitude are just repeats of the first coloumn
% dataset.
Amplitude = Amplitude(:,1);

% ******** check: ********

% Calculate the sum of amplitudes and round the sum to the order of a.u.
S = roundn((sum(Amplitude).* dt),-20);

% The Gaussian tone burst is supposed to be antisymmetric, hence the
% integral of G*dt for all time value is supposed to equal to
% zero.If the condiction is working fine, display 'The
% comp_Gaussian_tone_burst function is working fine'. Otherwise, display
% the error.

if S == 0
    disp('The comp_Gaussian_tone_burst function is working fine');
end

if S ~= 0
    error(' something went wrong.');
end

end



