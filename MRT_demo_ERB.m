% Example to illustrate the use of the MRT toolobox to perform a 
% Time-Frequency analysis by specifying ERB frequencies
%
% Copyright (c) 2012, Thomas Fillon,
% Institut Mines Telecom, Telecom-ParisTech, CNRS-LTCI
% All rights reserved.
%
%-------------------------------------------------------------------------
%
%     This file is part of the MRT Toolbox.
%
%     The MRT Toolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     The MRT Toolbox is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with the MRT Toolbox.  If not, see <http://www.gnu.org/licenses/>.
%
%-------------------------------------------------------------------------
% See also : MRT.m, MRT_demo.m, MRT_demo_CQT.m, MRT_demo_ERB.m

clear
close all
% Load Sound file
[x,Fs] = wavread('audio/33711__acclivity__excessiveexposure.wav');
x = sum(x,2); % Conversion to mono

%% Specification
% Specify the frequencies to analyze
f_min = 0; % Minimum frequency in Hz
f_max =  16000; % Maximum frequency in Hz

erb_min = MRT.freq2erb(f_min); % Minimum frequency in ERB scale
erb_max = MRT.freq2erb(f_max); % Maximum frequency in ERB scale

nb_bins = 80; % Number of bins to analyze

erb_freq = linspace(erb_min,erb_max,nb_bins);

% Specify the resolution of each bins
ERB_RES = 1; % Frequency resolution in ERB
erb_resolution = ERB_RES * ones(1,nb_bins); % nb_bins bins with a uniform frequency resolution on the erb scale

%% MRT object creation
% Initialize MRT object with the :
% - Fs, the sampling frequency (Fs)
% - 'ERB' the frequency unit
mrt_erb = MRT(Fs,'ERB');

% Set the frequency for each bin in ERB
mrt_erb.frequencies = erb_freq;
% Set the resolution in bins in ERB
mrt_erb.delta_frequencies = erb_resolution;
% Set the step-size
mrt_erb.R = floor(mrt_erb.Nmin/2);
% Display the MRT transform matrix
imagesc(20*log10(abs(mrt_erb.MRT_mat))')
set(gca,'YDir','normal');
xlabel('Time');
ylabel('Frequency');
title('Transform matrix / Time-Frequency kernels')

%% Compute the MRT transform on signal 'x'
y = mrt_erb.transform(x(:,1));
mrt_erb.display_MRT(y);