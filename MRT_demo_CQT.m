% Example to illustrate the use of the MRT toolobox to perform a CQTanalysis 
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
% See also : MRT.m, MRT_demo.m, MRT_demo_MIDI.m, MRT_demo_ERB.m

clear
close all
% Load Sound file
% Example #1
%[x,Fs] = wavread('audio/E-Core - Pingouin-Banquise_45s.wav');
% Example #2
[x,Fs] = wavread('audio/KIMIKO ISHIZAKA - Goldberg Variations BWV 988 - 01 - Aria_45s.wav');

x=sum(x,2); % Conversion from stereo to mono
x=x(1:10*Fs); % Select the first 10 seconds
%% CQT Specification
f_min = 55; % Minimum frequency in Hz
nb_octave = 7; % Number of octave
B = 24; % Number of bins per octave


%% MRT Specification (derived from CQT specification) 
cqt_freq = f_min*2.^((0:nb_octave*B-1)/B); % Frequency set
Q_factor = 1/(2.^(1/B)-1);                 % Q-factor to specified the frequency resolution
% Note : the Q-factor specified here is not the real Q-factor unless the
% resolution criteria is '3dB'

%% MRT object creation
mrt_cqt = MRT(Fs,'Hz');

% Set the frequency for each bin in Hz
mrt_cqt.frequencies = cqt_freq;
% Set the resolution for each bin through the Q-factor
mrt_cqt.Q = Q_factor;

% Set the step-size
mrt_cqt.R = floor(mrt_cqt.Nmin/2);

%% Display the MRT transform matrix
% With the Time-Frequency units
mrt_cqt.display_kernels('TF');
% Without the Time-Frequency units
%mrt_cqt.display_kernels('MATRIX'); 

%% Compute the MRT transform on signal 'x'
y = mrt_cqt.transform(x);
% Display
mrt_cqt.display_MRT(y);


%% Mixed CQT
% The MRT enable to set a arbitrary frequency resolution for each bin.
% This can be use to enhanced the time resolution in the low frequency bin
% of the CQT and have a hybrid resolution (linear in the low frequency, constant Q in the medium and high frequencies).

% Copy the MRT CQT
mrt_cqt_hybrid = mrt_cqt;
% Keep the same resolution for the high part of the spectrum
% and change the first 40 bins frequency resolution to a uniform resolution
mrt_cqt_hybrid.deltafreqHz(1:40) = mrt_cqt_hybrid.deltafreqHz(40);
% Compute the hybrid CQT transform on signal 'x'
y_hybrid = mrt_cqt_hybrid.transform(x);
% Display
mrt_cqt_hybrid.display_MRT(y_hybrid);
