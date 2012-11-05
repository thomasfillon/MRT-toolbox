% Example to illustrate the use of the MRT toolobox to perform a CQT
% analysis by specifying MIDI frequencies
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
%% Load Sound file
% Example #1
[x,Fs] = wavread('audio/E-Core - Pingouin-Banquise_45s.wav');
% Example #2
%[x,Fs] = wavread('audio/KIMIKO ISHIZAKA - Goldberg Variations BWV 988 - 01 - Aria.wav');

x=sum(x,2);    % Conversion from stereo to mono
x=x(1:10*Fs);  % Select the first 10 seconds
%% CQT Specification
f_min = 55; % Minimum frequency in Hz
midi_min = MRT.freq2midi(f_min);
nb_octave = 7; % Number of octave
B = 24; % Number of bins per octave


%% MRT Specification (derived from CQT specification) 
midi_freq = midi_min +(0:nb_octave*B-1)/B*12; % Frequency set
midi_resolution = 12/B;                     % Constant frequency resolution on the MIDI scale

%% MRT object creation
mrt_cqt_midi = MRT(Fs,'MIDI');

% Set the frequency for each bin in MIDI
mrt_cqt_midi.frequencies = midi_freq;
% Set the resolution for each bin in MIDI
mrt_cqt_midi.delta_frequencies = midi_resolution;

% Set the step-size
mrt_cqt_midi.R = floor(mrt_cqt_midi.Nmin/2);

%% Display the MRT transform matrix
% With the Time-Frequency units
%mrt_cqt_midi.display_kernels('TF');
% Without the Time-Frequency units
mrt_cqt_midi.display_kernels('MATRIX'); 

%% Compute the MRT transform on signal 'x'
y = mrt_cqt_midi.transform(x);
% Display
mrt_cqt_midi.display_MRT(y);