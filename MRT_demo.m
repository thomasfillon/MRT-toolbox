% Demonstration script for the classe MRT
% This script can serve as a rough unitary test
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
% See also : MRT.m, MRT_demo_CQT.m, MRT_demo_MIDI, MRT_demo_ERB.m




clear
close all
% Load Sound file
[x,Fs] = wavread('audio/E-Core - Pingouin-Banquise_45s.wav');
x=sum(x,2); % Conversion from stereo to mono
x=x(1:20*Fs); % Select the first 10 seconds

%% MRT object creation
% Default MRT object creation
mrt_obj = MRT();
display(mrt_obj);

% MRT object creation with Fs
mrt_obj = MRT(Fs);
display(mrt_obj);

% MRT object creation with Fs and the frequency unit
mrt_obj = MRT(Fs, 'Hz');
display(mrt_obj);

mrt_obj = MRT(Fs, 'MIDI');
display(mrt_obj);

mrt_obj = MRT(Fs, 'ERB');
display(mrt_obj);

% Change Fs 
mrt_obj.Fs = 32000; display(mrt_obj);
mrt_obj.Fs = Fs; display(mrt_obj);

% Change the frequency unit 
mrt_obj.freqUnit = 'Hz'; display(mrt_obj);

%% Frequencies specification
% Set the frequency for each bin
freq = 0:1000:5000;
mrt_obj.freqHz = freq;
display(mrt_obj);

mrt_obj.frequencies = freq;
display(mrt_obj);

mrt_obj.freqUnit = 'ERB'; 
mrt_obj.frequencies = 0:1:10;
display(mrt_obj);

%% Frequency resolution specification
mrt_obj.deltafreqHz = 100; % Constant Resolution in Hz
display(mrt_obj);

mrt_obj.delta_frequencies = 1; % Constant Resolution in ERB
display(mrt_obj);

mrt_obj.delta_frequencies = 1:1:11; % Specify Resolution for each bin in ERB
display(mrt_obj);

mrt_obj.N = 500; % Constant kernel length for all bin
display(mrt_obj);

mrt_obj.N = [ones(1,6)*400 ones(1,5)*600]; % Specify every N value
display(mrt_obj);

mrt_obj.N(1:3) = 250; % Specify only few N values
display(mrt_obj);

% Resolution Criteria 
mrt_obj.resolutionCriteria = '3dB';  % Frequency resolution = -3dB  bandwidth
mrt_obj.delta_frequencies = 1; display(mrt_obj);

mrt_obj.resolutionCriteria = 'ML';   % Frequency resolution = Main lobe width
mrt_obj.delta_frequencies = 1; display(mrt_obj);

mrt_obj.resolutionCriteria = 'EN';   % Frequency resolution = Equivalent Noise bandwidth
mrt_obj.delta_frequencies = 1; display(mrt_obj);

% Specify the Q-factor
mrt_obj.freqHz(1) = 1;
mrt_obj.Q = 1;
display(mrt_obj);

% Try to specify Q with a null frequency
mrt_obj.freqHz(1) = 0;
try 
    mrt_obj.Q = 1;
catch exception
    if strcmp(exception.identifier,'MRT:Q:zero_frequency')
        display(exception.message);
    else
      throw(exception);
    end
end
display(mrt_obj);


%% Kernel Parameter
clear mrt_obj;
mrt_obj = MRT(Fs, 'ERB');
mrt_obj.frequencies = 0:1:10;
mrt_obj.delta_frequencies = 1; 

mrt_obj.display_kernels;

mrt_obj.display_kernels('MATRIX'); % Default

mrt_obj.display_kernels('TF'); % XY ticks in time-frequency units


% Specify Time alignement of the kernels
mrt_obj.windowAlignment = 'right';
mrt_obj.display_kernels;

mrt_obj.windowAlignment = 'left';
mrt_obj.display_kernels;

mrt_obj.windowAlignment = 'center';
mrt_obj.display_kernels;

% Specify the window function (cf. matlab 'window' function)
mrt_obj.windowFunction = @blackman;
mrt_obj.display_kernels;


mrt_obj.resolutionCriteria = '3dB';  % Frequency resolution = -3dB  bandwidth
mrt_obj.windowFunction = @gausswin;
mrt_obj.windowOption = 2.5; % Specify window option ('winopt')
mrt_obj.display_kernels; 



%% Transform

close all;
clear mrt_obj;
mrt_obj = MRT(Fs, 'ERB');
mrt_obj.frequencies = 0:1:30;
mrt_obj.delta_frequencies = 1; 
% mrt_obj.R = 256;
mrt_obj.R = floor(mrt_obj.Nmin/2); % Specify the Step-Size
y = mrt_obj.transform(x);
% Display
mrt_obj.display_MRT(y);