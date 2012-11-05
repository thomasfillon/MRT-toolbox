function [bw_3db, mlbw, enbw] = window_properties(fhandle,winopt)
% Estimate the common characteristics of an analysis window function
% This function estimate the common characteristics of analysis window
% function.
% [bw_3db, mlbw, enbw] = window_properties(fhandle,winopt)
% window_properties rely on the Signal Processign Toolbox function WINDOW
% and share the same syntax minus the parameter N (the length of the analysis
% window)
% fhandle can be any valid window function name
% winopt, are the options for window function fhandle
% Example
% [bw_3db, mlbw, enbw] = window_properties(@hamming);
% [bw_3db, mlbw, enbw] = window_properties(@blackman, 'periodic');
% See also window, MRT
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
%     along with the MRT Toolbox .  If not, see <http://www.gnu.org/licenses/>.
%
%-------------------------------------------------------------------------

L = 2^16; % Length of the frequency analysis
N = 2^8;  % Analysis window length

% Compute window of length N
if nargin > 1,
    win = window(fhandle,N,winopt);
else 
    win = window(fhandle,N);
end
% Get frequency response over L-point
[h,w] = freqz(win,1,L);

%% -3db Bandwidth
val_3dB = 10^(3/20);

h_max = abs(h(1));      % Main lobe value at f=0;
h_3dB = h_max/val_3dB;  % -3dB relative value

I_3db = find(abs(h)<h_3dB,1,'first');     % First estimate of the -3dB location

% Zoom frequency response around -3dB area
[h_3db,w_3db] = freqz(win,1,linspace(w(I_3db-1),w(I_3db),N));
% Estimate -3dB location
I = find(abs(h_3db)<(h_max/val_3dB),1,'first');
% Get -3dB bandwidth
bw_3db = w_3db(I)*N/pi;

%% Main lobe width

diff_h = diff(abs(h(I_3db:end)));
I = find(diff_h>0,1,'first');

mlbw = w(I+I_3db-1)*N/pi; 

%% Equivalent Noise Bandwidth
enbw = (1/N*sum(win.^2))/((1/N*sum(win)).^2);

end