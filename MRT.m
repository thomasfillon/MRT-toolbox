classdef MRT
    %MRT Matlab Class for the Multi-Resolution time-frequency tranform
    %
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
    % See also : , MRT_demo.m, MRT_demo_CQT.m, MRT_demo_ERB.m, MRT_demo_MIDI.m
    %
    % References :
    %
    % [1] Thomas Fillon and Jacques Prado, "A Flexible Multi-Resolution
    % Time-Frequency Analysis Framework for Audio Signals", in proceedings
    % of the 11th International Conference On Information Science, Signal
    % Processing and their Applications (ISSPA) 3-5 july 2012, Montreal
    % (Canada)
    
    
    %% PROPERTIES
    properties
        % Frequency specification
        freqHz     = [];             % Frequencies in Hertz
        freqUnit   = 'Hz';           % Frequency unit (Hz, Midi or ERB)
        Fs         = 1;              % Sampling Frequency
        % Kernel and Frequency Resolution specification
        N          = [] ;            % Window lengths in samples
        resolutionCriteria  = 'ML'   % Criteria use for the computation of N and delta_frequencies
                                     % - '3dB' for -3dB bandwidth,
                                     % - 'ML' for the main lobe width (default),
                                     % - 'EN' for the Equivalent Noise bandwidth)
        windowFunction = @hamming;       % Analysis Window (default = @hamming)
        windowOption = [];
        windowAlignment = 'center';  % Alignment of the Time-Frequency kernels in time
        
        
        % Transform Specification
        R          = [];             % Stepsize of the MRT transform
        
    end
    properties(Dependent, Hidden)
        bandwidth              % Bandwith of the analysis window depending on the resolution criteria
        tau                    % Time delay of the kernel (to enable proper alignment)
        MRT_mat                % Transform Matrix (dimension = [M x K])
        %iMRT_mat              % Inverse Transform Matrix
    end
    properties(Dependent)
        frequencies            % Center frequencies     (in the unit given by freqUnit)
        delta_frequencies      % Frequencies resolution (in the unit given by freqUnit)
        deltafreqHz    =[];    % Frequency Resolutions in Hertz
        Q                      % Q factors
        
    end
    properties (Dependent, SetAccess = private)
        
        M                      % Frame length
        K                      % Number of bins of the transform
        Nmax
        Nmin
    end
    
    %% METHODS
    methods
        function mrt_obj = MRT(Fs,freqUnit)
            if nargin >= 1
                mrt_obj.Fs = Fs;
            end
            if nargin == 2
                mrt_obj.freqUnit = freqUnit;
            end
        end
        
        function mrt_obj = set.freqUnit(mrt_obj, unit)
            % Set the frequency Unit.
            % The input value should be a string among the following list :
            % - 'Hertz' or 'Hz'
            % - 'MIDI'
            % - 'ERB'
            
            switch unit
                case {'MIDI', 'midi', 'Midi'}
                    mrt_obj.freqUnit = 'MIDI';
                case {'ERB', 'erb'}
                    mrt_obj.freqUnit = 'ERB';
                case {'HERTZ', 'Hertz', 'hertz', 'Hz'}
                    mrt_obj.freqUnit = 'Hz';
                otherwise
                    error('Wrong Frequency Unit. It must be Hertz, MIDI or ERB');
            end
        end
        function mrt_obj = set.resolutionCriteria(mrt_obj,Criteria)
            % Criteria to take into account for the frequency resolution (could be :
            % - '3dB' for -3dB bandwidth,
            % - 'ML' for the main lobe width (default),
            % - 'EN' for the Equivalent Noise bandwidth)
            
            if any(strcmp(Criteria,{'3dB','ML','EN'}))
                mrt_obj.resolutionCriteria = Criteria;
            else
                error('Wrong input for the "resolutionCriteria" value');
            end
        end
        
        function freq = get.frequencies(mrt_obj)
            % Get the frequency values in the unit given by 'freqUnit'
            switch mrt_obj.freqUnit
                case 'MIDI'
                    freq = MRT.freq2midi(mrt_obj.freqHz);
                case 'ERB'
                    freq = MRT.freq2erb(mrt_obj.freqHz);
                case 'Hz'
                    freq = mrt_obj.freqHz;
            end
        end
        
        function deltaFreq = get.delta_frequencies(mrt_obj)
            if ~isempty(mrt_obj.deltafreqHz)
                freqHz_low  = mrt_obj.freqHz - mrt_obj.deltafreqHz/2;
                freqHz_high = mrt_obj.freqHz + mrt_obj.deltafreqHz/2;
                
                switch mrt_obj.freqUnit
                    case 'MIDI'
                        deltaFreq = MRT.freq2midi(freqHz_high) - MRT.freq2midi(freqHz_low);
                    case 'ERB'
                        deltaFreq = MRT.freq2erb(freqHz_high) - MRT.freq2erb(freqHz_low);
                    case 'Hz'
                        deltaFreq = mrt_obj.deltafreqHz;
                end
            else deltaFreq = [];
            end
            
        end
        
        function deltaFreqHz = get.deltafreqHz(mrt_obj)
            switch mrt_obj.resolutionCriteria
                case '3dB'
                    alpha = mrt_obj.bandwidth.m3dB;
                case 'ML'
                    alpha = mrt_obj.bandwidth.ML;
                case 'EN'
                    alpha = mrt_obj.bandwidth.EN;
            end
            
            deltaFreqHz = alpha./mrt_obj.N*mrt_obj.Fs;
        end
        
        function mrt_obj = set.N(mrt_obj,N_values)
            % Set N the analysis window length for each bin
            % and consequently update the frequency resolution values
            
            mrt_obj.N = N_values(:)';
        end
        
        function Q = get.Q(mrt_obj)
            if ~isempty(mrt_obj.deltafreqHz)
                % Get the Q-factor
                % This Q-factor is the real Q-factor only if the
                % resolution criteria is '3dB' otherwise deltafreqHz is not
                % the -3dB bandwith
                
                Q = mrt_obj.freqHz./mrt_obj.deltafreqHz;
                if ~strcmp(mrt_obj.resolutionCriteria,'3dB')
                    
                    warning('MRT:Qvalue',...
                        'The Q-factor given here is not the proper Q-factor because the resolution criteria is not ''3dB''');
                end
            else Q = [];
            end
        end
        
        function mrt_obj = set.Q(mrt_obj,Q_values)
            % Set Q the Q-factor for each bin
            % and consequently update the frequency resolution values
            % This Q-factor is the real Q-factor only if the
            % resolution criteria is '3dB' otherwise deltafreqHz is not
            % the -3dB bandwith
            
            % Update the frequency resolution values
            mrt_obj.deltafreqHz = mrt_obj.freqHz./Q_values;
            if any(mrt_obj.freqHz==0)
                error('MRT:Q:zero_frequency',...
                    'The Q-factor cannot be specidied if the null frequency is in the frequency list');
            end
            if ~strcmp(mrt_obj.resolutionCriteria,'3dB')
                
                warning('MRT:Qvalue',...
                    'The Q-factor given here is not the proper Q-factor because the resolution criteria is not ''3dB''');
            end
        end
        
        function Nmax = get.Nmax(obj)
            Nmax = max(obj.N);
        end
        
        function Nmin = get.Nmin(obj)
            Nmin = min(obj.N);
        end
        % Frame length
        function M = get.M(mrt_obj)
            M = mrt_obj.Nmax;
        end
        % K : Number of bins
        function K = get.K(mrt_obj)
            if not(isempty(mrt_obj.freqHz) | isempty(mrt_obj.N))
                K = length(mrt_obj.freqHz);
            else
                K = 0;
            end
        end
        
        function tau = get.tau(mrt_obj)
            switch mrt_obj.windowAlignment
                case 'center'
                    tau = round((mrt_obj.M-mrt_obj.N)/2);
                case 'right'
                    tau = mrt_obj.M-mrt_obj.N;
                    
                case 'left'
                    tau = zeros(1,mrt_obj.K);
            end
        end
        
        
        function mrt_obj = set.frequencies(mrt_obj, freqValues)
            % Set the frequencies set for the Multi Resolution Transform
            % The frequency could be either specified in Hertz, MIDI or ERB
            % as indicated in freqUnit
            
            switch mrt_obj.freqUnit
                case 'MIDI'
                    mrt_obj.freqHz = MRT.midi2freq(freqValues);
                case 'ERB'
                    mrt_obj.freqHz = MRT.erb2freq(freqValues);
                case 'Hz'
                    mrt_obj.freqHz = freqValues;
            end
        end
        
        function mrt_obj = set.delta_frequencies(mrt_obj,resolutionValues)
            % Set the frequency resolutions set for the Multi Resolution Transform
            % The frequency resolutions could be either specified in Hertz, MIDI or ERB
            % as indicated in freqUnit
            
            switch mrt_obj.freqUnit
                case 'MIDI'
                    deltaFreqHz = MRT.midi_resolution(mrt_obj.freqHz,resolutionValues);
                case 'ERB'
                    deltaFreqHz = MRT.ERB_resolution(mrt_obj.freqHz,resolutionValues);
                case 'Hz'
                    if isscalar(resolutionValues),
                        resolutionValues = resolutionValues * ones(1,length(mrt_obj.freqHz));
                    end
                    deltaFreqHz = resolutionValues;
            end
            
            mrt_obj.deltafreqHz = deltaFreqHz;
        end
        
        function mrt_obj = set.deltafreqHz(mrt_obj,deltaFreq)
            switch mrt_obj.resolutionCriteria
                case '3dB'
                    alpha = mrt_obj.bandwidth.m3dB;
                case 'ML'
                    alpha = mrt_obj.bandwidth.ML;
                case 'EN'
                    alpha = mrt_obj.bandwidth.EN;
            end
            mrt_obj.N = round(alpha./deltaFreq*mrt_obj.Fs);
        end
        function mrt_obj = set.windowAlignment(mrt_obj,alignment)
        %% Define Alignment
            switch alignment
                case {'center', 'Center'}
                    mrt_obj.windowAlignment = 'center';
                case {'right', 'Right'}
                    mrt_obj.windowAlignment = 'right';
                case {'left', 'Left'}
                    mrt_obj.windowAlignment = 'left';
                otherwise
                    error('Wrong input for windowlAlignment. It must be ''center''(default), ''right'' or ''left''')
            end
        end
        function MRT_mat = get.MRT_mat(mrt_obj)
            % Set the transform matrix of a MRT object
            
            MRT_mat = zeros(mrt_obj.M,mrt_obj.K);
            
            for k=1:mrt_obj.K,
                if isempty(mrt_obj.windowOption)
                    w = window(mrt_obj.windowFunction,mrt_obj.N(k));
                else
                    w = window(mrt_obj.windowFunction,mrt_obj.N(k),mrt_obj.windowOption);
                end
                
                gamma_k = sum(w); % Normalization factor
                
                % zero before kernel
                %if mrt_obj.tau(k)>1,
                %    MRT_mat(1:mrt_obj.tau(k),k) = 0;
                %end;
                % Time-Frequency kernel
                MRT_mat(mrt_obj.tau(k)+(1:mrt_obj.N(k)),k) = 1/gamma_k*w.*exp(1i*2*pi*(0:mrt_obj.N(k)-1).'*mrt_obj.freqHz(k)/mrt_obj.Fs);
                % zeros after kernel
                %if 1+mrt_obj.N(k)+mrt_obj.tau(k) < mrt_obj.M,
                %    MRT_mat(1+mrt_obj.N(k)+tau(k):mrt_obj.M,k) = 0 ;
                %end;
            end
        end
        
        function bandwidth = get.bandwidth(mrt_obj)
            %% Bandwidth estimation upon Resolution criteria
            if isempty(mrt_obj.windowOption)
                [bw_3dB, ml_bw, enbw] = window_properties(mrt_obj.windowFunction);
            else
                [bw_3dB, ml_bw, enbw] = window_properties(mrt_obj.windowFunction,mrt_obj.windowOption);
            end
            bandwidth = struct('m3dB', bw_3dB, 'ML', ml_bw, 'EN', enbw);
        end
    end
    methods(Static)

        % Frequency Conversion functions
        function f_midi = freq2midi(f_Hz)
            f_midi = (12*log2(f_Hz/440))+69;
        end
        function f_Hz = midi2freq(f_midi)
            f_Hz = 440 * exp ((f_midi-69) * log(2)/12);
        end
        function f_erb = freq2erb(f_Hz)
            f_erb = sign(f_Hz).* (21.4*log10(abs(f_Hz)*4.37/1000+1));
        end
        function f_Hz = erb2freq(f_erb)
            f_Hz  = (10.^(f_erb./21.4) - 1) * 1000/4.37;
        end
        function ERB_width = ERBwidth(f_Hz)
            ERB_width	= 24.7.*(4.37*f_Hz/1000 + 1);
        end
        function delta_fHz = ERB_resolution(fc_Hz, delta_ERB)
            alpha = 10.^(delta_ERB/21.4);
            delta_fHz = 2*(alpha-1).*(1+fc_Hz*4.37/1000)./((1+alpha)*4.37/1000);
        end
        function delta_fHz = midi_resolution(fc_Hz, delta_midi)
            alpha = 2.^(delta_midi/12);
            delta_fHz = 2*fc_Hz.*(alpha-1)./(alpha+1);
        end
        
        
    end
    methods(Hidden)
        function refresh_Yticks(mrt_obj,axes_handle)
            % Provide Y-ticks update with right frequency values when zooming
            set(axes_handle,'YTickMode','Auto');
            Ytick = get(axes_handle,'Ytick');
            Ytick = unique(round(Ytick));
            set(axes_handle,'Ytick',Ytick);
            set(axes_handle,'YTickLabel',mrt_obj.frequencies(Ytick));
        end
        function display_zoom_ticks(mrt_obj,~,evd)
            % Call back used for zooming inside a MRT display
            mrt_obj.refresh_Yticks(evd.Axes);
        end
    end
    %% Transform
    methods
        function  X_MRT = transform(mrt_obj,x)
            
            MAX_BLOCK = 2^20; % Block size for piecewise processing of a large signal
            
            % Check x dimension
            if ~isvector(x)
                error('Input data must be a vector')
            end
            x = x(:);
            
            % Check the existence of the stepsize R
            if isempty(mrt_obj.R)
                error('MRT:Transform:R','The step-size value ''R'' is not defined');
            end
            
            % Round MAX_BLOCK upon the stepsize R
            block_length = MAX_BLOCK - rem(MAX_BLOCK-(mrt_obj.M-mrt_obj.R),mrt_obj.R);
            
            
            %% Time implementation
            zeropad = zeros(ceil(mrt_obj.M/2)-1,1);
            x = [zeropad; x ; zeropad];
            
            % number of MRT frames for the whole signal
            nb_frames = ceil((length(x)-(mrt_obj.M-mrt_obj.R))/mrt_obj.R);
            X_MRT = zeros(mrt_obj.K,nb_frames);
            
            
            block_overlap = mrt_obj.M-mrt_obj.R;
            block_stepsize = block_length - block_overlap;
            nb_blocks =  ceil((length(x)-block_overlap)/block_stepsize);
            
            nb_MRT_blocks =  ceil((block_length-block_overlap)/mrt_obj.R);
            
            MRT_mat = mrt_obj.MRT_mat';
            
            for i_block =1:nb_blocks,
                i_start = 1 + (i_block-1)*block_stepsize;
                i_end = min(i_start-1+block_length, length(x));
                
                i_MRT_start = 1 + (i_block-1)*nb_MRT_blocks;
                i_MRT_end = min(i_MRT_start-1+nb_MRT_blocks,nb_frames);
                X_MRT(:,i_MRT_start:i_MRT_end) = mrt_transform(x(i_start:i_end));
            end
            
            
            % Nested Transform function
            function X_MRT_block = mrt_transform(x_block)
                % Time overlapping buffer
                X_block = buffer(x_block,mrt_obj.M,mrt_obj.M-mrt_obj.R,'nodelay');
                % MRT Transform
                X_MRT_block = MRT_mat*X_block;
            end
        end
        
        
        function display_MRT(mrt_obj,X_MRT)
            % Display the result of an MRT transform as return by MRT.transform
            
            nb_frames = size(X_MRT,2);
            time_in_s = (0:nb_frames-1)*mrt_obj.R/mrt_obj.Fs;
            frequency_index = 1:mrt_obj.K;
            
            fig_handle = figure;
            imagesc(time_in_s,frequency_index,20*log10(abs(X_MRT)));
            set(gca,'YDir','normal');
            xlabel('Time in seconds');
            ylabel(sprintf('Frequency in %s',mrt_obj.freqUnit));
            % Set Y-ticks in frequency units
            mrt_obj.refresh_Yticks(gca);
            % Provide Y-ticks update when zooming
            zoom_handle = zoom(fig_handle);
            set(zoom_handle,'ActionPostCallback',@mrt_obj.display_zoom_ticks);
            set(zoom_handle,'Enable','on');
            
        end
        
        function display_kernels(mrt_obj,display_mode)
            % Display the MRT transform matrix
            %
            % This function will display the transform matrix according to
            % the optionnal input parameter 'display_mode'.
            % For display_mode =
            % - 'MATRIX' (default) the x and y labels are given in samples
            % and bin index respectively
            % - 'TF' the x and y labels are given in seconds and in the specified frequency
            % unit respectively
            
            if nargin==1,
                display_mode = 'MATRIX';
            else
                if ~any(strcmp(display_mode,{'MATRIX','TF'}))
                    display_mode = 'MATRIX';
                    warning('MRT:DISP_KERNELS:MODE',...
                        'Wrong input parameter ''display_mode'' . It should be ''MATRIX'' (default) or ''TF''');
                end
            end
            
            L = size(mrt_obj.MRT_mat,1);
            time_index = 1:L;
            frequency_index = 1:mrt_obj.K;
            if strcmp(display_mode,'TF')
                time_index = (time_index-1)/mrt_obj.Fs;
                xlabel_str = 'Time in s';
                ylabel_str = sprintf('Frequency in %s',mrt_obj.freqUnit);
            else
                xlabel_str = 'Time in samples';
                ylabel_str = 'Frequency bin index';
            end
            
            fig_handle = figure;
            imagesc(time_index,frequency_index,20*log10(abs(mrt_obj.MRT_mat))');
            
            set(gca,'YDir','normal');
            xlabel(xlabel_str);
            ylabel(ylabel_str);
            title('MRT transform matrix / Time-Frequency kernels')
            if strcmp(display_mode,'TF')
                % Set Y-ticks in frequency units
                mrt_obj.refresh_Yticks(gca);
                % Provide Y-ticks update when zooming
                zoom_handle = zoom(fig_handle);
                set(zoom_handle,'ActionPostCallback',@mrt_obj.display_zoom_ticks);
                set(zoom_handle,'Enable','on');
            end
        end
    end
    
    
end % classdef

