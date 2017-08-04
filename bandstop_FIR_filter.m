function [x, FiltSpec, Messages] = bandstop_FIR_filter(x, Fs, notch_freq_list,transition_band,attenuation, isMirror, Function)
% BST_BANDPASS_HFILTER Linear phase FIR bandpass filter.
%
% USAGE:  [x, FiltSpec, Messages] = bst_bandpass_hfilter(x,  Fs, HighPass, LowPass, isMirror=0, isRelax=0, Function=[detect], Rolloff=[])
%         [~, FiltSpec, Messages] = bst_bandpass_hfilter([], Fs, HighPass, LowPass, isMirror=0, isRelax=0, Function=[detect], Rolloff=[])
%                               x = bst_bandpass_hfilter(x,  Fs, FiltSpec)
%
% DESCRIPTION:
%    - A linear phase FIR filter is created.
%    - Function "kaiserord" and "kaiser" are used to set the necessary order for fir1.
%    - The transition band is hard-coded.
%    - Requires Signal Processing Toolbox for the following functions:
%      kaiserord, kaiser, fir1, fftfilt. If not, using Octave-based alternatives.
%
% INPUT:
%    - x        : [nChannels,nTime] input signal  (empty to only get the filter specs)
%    - Fs       : Sampling frequency
%    - HighPass : Frequency below this value are filtered in Hz (set to 0 for low-pass filter only)
%    - LowPass  : Frequency above this value are filtered in Hz (set to 0 for high-pass filter only)
%    - isMirror : isMirror (default = 0 no mirroring)
%    - isRelax  : Change ripple and attenuation coefficients (default=0 no
%                 relaxation (Ripple = 10^(-3); Atten = 10^(-3) Equals 60db);
%                 1 relaxation (Ripple = 10^(-2); Atten  = 10^(-2) Equals 40db)
%    - Function : 'fftfilt', filtering in frequency domain (default)
%                 'filter', filtering in time domain
%                 If not specified, detects automatically the fastest option based on the filter order
%    - Rolloff  : Width of the transition band in Hz
%
% OUTPUT:
%    - x        : Filtered signals
%    - FiltSpec : Filter specifications (coefficients, length, ...)
%    - Messages : Warning messages, if any

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Hossein Shahabi, Francois Tadel, John Mosher, Richard Leahy, 2016


%% ===== PARSE INPUTS =====
% Filter is already computed
if (nargin == 3) && isstruct(notch_freq_list)
    FiltSpec = notch_freq_list;
    % Default filter options
else
    if (nargin < 4) || isempty(transition_band)
        transition_band = [3 3];
    end
    if (nargin < 5) || isempty(attenuation)
        attenuation = 60;%dB
    end
    if (nargin < 6) || isempty(isMirror)
        isMirror = 0;
    end
    if (nargin < 7) || isempty(Function)
        Function = [];  % Auto-detection based on the filter order later in the code
    end
    
    FiltSpec = [];
end
Messages = [];


%% ===== CREATE FILTER =====
if isempty(FiltSpec)
    % ===== FILTER SPECIFICATIONS =====
    Nyquist = Fs/2;
    notch_freq_list(notch_freq_list==0 & notch_freq_list >= Nyquist) = [];
    % High-pass filter
    if isempty(notch_freq_list)
        Messages = ['No frequency band in input.' 10];
        return;
    end
    
    [b,a,n] = notchFIRfilter_design(notch_freq_list/Nyquist,transition_band/Nyquist,attenuation);
    
    % Filtering function: Detect the fastest option, if not explicitely defined
    if isempty(Function)
        % The filter() function is a bit faster for low-order filters, but much slower for high-order filters
        if any(n > 800)  % Empirical threshold
            Function = 'fftfilt';
        else
            Function = 'filter';
        end
    end
    
    % Output structure
    FiltSpec.b              = b;
    FiltSpec.a              = a;
    FiltSpec.order          = n;
    FiltSpec.trans_width    = transition_band;
    FiltSpec.f_notch_freq     = notch_freq_list/Nyquist;
    FiltSpec.f_notch_freq_Hz  = notch_freq_list ;   % Stop and pass bands in Hz (instead of normalized)
    FiltSpec.function       = Function;
    FiltSpec.mirror         = isMirror;
    % If empty input: just return the filter specs
    if isempty(x)
        return;
    end
end

    %% ===== FILTER SIGNALS =====
    % Transpose signal: [time,channels]
    [nTime,nChan] = size(x);
    for filter_idx = 1:numel(FiltSpec.b)% Half of filter length
        M = FiltSpec.order(filter_idx) / 2;
        % If filter length > 10% of data length
        filter_coverage = nTime/M;
        if (filter_coverage < 3)
            Messages = [Messages, sprintf('The signal is too short for the order of the filter\n')];
        end
        
        % Mirroring requires the data to be longer than the filter
        if (FiltSpec.mirror) && (nTime < M)
            Messages = [Messages, 'Warning: Data is too short for mirroring. Option is ignored...' 10];
            FiltSpec.mirror = 0;
        end
        % Mirror signals
        if FiltSpec.mirror
            x = [ x; flipud(x(end-M+1:end,:))];
            % Zero-padding
        else
            x = [ x; zeros(M,nChan)] ;
        end
        
        % Filter signals
        switch (FiltSpec.function)
            case 'fftfilt'
                x = fftfilt(FiltSpec.b{filter_idx}, x);
            case 'filter'
                x = filter(FiltSpec.b{filter_idx}, FiltSpec.a{filter_idx}, x);
        end
        
        % Remove extra data
        x = x(M+1:end,:);
    end
end

function [b,a,n] = notchFIRfilter_design(fo,tw,attenutation)

Ripple = 10^(-attenutation/20);

for f_idx = 1:numel(fo)
    f_highstop = fo(f_idx)-tw(1)/2-tw(2);
    f_highpass = fo(f_idx)-tw(1)/2;
    f_lowpass = fo(f_idx)+tw(1)/2;
    f_lowstop = fo(f_idx)+tw(1)/2+tw(2);

    fcuts = [f_highstop, f_highpass, f_lowpass, f_lowstop];
    mags  = [1 0 1];               % filter magnitudes
    devs  = [Ripple Ripple Ripple];  % deviations
    fcuts = max(0,fcuts);      % Can't go below zero
    fcuts = min(1-eps, fcuts); % Can't go above or equal to 1
    
    [n(f_idx),Wn,beta,ftype] = kaiserord(fcuts, mags, devs, 2);
    n(f_idx) = n(f_idx) + rem(n(f_idx),2);  % ensure even order
    b{f_idx} = fir1(n(f_idx), Wn, ftype, kaiser(n(f_idx)+1,beta), 'noscale');
    a{f_idx} = 1;
end
end
