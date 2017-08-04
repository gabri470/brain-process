function [x, FiltSpec, Messages] = bandpass_FIR_filter(x, Fs, HighPass,LowPass,transitionBand,attenuation,removeDC, isMirror, Function)
% BST_BANDPASS_HFILTER Linear phase FIR bandpass filter.
%
% USAGE:  [x, FiltSpec, Messages] = bstBandpass_hfilter(x,  Fs, HighPass, LowPass, transitionBand,attenuation=60,removeDC=1, isMirror=0, Function)
%         [~, FiltSpec, Messages] = bstBandpass_hfilter([], Fs, HighPass, transitionBand, attenuation=60, removeDC=1, isMirror=0, Function)
%                               x = bstBandpass_hfilter(x,  Fs, FiltSpec)
%
% DESCRIPTION:
%    - A linear phase FIR filter is created.
%    - Requires Signal Processing Toolbox for the following functions:
%      kaiserord, kaiser, fir1, fftfilt. If not, using Octave-based alternatives.
%
% INPUT:
%    - x        : [nChannels,nTime] input signal  (empty to only get the filter specs)
%    - Fs       : Sampling frequency
%    - HighPass : Frequency below this value are filtered in Hz (set to 0 for low-pass filter only) 
%    							if HighPass is a 2D vector of [HP,LP], the function creates a band-pass filter
%    - LowPass  : Frequency above this value are filtered in Hz (set to 0 for high-pass filter only)
%    - transitionBand  : blablalbalblalabalblalbla
%    - attenuation: dB
%    - removeDc : indovinaindovinello
%    - isMirror : isMirror (default = 0 no mirroring)
%    - Function : 'fftfilt', filtering in frequency domain (default)
%                 'filter', filtering in time domain
%                 If not specified, detects automatically the fastest option based on the filter order
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
% ==============================================================================
%
% Authors: Hossein Shahabi, Francois Tadel, John Mosher, Richard Leahy, 2016


%% ===== PARSE INPUTS =====
% Filter is already computed
if (nargin == 3) && isstruct(HighPass)
    FiltSpec = HighPass;
		removeDC = 1;
    % Default filter options
elseif nargin == 4 && isstruct(HighPass)
		FiltSpec = HighPass;
		removeDC = LowPass;
else
	
    joint_filter = numel(HighPass)-1;
    if (nargin < 5-joint_filter) || isempty(transitionBand)
        transitionBand = [];
        LowPass = HighPass(2);
        HighPass = HighPass(1);
    end
    if (nargin < 6-joint_filter) || isempty(attenuation)
        attenuation = 60;%dB
    end
    if (nargin < 7-joint_filter) || isempty(removeDC)
        removeDC = 1;
    end
    if (nargin < 8-joint_filter) || isempty(isMirror)
        isMirror = 0;
    end
    if (nargin < 9-joint_filter) || isempty(Function)
        Function = [];  % Auto-detection based on the filter order later in the code
    end
    
    FiltSpec = [];
end
Messages = [];


%% ===== CREATE FILTER =====
if isempty(FiltSpec)
    % ===== FILTER SPECIFICATIONS =====
    if(joint_filter)
        transitionBand = min(transitionBand);
    else
        if(numel(transitionBand)==1)
            transitionBand = transitionBand([1 1]);
        end
    end
    Nyquist = Fs/2;
    LowPass(not(isfinite(LowPass)) || LowPass <0) = 0;
    HighPass(not(isfinite(HighPass)) || HighPass <0) = 0;
    % High-pass filter
    if ~isempty(HighPass) && (HighPass ~= 0)
        f_highpass = HighPass / Nyquist;    % Change frequency from Hz to normalized scale (0-1)
        % Default transition band
        if isempty(transitionBand)
            f_highstop = max(0.2, 0.5*HighPass) / Nyquist;
            % Specified manually
        else
            f_highstop = max(0.01, HighPass - transitionBand(1)) / Nyquist;
        end
    else
        f_highpass = 0;
        f_highstop = 0;
    end
    % Low-pass filter
    if ~isempty(LowPass) && (LowPass ~= 0)
        f_lowpass = LowPass / Nyquist;
        % Default transition band
        if isempty(transitionBand)
            if f_highpass==0    % If this is a low-pass filter
                f_lowstop  = 1.05 * f_lowpass;
            else
                f_lowstop  = 1.15 * f_lowpass;
            end
            % Specified manually
        else
            f_lowstop = min(Fs/2 - 0.2, LowPass + transitionBand(2)) / Nyquist;
        end
    else
        f_lowpass  = 0;
        f_lowstop  = 0;
    end
    
    % If both high-pass and low-pass are zero
    if (f_highpass == 0) && (f_lowpass == 0)
        Messages = ['No frequency band in input.' 10];
        return;
        % Input frequencies are too high
    elseif (f_highpass >= 1) || (f_lowpass >= 1)
        Messages = sprintf('Cannot filter above %dHz.\n', Nyquist);
        return;
    end
    % Transition parameters
    Ripple = 10^(-attenuation/20);
    
    % ===== DESIGN FILTER =====
    % Build the general case first
    fcuts = [f_highstop, f_highpass, f_lowpass, f_lowstop];
    mags  = [0 1 0];               % filter magnitudes
    devs  = [Ripple Ripple Ripple];  % deviations
    % Now adjust for desired properties
    fcuts = max(0,fcuts);      % Can't go below zero
    fcuts = min(1-eps, fcuts); % Can't go above or equal to 1
    
    if not(joint_filter)
        % We have implicitly created a bandpass, but now adjust for desired filter
        if(not(f_lowpass == 0) && not(f_highpass == 0))
            [n,Wn,beta,ftype] = kaiserord(fcuts(1:2), mags(1:2), devs(1:2), 2);
            n_high = n + rem(n,2);  % ensure even order
            b_high = fir1(n, Wn, ftype, kaiser(n+1,beta), 'noscale');
            [n,Wn,beta,ftype] = kaiserord(fcuts(3:4), mags(2:3), devs(2:3), 2);
            n_low = n + rem(n,2);  % ensure even order
            b_low = fir1(n, Wn, ftype, kaiser(n+1,beta), 'noscale');
            n = [n_high n_low];
            b = {b_high,b_low};
            a = {1,1};
        else
            if (f_lowpass == 0)  % User didn't want a lowpass
                fcuts(3:4) = [];
                mags(3) = [];
                devs(3) = [];
            end
            if (f_highpass == 0)  % User didn't want a highpass
                fcuts(1:2) = [];
                mags(1) = [];
                devs(1) = [];
            end
            [n,Wn,beta,ftype] = kaiserord(fcuts, mags, devs, 2);
            n = n + rem(n,2);  % ensure even order
            b = fir1(n, Wn, ftype, kaiser(n+1,beta), 'noscale');
            a = {1};
        end
    else
        [n,Wn,beta,ftype] = kaiserord(fcuts, mags, devs, 2);
        n = n + rem(n,2);  % ensure even order
        b = fir1(n, Wn, ftype, kaiser(n+1,beta), 'noscale');
        a = {1};
    end
    
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
    FiltSpec.f_highpass     = f_highpass;
    FiltSpec.f_lowpass      = f_lowpass;
    FiltSpec.fcuts          = fcuts * Nyquist ;   % Stop and pass bands in Hz (instead of normalized)
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
    % Remove the mean of the data before filtering
    if removeDC
        xmean = mean(x);
        x = bsxfun(@minus, x, xmean);
    end
    
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
    % Restore the mean of the signal (only if there is no high-pass filter)
    if (FiltSpec.f_highpass == 0) && removeDC
        x = bsxfun(@plus, x, xmean);
    end
    
