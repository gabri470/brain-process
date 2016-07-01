function varargout = process_resample_new( varargin )
% PROCESS_RESAMPLE: Resample matrix with a new sampling frequency.
%
% USAGE:      sProcess = process_resample('GetDescription')
%               sInput = process_resample('Run', sProcess, sInput, method)
%               sInput = process_resample('Run', sProcess, sInput)
%        [x, time_out] = process_resample('Compute', x, time_in, NewRate, method)
%        [x, time_out] = process_resample('Compute', x, time_in, NewRate)
%        [x,Pfac,Qfac] = process_resample('ResampleCascade', x, NewRate, OldRate, Method)

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2015 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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
% Authors: Francois Tadel, 2010-2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Resample interp';
    sProcess.FileTag     = '| resampinterp';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Pre-process';
    sProcess.Index       = 2000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results', 'matrix'};
    sProcess.OutputTypes = {'data', 'results', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Default values for some options
    sProcess.processDim  = 1;    % Process channel by channel
    sProcess.isSeparator = 1;
    
    % Definition of the options
    % === Resample frequency
    sProcess.options.freq.Comment = 'New frequency:  ';
    sProcess.options.freq.Type    = 'value';
    sProcess.options.freq.Value   = {1000,'Hz',2};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sprintf('Resample: %dHz', round(sProcess.options.freq.Value{1}));
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput, method) %#ok<DEFNU>
    % Get method name
    if (nargin < 3)
        method = [];
    end
    % Output frequency
    NewFreq = sProcess.options.freq.Value{1};
    % Check output frequency
    OldFreq = 1 ./ (sInput.TimeVector(2) - sInput.TimeVector(1));
    if (abs(NewFreq - OldFreq) < 0.05)
        bst_report('Error', sProcess, [], 'Sampling frequency was not changed.');
        sInput = [];
        return;
    end
    % Check for Signal Processing toolbox
    if ~bst_get('UseSigProcToolbox')
        bst_report('Warning', sProcess, [], [...
            'The Signal Processing Toolbox is not available. Using the EEGLAB method instead (results may be much less accurate).' 10 ...
            'This method is based on a fft-based low-pass filter, followed by a spline interpolation.' 10 ...
            'Make sure you remove the DC offset before resampling; EEGLAB function does not work well when the signals are not centered.']);
    end
    % Resample
    [sInput.A, sInput.TimeVector] = Compute(sInput.A, sInput.TimeVector, NewFreq);
    % Update file
    sInput.FileTag = sprintf('| resample(%dHz)', round(NewFreq));
    sInput.HistoryComment = sprintf('Resample from %0.2f Hz to %0.2f Hz (%s)', OldFreq, NewFreq, method);
end


%% ===== EXTERNAL CALL =====
% USAGE: [x, time_out] = process_resample('Compute', x, time_in, NewFreq)
% INPUT:
%     - x       : Signal to process [nChannels x nTime]
%     - Time    : Original time vector
%     - NewFreq : New sampling frequency (Hz)
% OUTPUT:
%     - x    : Resampled signal
%     - Time : New time vector
function [x, Time] = Compute(x, t, NewFreq)
 
    % Check output frequency
    OldFreq = 1 ./ (t(2) - t(1));
    if (abs(NewFreq - OldFreq) < 0.05)
        return;
    end

    Time = t(1):1/NewFreq:t(end);
    x=interp1(t,x',Time)';

end






