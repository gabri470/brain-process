function varargout = process_extract_time_raw( varargin )
% PROCESS_EXTRACT_TIME: Extract blocks of data from a set of files.

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Extract time';
    sProcess.FileTag     = '| extract';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'My-Pre-process';
    sProcess.Index       = 353;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/ExploreRecordings?highlight=%28Extract+time%29#Time_selection';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.processDim  = 1;    % Process channel by channel
    sProcess.isSeparator = 1;
    
    % Definition of the options
    % === TIME WINDOW
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    
    sProcess.options.overwrite.Comment = 'Overwrite';
    sProcess.options.overwrite.Type = 'checkbox';
    sProcess.options.overwrite.Value = 0;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = [sProcess.Comment, ': [', GetTimeString(sProcess), ']'];
end

%% ===== GET TIME STRING =====
function strTime = GetTimeString(sProcess, sInput)
    % Get time window
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        time = sProcess.options.timewindow.Value{1};
    elseif (nargin >= 2) && isfield(sInput, 'TimeVector') && ~isempty(sInput.TimeVector)
        time = sInput.TimeVector([1 end]);
    else
        time = [];
    end
    % Print time window
    if ~isempty(time)
        if any(abs(time) > 2)
            if (time(1) == time(2))
                strTime = sprintf('%1.3fs', time(1));
            else
                strTime = sprintf('%1.3fs,%1.3fs', time(1), time(2));
            end
        else
            if (time(1) == time(2))
                strTime = sprintf('%dms', round(time(1)*1000));
            else
                strTime = sprintf('%dms,%dms', round(time(1)*1000), round(time(2)*1000));
            end
        end
    else
        strTime = 'all';
    end
end



%% ===== RUN =====
function Outputfiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Get options
    Outputfiles = {};
    sProcess.Function = @process_extract_time;

    % ===== FILTER DATA =====
    for iInput=1:length(sInputs)
        Outputfiles{iInput} = ProcessFilter(sProcess, sInputs(iInput));
    end
end




