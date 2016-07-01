function varargout = process_cleanWrongFilterSettings( varargin )
% PROCESS_EXAMPLE_CUSTOMAVG: Example file that reads all the data files in input, and saves the average.

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
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
% Authors: Francois Tadel, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Clear Filter Settings';
    sProcess.FileTag     = '| cleaned';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Pre-process';
    sProcess.Index       = 123;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq','data'};
    sProcess.OutputTypes = {'timefreq','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};

	% BRUTAL HACKING %	
	DataMat1 = in_bst_timefreq(sInputs(1).FileName);	
	DataMat2 = in_bst_timefreq(sInputs(2).FileName);	
	DataMat3 = in_bst_timefreq(sInputs(3).FileName);	

	data1 = squeeze(log10(DataMat1.TF));
	data2 = squeeze(log10(DataMat2.TF));
	data3 = squeeze(log10(DataMat3.TF));

	data = data3 - data1 + data2;

	DataMat3.TF = permute(10.^data,[1 3 2]);

    OutputFiles{1} = bst_process('GetNewFilename',...
		fileparts(sInputs(3).FileName), 'timefreq_psd');

    % Save on disk
    save(OutputFiles{1}, '-struct', 'DataMat3');

    % Register in database
    db_add_data(sInputs(3).iStudy, OutputFiles{1}, DataMat3);

end
