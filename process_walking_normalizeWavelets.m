function varargout = process_walking_normalizeWavelets( varargin )
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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>

    % Description the process
    sProcess.Comment     = 'Normalize Wavelet';
    sProcess.FileTag     = '| normalized';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Walking';
    sProcess.Index       = 800;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
		% Sensor types
		sProcess.options.sensortypes.Comment = 'Sensor types or names: ';
		sProcess.options.sensortypes.Type    = 'text';
		sProcess.options.sensortypes.Value   = 'SEEG';
		sProcess.options.normalizationFactor.Comment = 'Normalize over: ';
		sProcess.options.normalizationFactor.Type    = 'text';
		sProcess.options.normalizationFactor.Value   = 'standing';
		sProcess.options.normalizationTypes.Comment = 'Normalize Type: ';
		sProcess.options.normalizationTypes.Type    = 'text';
		sProcess.options.normalizationTypes.Value   = 'zscore';


end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

	nFiles = numel(sInputs);

	fileOrder = 1:nFiles;

	% we should divide files in order to isolate
	% walking trials from standing trials.

	conditionOrder = {sInputs.Condition};

	normalizationCondition =sProcess.options.normalizationTypes.Value ;
	standingIdx = find(~cellfun(@isempty,regexp(conditionOrder,normalizationCondition )));

	baselineStruct = in_bst_timefreq(sInputs(standingIdx).FileName);

	% channels x F
	meanBaseline = mean(baselineStruct.TF,2);
	stdBaseline  = std(baselineStruct.TF,2);

	walkingIndexes = setdiff(fileOrder,standingIdx);

	for fileIdx = walkingIndexes(:)'

			walkingStruct = in_bst_timefreq(sInputs(fileIdx).FileName); 
			walking = walkingStruct.TF;

			if strcmp(sProcess.options.normalizationTypes.Value,'zscore')
				walking = bsxfun(@rdivide,bsxfun(@minus,walking,meanBaseline),stdBaseline);
			else
				walking = bsxfun(@rdivide,bsxfun(@minus,walking,meanBaseline),meanBaseline);
			end

			% Get the output study (pick the one from the first file)
			iStudy = sInputs(fileIdx).iStudy;

			% Create a new data file structure
			DataMat 						= walkingStruct;
			DataMat.TF          = walking;
			DataMat.DataType    = 'data';
			DataMat.Comment			= strcat('n_',normalizationCondition);
			
			% Create a default output filename 
			OutputFiles{fileIdx} = bst_process('GetNewFilename', ...
					fileparts(sInputs(fileIdx).FileName), 'timefreq');

			% Save on disk
			save(OutputFiles{fileIdx}, '-struct', 'DataMat');
			% Register in database
			db_add_data(iStudy, OutputFiles{fileIdx}, DataMat);
	end


end
