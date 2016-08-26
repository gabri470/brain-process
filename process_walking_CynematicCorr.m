function varargout = process_walking_CinematicCorr( varargin )
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

varargout = {};

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Correlation w\ cynematics';
    sProcess.FileTag     = '__';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Walking';
    sProcess.Index       = 801;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
		% Sensor types


		sProcess.options.outputLabel.Type = 'text';
		sProcess.options.outputLabel.Comment = ' Output string ';
		sProcess.options.outputLabel.Value = '';


   


end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
%
%	nFiles = numel(sInputs);
%
	DATA_FOLDER = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','info');
  trialRejectionFile = fullfile(DATA_FOLDER,'trialRejection.csv');
  [patNames,trialStrings,stepIds] = textread(trialRejectionFile,...
																				'%s %*s %*s %*s%s%d%*s','delimiter',',');

	

		strOut = sProcess.options.outputLabel.Value;

		fnameOut = fullfile(DATA_FOLDER,strOut);
		fid = fopen(fnameOut,'w');

    % Initialize returned list of files
    OutputFiles = {};

		for file = 1:numel(sInputs)

				% read data in
				DataMat = in_bst_timefreq(sInputs(file).FileName);
				parentStruct = bst_process('GetInputStruct',DataMat.DataFile);
				parentData = in_bst(parentStruct.FileName);


				% filter step rejection list for
				% a given subject
				subjMask = ismember(lower(patNames),lower(sInputs(file).SubjectName));
				trialString = regexp(sInputs(file).Condition,'trial\d+','match');
				trialMask = ~cellfun(@isempty,regexp(trialStrings,trialString));
				stepReject = stepIds(subjMask & trialMask);
				
				% extract event names
				evNames = {parentData.Events.label};
				
				% search event heel L and heel R
				leftHeelEvIdx = strcmp(evNames,'heelcontact_L');
				rightHeelEvIdx = strcmp(evNames,'heelcontact_R');

				fs = round(1/mean(diff( parentData.Time )));

				% compute offset
				offset = round(parentData.Time(1)*fs);

				leftHeelEvSamples = parentData.Events(leftHeelEvIdx).samples;
				rightHeelEvSamples = parentData.Events(rightHeelEvIdx).samples;

				% we create a lookup table for foot order
				% hc_L marks a right foot swing 
				% on the other hand, hc_R marks a left foot swing 
				% thus we switch labels at this level
				footLabels = [ repmat({'R'},1,numel(leftHeelEvSamples)), ...
						repmat({'L'},1,numel(rightHeelEvSamples))];

				[heelEvSamples, ord] = sort([leftHeelEvSamples(:);rightHeelEvSamples(:)]);

				% since Brainstorm saves samples from the beginning of the original file
				% even though we extracted a portion of it, we need recompute the exact
				% samples within the time window of interest
				heelEvSamples = heelEvSamples - offset;

				footLabels = [footLabels{ord}];

				stringTokens = regexp(parentStruct.Condition,'_','split');

				betaBand = DataMat.Freqs >= 13 & DataMat.Freqs < 30;
				gammaBand = DataMat.Freqs >= 30 & DataMat.Freqs < 60;

				nSwings = numel(leftHeelEvSamples) + numel(rightHeelEvSamples)-1;

				betaPower = squeeze(mean(DataMat.TF(:,:,betaBand),3));
				gammaPower = squeeze(mean(DataMat.TF(:,:,gammaBand),3));

				swingIndicies = setdiff(1:nSwings,stepReject);

				for iSwing = swingIndicies

						betaInStep = sum(betaPower(:,heelEvSamples(iSwing):heelEvSamples(iSwing+1)),2).*1e12;
						gammaInStep = sum(gammaPower(:,heelEvSamples(iSwing):heelEvSamples(iSwing+1)),2).*1e12;
						
						% subject, condition, speed, drug, trial, swingN, foot, betaR,betaL
						fprintf(fid,'%s %s %s %s %s %d %s %e %e %e %e\n',...
								stringTokens{1:5},iSwing,footLabels(iSwing),betaInStep(:),gammaInStep(:));

				end
			

		end
    
		fclose(fid);

	
end % function
