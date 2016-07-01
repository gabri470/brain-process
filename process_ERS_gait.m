function varargout = process_customavg( varargin )
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
%
macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'ERS/ERD average for gait ini';
    sProcess.FileTag     = '| imported';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Pre-process';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
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

	
		conditions   = unique([ {sInputs.Condition} ]);
		subjects   	 = unique([ {sInputs.SubjectName} ]);

		erdControlat = zeros(1,550,40);
		erdIpsilat 	 = zeros(1,550,40);

		% i have to read the foot from csv
		% WUE, MEDS, TRIAL, FOOT, MOVONSET(ignore)
		[subjectNames, meds, trials, foots] = textread('/home/gabri/Gait.csv',...
				'%s %s %d %s %*f','delimiter',',','headerlines',1);

		chControlatIndices = ones(size(trials));
		chIpsilatIndices = ones(size(trials));

		chControlatIndices(cellfun(@strcmp,foots,repmat({ 'right' },size(foots)))) = 3;
		chIpsilatIndices(cellfun(@strcmp,foots,repmat({ 'left' },size(foots)))) = 3;

		% reshape as trial, med, subj
		chControlatIndices = reshape(chControlatIndices,[4,2,3]);
		chIpsilatIndices = reshape(chIpsilatIndices,[4,2,3]);

		f1 = figure();

		for conditionIdx = 1:numel(conditions)
			n = 0;
			for subjectIdx = 1:numel(subjects)

				theCondition = find(strcmp([{sInputs.Condition}],...
						conditions{conditionIdx}));

				theSubject = find(strcmp([{sInputs.SubjectName}],...
						subjects{subjectIdx}));

				fileIndices = intersect(theCondition,theSubject);

				n = n +numel(fileIndices);

				for file = fileIndices
						
						DataMat = in_bst_timefreq(sInputs(file).FileName);
						parentStruct = bst_process('GetInputStruct',DataMat.DataFile);

						trial 	= regexp(parentStruct.Comment,'\d+','match');
						trial		= str2num(trial{:});

						chControlatIdx = chControlatIndices(trial,conditionIdx,subjectIdx);
						chIpsilatIdx = chIpsilatIndices(trial,conditionIdx,subjectIdx);

						erdControlat = erdControlat + DataMat.TF(chControlatIdx,:,:);
						erdIpsilat 	 = erdIpsilat + DataMat.TF(chIpsilatIdx,:,:);

				end
			end

			erdControlat = erdControlat ./n;
			erdIpsilat = erdIpsilat ./n;

			freq = 1:2:80;
			time=linspace(-300,1000,550);

			h(conditionIdx) = subplot(2,2,conditionIdx);
			imagesc(time,freq, squeeze(erdControlat)');
			title(strcat('Controlat ',conditions{conditionIdx}));
			axis xy;
			hold on;
			plot([0 0],[0 80],'k--','LineWidth',1);
			plot([300 300],[0 80],'k--','LineWidth',1);

			h(conditionIdx+2) =	subplot(2,2,conditionIdx+2);
			imagesc(time, freq, squeeze(erdIpsilat)');
			title(strcat('Ipsilat ',conditions{conditionIdx}));
			axis xy;
			hold on;
			plot([0 0],[0 80],'k--','LineWidth',1);
			plot([300 300],[0 80],'k--','LineWidth',1);


		end

		clim = get(h,'clim');
		clim = cat(1,clim{:});
		m = min(clim(:));
		M = min(clim(:,2));

		set(h,'clim',[m M]);
end
