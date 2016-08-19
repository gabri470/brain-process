function varargout = process_walking_StandingVsWalking( varargin )
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
varargout =  {};
macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Standing vs Walking';
    sProcess.FileTag     = '__';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Walking';
    sProcess.Index       = 801;
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

	DATA_FOLDER = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','info');

	sideFile = fullfile(DATA_FOLDER,'patientSides.csv');
	[nameSubjects, mostAffSides] = textread(sideFile,'%s %s\n','delimiter',',');

	subjectNames = unique({sInputs.SubjectName});
	nSubjects = numel(subjectNames);

	mostAffSides = mostAffSides(ismember(nameSubjects,subjectNames));

	% we need to find files that have to be processed and group 
	% them for subjects and conditions
	% First group all sInput.Comment together
	conditionStrings = {sInputs.Comment};
	
	standingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(s|S)tanding'));
	walkingConMask = ~cellfun(@isempty,regexp(conditionStrings,'Original.+SEEG'));

	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
							 'paperunits','normalized','paperorientation',...
								'portrait','Visible','on');


	OutputFiles = {};

	for subjIdx = 1:nSubjects
			% for each subject separately we pick standing condition
			subjectMask = ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));

			standingFileIdx = find(subjectMask & standingConMask);
			walkingFileIdx = find(subjectMask & walkingConMask);

			% read standing time-freq data
			standingStruct = in_bst_timefreq(sInputs(standingFileIdx).FileName);

			% we should cycle through walking trials
			% compute the min length in order to do 
			% the statistics lately
			for walkTrialIdx = 1:numel(walkingFileIdx)

					walkingStruct = in_bst_timefreq(sInputs(walkingFileIdx(walkTrialIdx)).FileName);
					parentStruct 	= bst_process('GetInputStruct',walkingStruct.DataFile);
					parentData 		= in_bst(parentStruct.FileName);
		
					% compute sampling frequency
					fs = round(1/mean(diff( parentData.Time )));

					% filter cardiac and peakVelocity events from gait-related events
					evGroupNames = {parentData.Events.label};
					gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel|toe)'));

					% concat all heel contact events in order to have
					% a vector of latencies of this form: e.g.
					% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
					eventSamples = sort([parentData.Events(gaitEventGroups).samples]);

					% get the maximum duration from first heel contact 
					% to the last heel contact 
					walkingLength(walkTrialIdx) = max(eventSamples)-min(eventSamples);

			end % walking loop

			% we then get the min windows
		  walkingRefLength = min(walkingLength);

			standingLength = numel(standingStruct.Time);

			% how many windows of walking length can 
			% we extract from standing data 
			% in order to make test for difference?
			nWindows = floor(standingLength/walkingRefLength);

			% we separate chucks of standing of walkingRefLength in order to 
			% make periods with the same time lengths 
			standData = reshape(standingStruct.TF(:,1:(nWindows*walkingRefLength),:),...
																[2,nWindows,walkingRefLength,90]);


			walkData = nan(2,walkingRefLength,90,numel(walkingFileIdx));

			for walkIdx = 1:numel(walkingFileIdx)
					walkingStruct = in_bst_timefreq(sInputs(walkingFileIdx(walkTrialIdx)).FileName);
					parentStruct 	= bst_process('GetInputStruct',walkingStruct.DataFile);
					parentData 		= in_bst(parentStruct.FileName);
		
					% filter cardiac and peakVelocity events from gait-related events
					evGroupNames = {parentData.Events.label};
					gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel|toe)'));

					% concat all heel contact events in order to have
					% a vector of latencies of this form: e.g.
					% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
					eventSamples = sort([parentData.Events(gaitEventGroups).samples]);

					% we have to correct the event adding the offset
					% since they are referred to the 0 of the raw data 
					evOffset 			= round(walkingStruct.Time(1)*fs);
					eventSamples  = eventSamples - evOffset;

					% we then pack trails together in order to have a matrix 2 x walkingRefLength x f x trials
					% that we will rotate to match the form 2 x windows x walkingRefLength x f as 
					% standing condition data are represented in
					walkData(:,:,:,walkIdx) = walkingStruct.TF(:,eventSamples(1)+(1:walkingRefLength),:);

		  end % walking trial loop

			% permute dim order to match 2 x #windows x length x f
			walkData = permute(walkData,[1 4 2 3]);

			% compute the PSD as mean over time of the morlet coefficients
			% for stand which will result in a 2 x nWindows x 90 matrix
			standPsd = squeeze(mean(standData,3));
			% and walking that will result in a 2 x trials x 90 points
			walkPsd = squeeze(mean(walkData,3));

			% we should normalize 'em all as metallica are singing now
			standPsd = bsxfun(@rdivide,standPsd,sum(standPsd,3));
			walkPsd = bsxfun(@rdivide,walkPsd,sum(walkPsd,3));

			% then compute a relative change to test for the hypothesis
			% that standing would have more power in beta band then walking
			% since walking should be suppressing beta and gamma during strides
			relChange = (squeeze(mean(walkPsd,2)) - squeeze(mean(standPsd,2)))./squeeze(mean(standPsd,2));
			
			% then separately for STN- and STN+
			% compute the difference and do a permutation test
			% we should of course correct for multiple comparisons
			pvalue(1,:) = runPermutationTest(relChange(1,:),walkPsd(1,:,:),standPsd(1,:,:),100,0.05);
			pvalue(2,:) = runPermutationTest(relChange(2,:),walkPsd(2,:,:),standPsd(2,:,:),100,0.05);
			
			
			% sort rows (ie. channels) to match STN- / STN+ order
			% that is fixed thorought the paper/analyses
			if(strcmp(mostAffSides(subjIdx),'L'))
					relChange = relChange([2 1],:);
			end

			significanceMask = nan(2,90);
			significanceMask(pvalue < 0.05) = 3;

			subplot(nSubjects,4,4*(subjIdx-1)+1,'NextPlot','add')
			plot(1:90,relChange(1,:))
			plot(1:90, significanceMask(1,:),'r.','MarkerSize',10);%,'EdgeColor',[1 0 0]);
			xlim([6 50]);
			ylim([-1 4]);

			subplot(nSubjects,4,4*(subjIdx-1)+2,'NextPlot','add')
			plot(1:90,squeeze(mean(walkPsd(1,:,:),2)),'r'),
			plot(1:90,squeeze(mean(standPsd(1,:,:),2)),'b');
			xlim([6 50]);
			
			subplot(nSubjects,4,4*(subjIdx-1)+3,'NextPlot','add')
			plot(1:90,relChange(2,:))
			plot(1:90, significanceMask(2,:),'r.','MarkerSize',10);%,'EdgeColor',[1 0 0]);
			xlim([6 50]);
			ylim([-1 4]);

			subplot(nSubjects,4,4*(subjIdx-1)+4,'NextPlot','add')
			plot(1:90,squeeze(mean(walkPsd(2,:,:),2)),'r'),
			plot(1:90,squeeze(mean(standPsd(2,:,:),2)),'b');
			xlim([6 50]);


	end % subject loop
	clearvars walkData standData
end % function


function [pvalue, unCorrpvalue] = runPermutationTest(dataObs,dataA,dataB,nPermutation,alpha)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%
		pvalue  	 = zeros(1,90);
		nStanding  = size(dataB,2);
		pooledData = cat(2,dataA,dataB);

		% we perform a permutation test for each STN separatelly
		for permIdx = 1:nPermutation

			riffledIndices = randperm(size(pooledData,2));

			standPsd = pooledData(:,riffledIndices(1:nStanding),:);
			walkPsd	= pooledData(:,riffledIndices(nStanding+1:end),:);
			
			% compute permutated statistics
			dataPerm = (squeeze(mean(walkPsd,2)) - squeeze(mean(standPsd,2)))./squeeze(mean(standPsd,2));
		
			% compute pvalues for all frequencies and all time points.
			pvalue = pvalue + double(dataPerm' > dataObs)./nPermutation;

		end
		unCorrpvalue = pvalue;
		pvalue = fdrCorrection(pvalue,alpha);

end


function pvalue = fdrCorrection(pvalue, alpha)
%FDRCORRECTION Description
%	PVALUE  = FDRCORRECTION() Long description
%

	tmpPvalue 	= sort(pvalue(:));
	N 					= numel(pvalue);
	FDR 				= alpha.*(1:N)./N;
	thr 				= FDR(find(tmpPvalue <= FDR',1,'last'));
	pvalue(pvalue >= thr) = 1;

end
