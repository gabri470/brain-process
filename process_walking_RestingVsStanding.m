function varargout = process_walking_RestingVsStanding( varargin )
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
    sProcess.Comment     = 'Resting Vs Standing';
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
	restingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(r|R)esting'));

	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
							 'paperunits','normalized','paperorientation',...
								'portrait','Visible','on');


	OutputFiles = {};

	for subjIdx = 1:nSubjects
			% for each subject separately we pick standing condition
			subjectMask = ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));

			standingFileIdx = find(subjectMask & standingConMask);
			restingFileIdx = find(subjectMask & restingConMask);

			% read standing time-freq data
			standingStruct 	= in_bst_timefreq(sInputs(standingFileIdx).FileName);
			parentStruct 		= bst_process('GetInputStruct',standingStruct.DataFile);
			standChannels 	= in_bst_channel(parentStruct.ChannelFile);
			standiChannels	= channel_find( standChannels.Channel,'SEEG');

			restingStruct 	= in_bst_timefreq(sInputs(restingFileIdx).FileName);
			parentStruct 		= bst_process('GetInputStruct',restingStruct.DataFile);
			restChannels 		= in_bst_channel(parentStruct.ChannelFile);
			restiChannels		= channel_find( standChannels.Channel,'SEEG');

			standingPsd 		= squeeze(mean(standingStruct.TF(standiChannels,:,:),2));
			restingPsd			= squeeze(mean(restingStruct.TF(restiChannels,:,:),2));

%			relChange				= bsxfun(@rdivide,bsxfun(@minus,standingPsd,mean(restingPsd)),mean(restingPsd));
			relChange 			= (restingPsd - standingPsd)./standingPsd;
			
			% sort rows (ie. channels) to match STN- / STN+ order
			% that is fixed thorought the paper/analyses
			if(strcmp(mostAffSides(subjIdx),'L'))
					relChange = relChange([2 1],:);
			end

			subplot(nSubjects,2,2*(subjIdx-1)+1,'NextPlot','add')
			plot(1:90,relChange(1,:),'r')
%			plot(1:90,significanceMask(1,:),'r.','MarkerSize',10);%,'EdgeColor',[1 0 0]);

			plot(1:90,relChange(2,:),'c')
%			plot(1:90, significanceMask(2,:),'c.','MarkerSize',10);%,'EdgeColor',[1 0 0]);
			xlim([6 40]);
%			ylim([-1 4]);

			restingPsd = bsxfun(@rdivide,restingPsd,sum(restingPsd,2));
			standingPsd = bsxfun(@rdivide,standingPsd,sum(standingPsd,2));

			subplot(nSubjects,2,2*(subjIdx-1)+2,'NextPlot','add')
			plot(1:90,restingPsd(1,:),'r'),
			plot(1:90,standingPsd(1,:),'r--');
			plot(1:90,restingPsd(2,:),'c'),
			plot(1:90,standingPsd(2,:),'c--');
			xlim([6 40]);


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
