function varargout = process_Montage_betaGammaRatio( varargin )
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
    sProcess.Comment     = 'Fig. 5';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Plot';
    sProcess.Index       = 1003;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq','data'};
    sProcess.OutputTypes = {'timefreq','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
	% Sensor types
	sProcess.options.normalize.Comment = 'Normalize PSDs';
	sProcess.options.normalize.Type = 'checkbox';
	sProcess.options.normalize.Value = 1;

	sProcess.options.chType.Comment = 'Channel Type: ';
	sProcess.options.chType.Value = 'MISC';
	sProcess.options.chType.Type = 'text';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
	OutputFiles = {};
	
	chType = sProcess.options.chType.Value;

	figure(1);
	colorOrder = jet(6);
	set(gcf,'DefaultAxesColorOrder',colorOrder);
%	sPlotOrder = reshape(1:9,[3 3])';
%	sPlotOrder = mat2cell(sPlotOrder,ones(9,1),1);
	sPlotOrder = 1:9;

	tmp = [2 3 4 5 6 7 9 10 11];
	subjectNames = arrayfun(@(x) sprintf('wue%02d',x),tmp,'uni',false);

	avg = zeros(2, 129);
	sqS = zeros(2, 129);

	nFiles = numel(sInputs);
	
	groupMat = zeros(numel(tmp),2,2);
	
    % Reading all the input files in a big matrix
    for i = 1:length(sInputs)

		betaGammaRatio = zeros(2,2);

		% we should check that we are effectively comparing 
		% the correct studies
		[parStruct, subjName, med, visit,...
			condTag,rConditionTag] = extractInfo(sInputs(i));

	    pConditionTag = regexprep(rConditionTag,'rest','ProSup');

		pCondition = bst_get('StudyWithCondition',strcat(subjName,...
			filesep,pConditionTag));
	
		if(~isempty(pCondition.Timefreq))

			subjIdx = find(strcmp(subjectNames,subjName));
			Channels = in_bst_channel(sInputs(i).ChannelFile,'Channel');
			[chLabels, chOrder] = sortLabels({Channels.Channel.Name});

			restMat = in_bst_timefreq(sInputs(i).FileName);
			prosMat = in_bst_timefreq(pCondition.Timefreq.FileName);

			rest = normalizePSD(restMat.TF);
			pros = normalizePSD(prosMat.TF);

			rest = filterTimeFreqMat(rest,Channels,chType,chOrder);
			pros = filterTimeFreqMat(pros,Channels,chType,chOrder);
			
			fAxis = restMat.Freqs;

			gammaBand = fAxis >= 31 & fAxis <= 80; 
			betaBand  = fAxis >= 13 & fAxis <= 30; 

			restRatio = mean(rest(:,betaBand),2) ./ mean(rest(:,gammaBand),2);
			prosRatio = mean(pros(:,betaBand),2) ./ mean(pros(:,gammaBand),2);
			
			groupMat(subjIdx,:,:) = [restRatio, prosRatio];	

			subplot(3,3,sPlotOrder(subjIdx))
			bar(1:2,[restRatio, prosRatio]);

			chLabelsT = [{Channels.Channel.Name}];
			cChLabels = chLabelsT(strcmp({Channels.Channel.Type},chType)); 
			cChLabels = chLabels(ismember(chLabels,cChLabels));

			cChLabels = cellfun(@(x) x(1:end-5),cChLabels,'UniformOutput',false);
			set(gca,'XTickLabel',cChLabels);

			title(subjName);			

		else
			continue
		end
    end
	avgRatio = reshape(permute(groupMat,[3 1 2]),[2 18]);
	stdRatio = reshape(permute(groupMat,[3 1 2]),[2 18]);
	avgRatio = mean(avgRatio,2);
	stdRatio = std(stdRatio,[],2)./sqrt(18);
	
	figure(2)
	bar(1:2,avgRatio);
	title('Population Average');
	set(gca,'XTickLabel',[{'rest'} {'pros'}]);

	hold on;
	errorbar(1:2,avgRatio,stdRatio,'LineStyle','none');
    

end

function [parentStruct, subjectName,med, visit,condTag, conditionTag] = extractInfo(in)
%EXTRACTINFO Description
%	[PARENTSTRUCT, SUBJECTNAME,MED, VISIT,CONDTAG] = EXTRACTINFO(IN) Long description
%
	parentStruct = bst_process('GetInputStruct',in.DataFile); 
	subjectName = parentStruct.SubjectName;
	med = regexp(parentStruct.Condition,'ON|OFF','match');
	visit = regexp(parentStruct.Condition,'visite_\d+','match');
	condTag = regexp(parentStruct.Condition,'rest|ProSup|walking','match');
	conditionTag = parentStruct.Condition;
end

function normPSD = normalizePSD(normPSD)
%NORMALIZEPSD Description
%	NORMPSD = NORMALIZEPSD(PSD) Long description
%
	normPSD = squeeze(normPSD);
	normPSD = normPSD ./ repmat(sum(normPSD,2),[1,size(normPSD,2)]);
end

function [labels, order]= sortLabels(labels)
%SORTLABELS in left-right groups, and in alphabetical order within each group
%	[LABELS,ORDER] = SORTLABELS(LABELS) Long description
%
	
	tmp = regexp(labels,'\d+','match');
	mask = sum(reshape(cellfun(@str2num,[tmp{:}]),[2, 12]));
	rightChannels = mask < 6;
	leftChannels = mask > 6;
	tmpA = sort(labels(rightChannels));
	tmpB = sort(labels(leftChannels));

	tmp = [tmpB tmpA];

	[~,order] = ismember(tmp, labels);
	labels = tmp;

end

function timefreq = filterTimeFreqMat(timefreq,Channels,chType,chOrder)
%FILTERTIMEFREQMAT Description
%	TIMEFREQ = FILTERTIMEFREQMAT(TIMEFREQ,CHANNELS,CHTYPE,CHORDER) Long description
%

	timefreq = squeeze(timefreq);

	Mask = strcmp({Channels.Channel.Type},chType);

	timefreq(~Mask,:) = NaN;

	timefreq = timefreq(chOrder,:);

	timefreq(sum(isnan(timefreq),2) > 1,:) = [];

end


function [avg sqsum] = computeStat(psd,avg,sqsum)
%COMPUTESTAT Description
%	[AVG STD] = COMPUTESTAT(PSD,AVG) Long description
%
	
	% remove the nan channels
	psd(sum(isnan(psd),2) > 1,:) = [];

	avg = avg + psd;
	sqsum  =  sqsum + psd.^2;
end
