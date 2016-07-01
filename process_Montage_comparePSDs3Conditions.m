function varargout = process_compare_PSDs3Conditions( varargin )
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
    sProcess.Comment     = 'Fig. 4';
    sProcess.FileTag     = '| compared';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Plot';
    sProcess.Index       = 702;

    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq','data'};
    sProcess.OutputTypes = {'timefreq','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Definition of the options
	sProcess.options.relativeChange.Comment = 'Relative Change';
	sProcess.options.relativeChange.Type    = 'checkbox';
	sProcess.options.relativeChange.Value   = 0;


	sProcess.options.chType.Comment 		= 'Channel Type: ';
	sProcess.options.chType.Value 			= 'MISC';
	sProcess.options.chType.Type 			= 'text';
%
%	% Definition of the options
%	SelectOptions = {...
%		'', ...                               % Filename
%		'', ...                               % FileFormat
%		'open', ...                           % Dialog type: {open,save}
%		'Select directory...', ...            % Window title
%		'ImportData', ...                    
%		'single', ...                         % Selection mode: {single,multiple}
%		'dirs', ...                           % Selection mode: {files,dirs}
%		bst_get('FileFilters', 'events'), ... % Get all the available file formats
%		'EventsIn'};                          % DefaultFormats: {ChannelIn,Da}
%
%	% Option: Event file
%	sProcess.options.chFile.Comment = 'Saving directory:';
%	sProcess.options.chFile.Type    = 'filename';
%	sProcess.options.chFile.Value   = SelectOptions;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
    
    % ===== LOAD THE DATA =====
    % Read the first file in the list, to initialize the loop
	

	nFiles 	= numel(sInputs);

	chType 	= sProcess.options.chType.Value;

	figure(1);
	colorOrder = jet(6);
	set(gcf,'DefaultAxesColorOrder',colorOrder);
	sPlotOrder = reshape(1:18,[2 9])';
	sPlotOrder = mat2cell(sPlotOrder,ones(9,1),2);

	tmp = [2 3 4 5 6 7 9 10 11];
	subjectNames = arrayfun(@(x) sprintf('wue%02d',x),tmp,'uni',false);

	prosAvg = zeros(2,129);
	walkAvg = zeros(2,129);
	restAvg = zeros(2,129);

	prosSqSum = zeros(2,129);
	walkSqSum = zeros(2,129);
	restSqSum = zeros(2,129);
	for el = 1:nFiles	


		% we should check that we are effectively comparing 
		% the correct studies
		[parStruct, subjName, med, visit,...
			condTag,rConditionTag] = extractInfo(sInputs(el));


		% we now have to pick the psds of corresponding conditions ProSup and 
		% walking
									
	    pConditionTag = regexprep(rConditionTag,'rest','ProSup');
	    wConditionTag = regexprep(rConditionTag,'rest','walking');

		pCondition = bst_get('StudyWithCondition',strcat(subjName,...
			filesep,pConditionTag));
		wCondition = bst_get('StudyWithCondition',strcat(subjName,...
			filesep,wConditionTag));
	
		if(~isempty(pCondition.Timefreq) & ~isempty(wCondition.Timefreq))

				subjIdx = find(strcmp(subjectNames,subjName));

				Channels = in_bst_channel(sInputs(el).ChannelFile,'Channel');

				% the ChannelType is fixed to MISC for those
				% that have been selected for further analyses 
				% after visual inspection of beta peaks among the 
				% significant ones. Thus we consider only the MISC channels
				% and discarding all the SEEG channels remaning
				[chLabels, chOrder] = sortLabels({Channels.Channel.Name});

				restMat = in_bst_timefreq(sInputs(el).FileName);
				prosMat = in_bst_timefreq(pCondition.Timefreq.FileName);
				walkMat = in_bst_timefreq(wCondition.Timefreq.FileName);

				rest = normalizePSD(restMat.TF);
				pros = normalizePSD(prosMat.TF);
				walk = normalizePSD(walkMat.TF);

				rest = filterTimeFreqMat(rest,Channels,chType,chOrder);
				pros = filterTimeFreqMat(pros,Channels,chType,chOrder);
				walk = filterTimeFreqMat(walk,Channels,chType,chOrder);
				
				% increase the average
				[prosAvg prosSqSum] = computeStat(pros,prosAvg,prosSqSum);
				[restAvg restSqSum] = computeStat(rest,restAvg,restSqSum);
				[walkAvg walkSqSum] = computeStat(walk,walkAvg,walkSqSum);

				xAxisBounds = [5 60];

				if sProcess.options.relativeChange.Value

					pros = (pros - rest) ./ rest;
					walk = (walk - rest) ./ rest;

					subplot(3,6,sPlotOrder{subjIdx}(1));
						plot(restMat.Freqs, pros(1:6,:),'--',...
							restMat.Freqs, walk(1:6,:));
						set(gca,'XMinorTick','on');
						title(strcat(subjName,' L'));
						xlim(xAxisBounds);
						axis square;
						box off

					subplot(3,6,sPlotOrder{subjIdx}(2));
						plot(restMat.Freqs, pros(7:end,:),'--',...
							restMat.Freqs, walk(7:end,:));
						set(gca,'XMinorTick','on');
						title(strcat(subjName,' R'));
						xlim(xAxisBounds);
						axis square;
						box off

				else

					subplot(3,6,sPlotOrder{subjIdx}(1));
						
						plot(restMat.Freqs, pros(1:6,:),...
							restMat.Freqs, walk(1:6,:),'--',...
							restMat.Freqs, rest(1:6,:),'-.');
						set(gca,'XMinorTick','on');
						title(strcat(subjName,' L'));
						xlim(xAxisBounds);
						axis square;
						box off

					subplot(3,6,sPlotOrder{subjIdx}(2));
					plot(restMat.Freqs, pros(7:end,:),...
						restMat.Freqs, walk(7:end,:),'--',...
						restMat.Freqs, rest(7:end,:),'-.');

						set(gca,'XMinorTick','on');
						title(strcat(subjName,' R'));
						xlim(xAxisBounds);
						axis square;
						box off

				end
		else
			bst_report('Error',sProcess,sInputs,'Wrong data pairs');
		end

	end; clear el;

	prosStd = sqrt( sum(prosAvg,1).^2 + sum(prosSqSum,1) ) ./(nFiles*2);
	prosAvg = sum(prosAvg,1)./(nFiles*2);

	restStd = sqrt( sum(restAvg,1).^2 + sum(restSqSum,1) ) ./(nFiles*2);
	restAvg = sum(restAvg,1)./(nFiles*2);

	walkStd = sqrt( sum(walkAvg,1).^2 + sum(walkSqSum,1) ) ./(nFiles*2);
	walkAvg = sum(walkAvg,1)./(nFiles*2);



	figure,
	subplot(2,1,1),
	plot(restMat.Freqs,restAvg,'r',...
		restMat.Freqs,walkAvg,'g');
	fill_between(restMat.Freqs,restAvg-(restStd./sqrt(nFiles*2)),...
		restAvg+(restStd./sqrt(nFiles*2)),[],...
			'EdgeColor',[1 0 0],'FaceColor','none');
	fill_between(walkMat.Freqs,walkAvg-(walkStd./sqrt(nFiles*2)),...
		walkAvg+(walkStd./sqrt(nFiles*2)),[],...
			'EdgeColor',[0 1 0],'FaceColor','none');


	set(gca,'XMinorTick','on');
	title('Population Average');
	legend([{'rest'},{'walk'}]);
	xlim(xAxisBounds);
	box off

	subplot(2,1,2)
	plot(restMat.Freqs,prosAvg,...
		restMat.Freqs,restAvg,'r');
	fill_between(restMat.Freqs,restAvg-(restStd./sqrt(nFiles*2)),...
		restAvg+(restStd./sqrt(nFiles*2)),[],...
			'EdgeColor',[1 0 0],'FaceColor','none');
	fill_between(prosMat.Freqs,prosAvg-(prosStd./sqrt(nFiles*2)),....
		prosAvg+(prosStd./sqrt(nFiles*2)),[],...
			'EdgeColor',[0 0 1],'FaceColor','none');



	legend([{'prosup'},{'rest'}]);
	set(gca,'XMinorTick','on');
	title('Population Average');
	xlim(xAxisBounds);
	box off


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

end

function normPSD = normalizePSD(normPSD)
%NORMALIZEPSD Description
%	NORMPSD = NORMALIZEPSD(PSD) Long description
%
	normPSD = squeeze(normPSD);
	normPSD = normPSD ./ repmat(sum(normPSD,2),[1,size(normPSD,2)]);
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
