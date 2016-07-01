function varargout = process_compare_PSDs2( varargin )
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
    sProcess.Comment     = 'Compare PSDs';
    sProcess.FileTag     = '| compared';
    sProcess.Category    = 'Filter2';
    sProcess.SubGroup    = 'Plot';
    sProcess.Index       = 702;

    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq','data'};
    sProcess.OutputTypes = {'timefreq','data'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;

    % Definition of the options
	sProcess.options.relativeChange.Comment = 'Relative Change';
	sProcess.options.relativeChange.Type    = 'checkbox';
	sProcess.options.relativeChange.Value   = 0;
%
%	sProcess.options.compAvg.Comment = 'Compute Group Average';
%	sProcess.options.compAvg.Type    = 'checkbox';
%	sProcess.options.compAvg.Value   = 0;

	sProcess.options.medA.Comment = 'Med(A) ON/OFF';
	sProcess.options.medB.Comment = 'Med(B) ON/OFF';
	sProcess.options.medA.Type = 'text';
	sProcess.options.medB.Type = 'text'; 
	sProcess.options.medA.Value = 'ON';
	sProcess.options.medB.Value = 'OFF'; 

	sProcess.options.condTagA.Comment = 'Cond (A)';
	sProcess.options.condTagB.Comment = 'Cond (B)';...
	sProcess.options.condTagA.Type = 'text';
	sProcess.options.condTagB.Type = 'text'; 

	sProcess.options.condTagA.Value = 'rest';
	sProcess.options.condTagB.Value = 'rest'; 

	sProcess.options.chType.Comment = 'Channel Type: ';
	sProcess.options.chType.Value = 'MISC';
	sProcess.options.chType.Type = 'text';
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
%
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA,sInputsB) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
    
    % ===== LOAD THE DATA =====
    % Read the first file in the list, to initialize the loop
	
	if numel(sInputsA) ~= numel(sInputsB)
		bst_report('Error',sProcess,sInputsA,...
			'sInputsA & sInputsB do not have the same amount of files');
		return;
	end

	nFiles = numel(sInputsA);
	
	aCond = sProcess.options.condTagA.Value;
	bCond = sProcess.options.condTagB.Value;

	aMed  = sProcess.options.medA.Value;
	bMed  = sProcess.options.medB.Value;

	chType = sProcess.options.chType.Value;

	figure(1);
	colorOrder = jet(6);
	set(gcf,'DefaultAxesColorOrder',colorOrder);
	sPlotOrder = reshape(1:18,[2 9])';
	sPlotOrder = mat2cell(sPlotOrder,ones(9,1),2);

	tmp = [2 3 4 5 6 7 9 10 11];
	subjectNames = arrayfun(@(x) sprintf('wue%02d',x),tmp,'uni',false);

%	avgA = zeros(2,129);
%	avgB = zeros(2,129);
%
%	sqSumA = zeros(2,129);
%	sqSumB = zeros(2,129);

	for el = 1:nFiles	

		DataMatA = in_bst(sInputsA(el).FileName, [], 0);
		DataMatB = in_bst(sInputsB(el).FileName, [], 0);

		% we should check that we are effectively com...paring 
		% the correct studies
		[parStructA, subjNameA,medA, visitA,condTagA] = extractInfo(sInputsA(el));
		[parStructB, subjNameB,medB, visitB,condTagB] = extractInfo(sInputsB(el));

		% Do these files belongs to the same subject?
		if (strcmp(subjNameB,subjNameA) & ... 
			strcmp(condTagA,aCond) & strcmp(condTagB,bCond) & ...
			strcmp(medA,aMed) & strcmp(medB,bMed) )
				subjIdx = find(strcmp(subjectNames,subjNameA));

				aChannels = in_bst_channel(sInputsA(el).ChannelFile,'Channel');
				bChannels = in_bst_channel(sInputsB(el).ChannelFile,'Channel');

				% the ChannelType is fixed to MISC for those
				% that have been selected for further analyses 
				% after visual inspection of beta peaks among the 
				% significant ones. Thus we consider only the MISC channels
				% and discarding all the SEEG channels remaning
				[aChLabels, aChOrder] = sortLabels({aChannels.Channel.Name});
				[bChLabels, bChOrder] = sortLabels({bChannels.Channel.Name});

				psdA = squeeze(DataMatA.TF);
				psdB = squeeze(DataMatB.TF);

				aMask = strcmp({aChannels.Channel.Type},chType);
				bMask = strcmp({bChannels.Channel.Type},chType);

				psdA(~aMask,:) = NaN;
				psdB(~bMask,:) = NaN;

				psdA = psdA(aChOrder,:);
				psdB = psdB(bChOrder,:);

				psdA = psdA ./ repmat(sum(psdA,2),[1,size(psdA,2)]);
				psdB = psdB ./ repmat(sum(psdB,2),[1,size(psdB,2)]);
%
%				[avgA sqSumA] = computeStat(psdA,avgA,sqSumA);
%				[avgB sqSumB] = computeStat(psdB,avgB,sqSumB);

				xAxisBounds = [5 60];

				if sProcess.options.relativeChange.Value

					relChange = (psdA - psdB) ./ psdB;
					subplot(3,6,sPlotOrder{subjIdx}(1));
						plot(DataMatA.Freqs, relChange(1:6,:));
						set(gca,'XMinorTick','on');
						title(strcat(subjNameA,' L'));
						xlim(xAxisBounds);
						axis square;
						box off

					subplot(3,6,sPlotOrder{subjIdx}(2));
						plot(DataMatA.Freqs, relChange(7:end,:));
						set(gca,'XMinorTick','on');
						title(strcat(subjNameA,' R'));
						xlim(xAxisBounds);
						box off
						axis square;

				else

					subplot(3,6,sPlotOrder{subjIdx}(1));
						semilogy(DataMatA.Freqs, psdA(1:6,:),'--',...
							DataMatB.Freqs, psdB(1:6,:));
						set(gca,'XMinorTick','on');
						title(strcat(subjNameA,' L'));
						xlim(xAxisBounds);
						box off

					subplot(3,6,sPlotOrder{subjIdx}(2));
						semilogy(DataMatA.Freqs, psdA(7:end,:),'--',...
						DataMatB.Freqs, psdB(7:end,:));
						set(gca,'XMinorTick','on');
						title(strcat(subjNameA,' R'));
						xlim(xAxisBounds);
						box off

				end
		else
			bst_report('Error',sProcess,sInputsA,'Wrong data pairs');
		end

	end; clear el;
%
%	stdA = sqrt( sum(avgA,1).^2 + sum(sqSumA,1) ) ./(nFiles*2);
%	avgA = sum(avgA,1)./(nFiles*2);
%
%	stdB = sqrt( sum(avgB,1).^2 + sum(sqSumB,1) ) ./(nFiles*2);
%	avgB = sum(avgB,1)./(nFiles*2);
%
%	figure(2);
%	plot(DataMatA.Freqs,avgA,'r','LineWidth',4);
%	fill_between(DataMatA.Freqs, avgA-stdA,avg+stdA,[],...
%			'EdgeColor',0.5.*[1 0 0],'FaceColor',0.5.*[1 0 0],'LineWidth',2);
%	
%	plot(DataMatB.Freqs,avgB,'r','LineWidth',4);
%	fill_between(DataMatB.Freqs, avgB-stdB,avg+stdB,[],...
%			'EdgeColor',0.5.*[1 0 0],'FaceColor',0.5.*[1 0 0],'LineWidth',2);
%
%	set(gca,'XMinorTick','on');
%	title('Population Average');
%	xlim(xAxisBounds);
%	box off

end




function [parentStruct, subjectName,med, visit,condTag] = extractInfo(in)
%EXTRACTINFO Description
%	[PARENTSTRUCT, SUBJECTNAME,MED, VISIT,CONDTAG] = EXTRACTINFO(IN) Long description
%
	parentStruct = bst_process('GetInputStruct',in.DataFile); 
	subjectName = parentStruct.SubjectName;
	med = regexp(parentStruct.Condition,'ON|OFF','match');
	visit = regexp(parentStruct.Condition,'visite_\d+','match');
	condTag = regexp(parentStruct.Condition,'rest|ProSup|walking','match');
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
%
%function [mask] = filterChannelFile(channel,chType)
%%FILTERCHANNELFILE Description
%%	[CHANNEL] = FILTERCHANNELFILE(CHANNEL,CHTYPE) Long description
%%
%
%	chLabels = {channel.Channel.Name};
%	chLabels = cellfun(@(x) x(1:end-5),chLabels,'UniformOutput',false);
%	[chLabels, chOrder] = sortLabels(chLabels);
%
%
%
%end

function [avg sqsum] = computeStat(psd,avg,sqsum)
%COMPUTESTAT Description
%	[AVG STD] = COMPUTESTAT(PSD,AVG) Long description
%
	
	% remove the nan channels
	psd(sum(isnan(psd),2) > 1,:) = [];

	avg = avg + psd;
	sqsum  =  sqsum + psd.^2;
end
