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

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Montage';
    sProcess.FileTag     = '| extracted';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Plot';
    sProcess.Index       = 1000;
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

	sProcess.options.Axis.Comment = 'X Axis Bounds';
    sProcess.options.Axis.Type = 'text';
	sProcess.options.Axis.Value = '0,100';


end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
	OutputFiles = {};
	
	xAxisBounds = textscan(sProcess.options.Axis.Value,'%f','delimiter',',');
	xAxisBounds = xAxisBounds{:};

	figure(1);
	colorOrder = jet(6);
	set(gcf,'DefaultAxesColorOrder',colorOrder);
	sPlotOrder = reshape(1:18,[2 9])';
	sPlotOrder = mat2cell(sPlotOrder,ones(9,1),2);

	tmp = [2 3 4 5 6 7 9 10 11];
	subjectNames = arrayfun(@(x) sprintf('wue%02d',x),tmp,'uni',false);
	
    % Reading all the input files in a big matrix
    for i = 1:length(sInputs)

		% check whether we are looking at visite_01 rest MedOFF
		parentStruct = bst_process('GetInputStruct',sInputs(i).DataFile); 
		subjectName = parentStruct.SubjectName;
		med = regexp(parentStruct.Condition,'ON|OFF','match');
		visit = regexp(parentStruct.Condition,'visite_\d+','match');
		condTag = regexp(parentStruct.Condition,'rest|ProSup|walking','match');

		if strcmp(med,'OFF') & strcmp(visit,'visite_01') & strcmp(condTag,'rest')
			logFid = fopen(fullfile('/home/gabri/',strcat(subjectName,'.log')),'w');
			subjIdx = find(strcmp(subjectNames,subjectName));

        	dataMat = in_bst_timefreq(sInputs(i).FileName);
			channelMat = in_bst_channel(parentStruct.ChannelFile,'Channel');
			chFlag = in_bst_data(parentStruct.FileName,'ChannelFlag');

			chLabels = {channelMat.Channel.Name};
			chLabels = cellfun(@(x) x(1:end-5),chLabels,'UniformOutput',false);
			[chLabels, chOrder] = sortLabels(chLabels);


			psds = squeeze(dataMat.TF);

			psds(chFlag.ChannelFlag==-1,:) = NaN;

			psds = psds(chOrder,:);

			freqs = dataMat.Freqs;

			if sProcess.options.normalize.Value
				psds = psds./repmat(sum(psds,2),[1 size(psds,2)]);
			end

			fprintf(logFid,'%s\n',strjoin(chLabels,','));
			fprintf(logFid,'\n');
				
			subplot(3,6,sPlotOrder{subjIdx}(1));
				plot(freqs,10.*log10(psds(1:6,:)));
				set(gca,'XMinorTick','on');
				xlabel('Frequency (Hz)');
				ylabel('Power (dB)');
				title(strcat(subjectName,' L'));
				xlim(xAxisBounds);
				ylim([-200 -160]);
				axis square;
				box off
			subplot(3,6,sPlotOrder{subjIdx}(2));
				plot(freqs,10.*log10(psds(7:end,:)));
				set(gca,'XMinorTick','on');
				xlabel('Frequency (Hz)');
				ylabel('Power (dB)');
				title(strcat(subjectName,' R'));
				xlim(xAxisBounds);
				ylim([-200 -160]);

				axis square;
				box off
				fclose(logFid);

			
		else
			continue
		end


    end
	
	fnameOut ='/media/gabri/My Passport/DBS/Fig1.eps';
	print(gcf, '-depsc2',fnameOut)
	fixPSlinestyle(fnameOut,fnameOut);

    
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
