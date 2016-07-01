function varargout = process_gaitIniImportToDB( varargin )
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
    sProcess.Comment     = 'GaitIni - import recordings';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Import recordings';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
    % === OPTION EXAMPLE
    sProcess.options.example1.Comment = {'Choice 1', 'Choice 2', 'Choice 3'};
    sProcess.options.example1.Type    = 'radio';
    sProcess.options.example1.Value   = 1;
    % === OPTION EXAMPLE
    sProcess.options.example2.Comment = 'Example option 2';
    sProcess.options.example2.Type    = 'checkbox';
    sProcess.options.example2.Value   = 1;
    % === OPTION EXAMPLE
    sProcess.options.example3.Comment = 'Example option 3: ';
    sProcess.options.example3.Type    = 'value';
    sProcess.options.example3.Value   = {5, 'units', 2};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
		
		basePath 	= '/media/gabri/My Passport/gaitIni';
		fname			= fullfile(basePath,'leadingLimb.csv');
		[subjectName, Trial, Drug, LeadingLimb] = textread(fname,'%s %d %s %s\n','delimiter',',');

		ControLat = zeros(numel(sInputs),2);
		IpsiLat = zeros(numel(sInputs),2);
    
    % ===== LOAD THE DATA =====
    for i = 1:length(sInputs)

				fprintf('Processing %d\n',i);

        data = in_bst(sInputs(i).FileName, [], 0);
				eventInfo = in_bst(sInputs(i).DataFile,'Events','Comment');
				trialIndex = str2double(regexp(eventInfo.Comment,'\d+','match'));

				subjectIdx = ismember(subjectName,sInputs(i).SubjectName);
				trialIdx = (Trial == trialIndex);
				drugIdx = ismember(Drug,sInputs(i).Condition);
				leadingLimb = LeadingLimb(logical(subjectIdx .* trialIdx .* drugIdx));

				chLabels = data.RowNames;
				tmp = (regexp(chLabels,'\d+','match'));
				tmp = [tmp{:}];
				tmp = reshape(tmp,[2 2]);
				tmp = sum(cellfun(@str2num,tmp));

				chSide( tmp >= 6 ) = {'sx'};
				chSide( tmp < 6 ) = {'dx'};

				betaBand = data.Freqs >= 13 & data.Freqs <= 30;

				for ch = 1:2

						[val,idx] = findpeaks(squeeze(sum(data.TF(ch,:,betaBand),2)));
						if isempty(val)
								betaPeaks(ch) = 5;
								
						elseif(numel(val) == 1)
								betaPeaks(ch) = idx;
						else
								[v,ii] = max(val);

								betaPeaks(ch) = idx(ii);
						end
				end

				betaPeakFrequency = data.Freqs(betaBand);
				betaPeakFrequency = betaPeakFrequency(betaPeaks);

				betaBands = bsxfun(@plus,betaPeakFrequency', (-5:5));

				for ch = 1:2
						waveletInterp(ch,:,:) = interp1(data.Freqs,squeeze(data.TF(ch,:,:))',betaBands(ch,:))';
				end

				% da t_fine_standing a t_heel_off_leading_limb
				APAonsetWindow 			= eventInfo.Events(1).samples:eventInfo.Events(2).samples;
				GaitIniWindow 			= eventInfo.Events(2).samples:eventInfo.Events(3).samples;
				StandingWindow 			= eventInfo.Events(1).samples-2*numel(APAonsetWindow)...
						:eventInfo.Events(1).samples-numel(APAonsetWindow)-1;

				% we should check for bad segments

				betaStandingWindow 	= mean(mean(waveletInterp(:,StandingWindow,:),2),3);
				betaAPAonset 				= (mean(mean(waveletInterp(:,APAonsetWindow,:),2),3) - ...
						betaStandingWindow)./betaStandingWindow;
				betaGaitIniWindow 	= (mean(mean(waveletInterp(:,GaitIniWindow,:),2),3) - ...
						betaStandingWindow)./betaStandingWindow;

				if strcmp(leadingLimb,'sx')

						ControLat(i,1) 	= betaAPAonset(strcmp(chSide, 'dx'));
						ControLat(i,2)  = betaGaitIniWindow(strcmp(chSide, 'dx'));
						IpsiLat(i,1) 		= betaAPAonset(strcmp(chSide, 'sx'));
						IpsiLat(i,2) 		= betaGaitIniWindow(strcmp(chSide, 'sx'));

				else
						ControLat(i,1) 	= betaAPAonset(strcmp(chSide, 'sx'));
						ControLat(i,2) 	= betaGaitIniWindow(strcmp(chSide, 'sx'));
						IpsiLat(i,1) 		= betaAPAonset(strcmp(chSide, 'dx'));
						IpsiLat(i,2)  	= betaGaitIniWindow(strcmp(chSide, 'dx'));

				end
		clear chSide waveletInterp betaPeaks
    end
		
		meanControlat = squeeze(mean(ControLat,1));
		meanIpsiLat = squeeze(mean(IpsiLat,1));
		stdControlat = squeeze(std(ControLat,1))/sqrt(numel(sInputs));
		stdIpsiLat = squeeze(std(IpsiLat,1))/sqrt(numel(sInputs));

		figure, bar([meanControlat',meanIpsiLat']);
		colormap([0 1 1;1 0 1]);
		title(sInputs(i).Condition);
		hold on
		legend([{'ControLat'},{'IpsiLat'}]);
		errorbar(-.15+(1:2),meanControlat,stdControlat,'LineStyle','none','Color','k');
		errorbar(+.15+(1:2),meanIpsiLat,stdIpsiLat,'LineStyle','none','Color','k');
		set(gca,'XTickLabel',[{'APA vs Stand'},{'GaitIni vs Stand'}])
		ylabel('Relative Change');
		
end
