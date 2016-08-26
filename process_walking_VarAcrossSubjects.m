function varargout = process_walking_VarAcrossSubjects( varargin )
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
    sProcess.Comment     = 'Stride Variability';
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
		sProcess.options.sensortypes.Comment = 'Sensor types or names: ';
		sProcess.options.sensortypes.Type    = 'text';
		sProcess.options.sensortypes.Value   = 'SEEG';

		sProcess.options.doWarping.Comment = 'time warping ';
		sProcess.options.doWarping.Type    = 'checkbox';
		sProcess.options.doWarping.Value   = true;
%
%		sProcess.options.saveOutput.Comment = 'Save output to brainstormDB';
%		sProcess.options.saveOutput.Type    = 'checkbox';
%		sProcess.options.saveOutput.Value   = false;
%
%		sProcess.options.method.Comment = 'Method:';
%		sProcess.options.method.Type    = 'text';
%		sProcess.options.method.Value   = 'zscore';
%
		sProcess.options.normalizeOnStride.Comment = 'Normalize On Stride';
		sProcess.options.normalizeOnStride.Type = 'checkbox';
		sProcess.options.normalizeOnStride.Value = false; 

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

	nFiles = numel(sInputs);

	DATA_FOLDER = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','info');
  trialRejectionFile = fullfile(DATA_FOLDER,'trialRejection.csv');
  [patNames,trialStrings,stepIds] = textread(trialRejectionFile,...
																				'%s %*s %*s %*s%s%d%*s','delimiter',',');

	sideFile = fullfile(DATA_FOLDER,'patientSides.csv');
	[subjectNames, mostAffSides] = textread(sideFile,'%s %s\n','delimiter',',');
%	
%	cynematicFile = fullfile('/home/lgabri/Desktop','walking','cynematics.csv');
%	[velocity, lenght] = textread(cynematicFile,'%f %f\n','delimiter',',');

	OutputFiles = {};

	nSubjects = numel(unique([{sInputs.SubjectName}]));
	stnResults = cell(nSubjects,2); % this var holds for each subjects the STN-/+ power for each stride
	currentSubject=[];
	subjectIdx = 0;
	subjectNameOrdered = cell(nSubjects,1);

	swing = cell(nSubjects,2);
	stance = cell(nSubjects,2);

	for fileIdx = 1:nFiles; 


			walkingStruct = in_bst_timefreq(sInputs(fileIdx).FileName); 

			parentStruct 	= bst_process('GetInputStruct',walkingStruct.DataFile);
			parentData 		= in_bst(parentStruct.FileName);
			channelData		= in_bst_channel(sInputs(fileIdx).ChannelFile);
			iChannels			= channel_find(channelData.Channel,'SEEG');
			signals				= parentData.F(iChannels,:);

			mostAffSide		= mostAffSides(strcmpi(subjectNames,(sInputs(fileIdx).SubjectName) ));

			if isempty(currentSubject) || ~strcmp(currentSubject,sInputs(fileIdx).SubjectName)
				currentSubject = sInputs(fileIdx).SubjectName;
				subjectIdx = subjectIdx + 1;
				subjectNameOrdered{subjectIdx} = sInputs(fileIdx).SubjectName;
			end

			fs = round(1/mean(diff( parentData.Time )));

			% filter cardiac events from gait-related events
			evGroupNames = {parentData.Events.label};
			gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel|toe)'));

			% concat all heel contact events in order to have
			% a vector of latencies of this form: e.g.
			% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
			[strideStart,ord]		= sort([parentData.Events(gaitEventGroups).samples]);

			% extract event names
			evLabels = cell(1,sum(gaitEventGroups));
			evLabelsIdx = 1;
			for evIdx = find(gaitEventGroups)
			
					evLabels{evLabelsIdx} = repmat({parentData.Events(evIdx).label},...
							[1 numel(parentData.Events(evIdx).samples)]);		
					evLabelsIdx = evLabelsIdx + 1;

			end			% re-order event names accordingly
			evNames 	= [evLabels{:}];
			evNames 	= evNames(ord);
	
			% count how many strides we have recorded
			nStrideLeft 	= sum(strcmp(evNames,'heelcontact_L'))-1; 
			nStrideRight	= sum(strcmp(evNames,'heelcontact_R'))-1; 
			nStrides 			= nStrideLeft + nStrideRight;

			% prepare event mask to reject artefactual or incomplete step
			trialString 	= regexp(sInputs(fileIdx).Condition,'trial\d+','match');
			subjMask 			= ismember(lower(patNames),lower(sInputs(fileIdx).SubjectName)); 
			trialMask 		= (ismember(lower(trialStrings),lower(trialString)));
			stepRej  			= stepIds(and(subjMask,trialMask));

			% stride mask each stride is composed by two steps
			strideIndexes = [1:nStrides;2:nStrides+1]';
			strideRej 		= find(sum(ismember(strideIndexes,stepRej),2));

			% prepare plot order
%			nRowsInPlot		= nStrideLeft+nStrideRight+1;
%			plotOrder 		= mat2cell(reshape(1:nRowsInPlot*4,4,nRowsInPlot)',ones(nRowsInPlot,1),[2 2]);

			% we have to correct the event adding the offset
			% since they are referred to the 0 of the raw data 
			evOffset 			= round(walkingStruct.Time(1)*fs);
			strideStart 	= strideStart - evOffset;

			% we create the normalized stride time vector
			referenceTimeVector = -1:1/fs:(1-1/fs);
			doubleSuppDur 			= floor(0.19*fs);
			singleSuppDur 			= floor(0.4*fs);

			referenceStance 		= 400 + [ -doubleSuppDur-singleSuppDur -singleSuppDur 0 ...
																			+doubleSuppDur +doubleSuppDur+singleSuppDur ];

			referenceVector			= [1 referenceStance 800];

			plotIdx = 1;

			strideCheck = evNames(bsxfun(@plus,(-2:2),(3:2:numel(evNames)-2)'));
			nStrides = size(strideCheck,1);
			matchingString = {'heelcontact_[L|R]', 'toeoff_[R|L]',...
						'heelcontact_[R|L]','toeoff_[L|R]','heelcontact_[L|R]'};
			orderCheck = nan(1,nStrides);

			for el = 1:nStrides
					orderCheck(el) = sum(cellfun(@isempty,cellfun(@regexp,strideCheck(el,:),...
															matchingString,'uni',false)));
			end

			for strideIdx = 3:2:nStrides*2+2


					% if the stride contains bad steps 
					% we skip it and continue to the next
					if ismember(plotIdx,strideRej) || orderCheck(plotIdx) > 0
							plotIdx = plotIdx + 1;
							continue;
					end

					%								[ idx ] 	
					% (hc_L) *tof_R  hc_R   *tof_L (hc_L) 
					%    ^ start-2		|t0				      ^ start + 2 == end stride
					timeWindow = strideStart(strideIdx) + (-399:400);
					freqMask 	 = walkingStruct.Freqs > 6;
					dataTF 		 = walkingStruct.TF(:,timeWindow,freqMask);

					strideRaw	 = signals(:,timeWindow)';
					f 				 = walkingStruct.Freqs(freqMask);	

					% then create the time-warping vector
					originalVector = [1 (strideStart((-2:1:2) + strideIdx) - timeWindow(1))  800];

					if sProcess.options.doWarping.Value

							% compute the mixing matrix that maps the single orignal
							% stance on the normalized stance
							mixingMatrix 	 = mytimewarp(referenceVector,originalVector,3);

							% apply warping at each channel separately for both TF 
							finalTF(1,:,:) = mixingMatrix * squeeze(dataTF(1,:,:));
							finalTF(2,:,:) = mixingMatrix * squeeze(dataTF(2,:,:));
					
							% and raw data
							finalRaw(1,:)	 = mixingMatrix * strideRaw(:,1);
							finalRaw(2,:)	 = mixingMatrix * strideRaw(:,2);

%							zLimit 		= [-2 2];
%							tEvAxis 	= repmat(referenceTimeVector(referenceStance([2 3 4])),2, 1);
%
							% tf_* -> hc_*
							stanceWnd	= referenceStance(2):referenceStance(3);
							% tf_* -> hc_*
							swingWnd	= referenceStance(4):referenceStance(5);
%							tAxis			= referenceTimeVector;

					else
							finalTF 	= dataTF;
							finalRaw	= strideRaw';
%							zLimit 		= [min(finalTF(:)) max(finalTF(:))];
%							tAxis			= (-399:400)/fs;
%							tEvAxis 	= repmat(originalVector([3 4 5])./400,2, 1);
					end

					footLabel 		 	= regexp(evNames(strideIdx),'[L|R]','match');
					footLabel				= footLabel{:}{:}; % this is the label of the central hc event


					% we need to switch foot Label since a left foot stride
					% will have a hcR event as central (t0) event
					% stnIdx are ordered as STN- first and STN+ next
					if strcmp(footLabel,'L')
							% right foot stride
							controLatIdx = 2;
							if strcmp(mostAffSide,'L')
								stnIdx					= 1;
							else 
								stnIdx					= 2;
							end
					else 

							% left foot stride
							controLatIdx = 1;
							if strcmp(mostAffSide,'L')
								stnIdx					= 2;
							else 
								stnIdx					=1;
							end
					end

					if sProcess.options.normalizeOnStride.Value
							normFactor 	= repmat(mean(finalTF,2),[1 numel(timeWindow) 1]);
							A 					= finalTF./normFactor;
							swingData 	= A(controLatIdx,swingWnd,:);
							stanceData	= A(controLatIdx,stanceWnd,:);

					else

							swingData 	= finalTF(controLatIdx,swingWnd,:);
							stanceData	= finalTF(controLatIdx,stanceWnd,:);

					end

					% for each subject and stn-(1)/+(2) we save the time-averaged swing and stance power
					swing{subjectIdx,stnIdx} = cat(1,swing{subjectIdx,stnIdx},swingData);
					stance{subjectIdx,stnIdx} = cat(1,stance{subjectIdx,stnIdx},stanceData);

					plotIdx = plotIdx + 1;
					clear finalTF;

			end % for stride

			clear finalTF;

	end % for sInputs files


	stanceTAvg 		= cellfun(@mean,stance,repmat({2},size(stance)),'UniformOutput',false);
	swingTAvg  		= cellfun(@mean,swing,repmat({2},size(swing)),'UniformOutput',false);
	stanceTAvg		= cellfun(@squeeze,stanceTAvg,'UniformOutput',false);
	swingTAvg		  = cellfun(@squeeze,swingTAvg,'UniformOutput',false);
	stanceStrTAvg = cellfun(@mean,stanceTAvg,'UniformOutput',false);
	swingStrTAvg  = cellfun(@mean,swingTAvg,'UniformOutput',false);
	stanceStrTStd = cellfun(@std,stanceTAvg,'UniformOutput',false);
	swingStrTStd  = cellfun(@std,swingTAvg,'UniformOutput',false);
	stanceStrTAvg = cellfun(@squeeze,stanceStrTAvg,'UniformOutput',false);
	swingStrTAvg  = cellfun(@squeeze,swingStrTAvg,'UniformOutput',false);
	stanceStrTStd = cellfun(@squeeze,stanceStrTStd,'UniformOutput',false);
	swingStrTStd  = cellfun(@squeeze,swingStrTStd,'UniformOutput',false);

	
	f2 = figure('papertype','a4','paperposition',[0 0 1 1],'paperunits','normalized',...
										'paperorientation','portrait','Visible','on');
	

	for ii = 1:nSubjects
			pvalue = ones(numel(f),2);

			for freqPoint = 1:numel(f)
				% for each STN-/+ and for each frequency we compute the wilcoxon ranksum test
				% to compare stance and swing periods in terms of beta power
				pvalue(freqPoint,1) = ranksum(stanceTAvg{ii,1}(:,freqPoint),...
																			swingTAvg{ii,1}(:,freqPoint));
				pvalue(freqPoint,2) = ranksum(stanceTAvg{ii,2}(:,freqPoint),...
																			swingTAvg{ii,2}(:,freqPoint));

			end

			significanceBar = nan(numel(f),2);
			significanceBar(pvalue < 0.05) = 1e-18;

			% this is the STN- 
			subplot(nSubjects,2,2*(ii-1)+1,'NextPlot','add')
			plot(f,stanceStrTAvg{ii,1},f,swingStrTAvg{ii,1},'LineWidth',2)
			plot(f,stanceStrTAvg{ii,1}+stanceStrTStd{ii,1},'b--',...
						f,stanceStrTAvg{ii,1}-stanceStrTStd{ii,1},'b--',...
						f,swingStrTAvg{ii,1}+swingStrTStd{ii,1},'g--',...
						f,swingStrTAvg{ii,1}-swingStrTStd{ii,1},'g--');

%			plot(f,significanceBar(:,1),'rx','MarkerSize',4);
			yLimValues = get(gca,'YLim');
			ylim([0 yLimValues(2)]);
			xlim([min(f) 40]);
			xlabel('Frequency (Hz)');


			subplot(nSubjects,2,2*(ii-1)+2,'NextPlot','add')
			plot(f,stanceStrTAvg{ii,2},f,swingStrTAvg{ii,2},'LineWidth',2)
			plot(f,stanceStrTAvg{ii,2}+stanceStrTStd{ii,2},'b--',...
						f,stanceStrTAvg{ii,2}-stanceStrTStd{ii,2},'b--',...
						f,swingStrTAvg{ii,2}+swingStrTStd{ii,2},'g--',...
						f,swingStrTAvg{ii,2}-swingStrTStd{ii,2},'g--');

			yLimValues = get(gca,'YLim');
			ylim([0 yLimValues(2)]);
			
			xlim([min(f) 40]);
%			plot(f,significanceBar(:,2),'rx','MarkerSize',4);
			xlabel('Frequency (Hz)');

	end

	annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
	annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');
	fname = fullfile('/','home','lgabri','Dropbox','Isaias_group','walking','figs',...
			'powerVariability.ps');
	print(f2,'-dpsc2',fname);
		
end % function
