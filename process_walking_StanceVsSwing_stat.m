function varargout = process_walking_timeWarp2( varargin )
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
    sProcess.Comment     = 'Stance vs Swing Stat';
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

		sProcess.options.saveOutput.Comment = 'Save output to brainstormDB';
		sProcess.options.saveOutput.Type    = 'checkbox';
		sProcess.options.saveOutput.Value   = false;

		sProcess.options.method.Comment = 'Method:';
		sProcess.options.method.Type    = 'text';
		sProcess.options.method.Value   = 'zscore';

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

  trialRejectionFile = fullfile('/home/lgabri/Desktop','walking','trialRejection.csv');
  [patNames,trialStrings,stepIds] = textread(trialRejectionFile,...
																				'%s %*s %*s %*s%s%d%*s','delimiter',',');

	sideFile = fullfile('/home/lgabri/Desktop','walking','patientSides.csv');
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

	swingTAvg = cell(nSubjects,2);
	stanceTAvg = cell(nSubjects,2);

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
			gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'[heel|toe]'));

			% concat all heel contact events in order to have
			% a vector of latencies of this form: e.g.
			% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
			[strideStart,ord]		= sort([parentData.Events(gaitEventGroups).samples]);

			% extract event names
			for evIdx = find(gaitEventGroups)
			
					evLabels{:,evIdx} = repmat({parentData.Events(evIdx).label},...
							[1 numel(parentData.Events(evIdx).samples)]);		

			end

			% re-order event names accordingly
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
			nRowsInPlot		= nStrideLeft+nStrideRight+1;
			plotOrder 		= mat2cell(reshape(1:nRowsInPlot*4,4,nRowsInPlot)',ones(nRowsInPlot,1),[2 2]);

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

			for strideIdx = 3:2:nStrides*2+2


					% if the stride contains bad steps 
					% we skip it and continue to the next
					if ismember(plotIdx,strideRej)
							plotIdx = plotIdx + 1;
							continue;
					end

					%								[ idx ] 	
					% (hc_L) *tof_R  hc_R   *tof_L (hc_L) 
					%    ^ start-2		|t0				      ^ start + 2 == end stride
					timeWindow = strideStart(strideIdx) + (-399:400);
					freqMask 	 = walkingStruct.Freqs > 6;
					dataTF 		 = walkingStruct.TF(:,timeWindow,freqMask);


					strideRaw	 = signals(:,timeWindow)'.*1e6;
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

							zLimit 		= [-2 2];
							tEvAxis 	= repmat(referenceTimeVector(referenceStance([2 3 4])),2, 1);

							% tf_* -> hc_*
							stanceWnd	= referenceStance(2):referenceStance(3);
							% tf_* -> hc_*
							swingWnd	= referenceStance(4):referenceStance(5);
							tAxis			= referenceTimeVector;

					else
							finalTF 	= dataTF.*1e12;
							finalRaw	= strideRaw';
							zLimit 		= [min(finalTF(:)) max(finalTF(:))];
							tAxis			= (-399:400)/fs;
							tEvAxis 	= repmat(originalVector([3 4 5])./400,2, 1);
					end

					footLabel 		 	= regexp(evNames(strideIdx),'[L|R]','match');
					footLabel				= footLabel{:}{:}; % this is the label of the central hc event

					plotOrderStride	= plotOrder(plotIdx,:);

					% we need to switch foot Label since a left foot stride
					% will have a hcR event as central (t0) event
					% stnIdx are ordered as STN- first and STN+ next
					if strcmp(footLabel,'L')
							% right foot stride
							footLabel 	 = 'R';
							controLatIdx = 2;
							if strcmp(mostAffSide,'L')
								plotOrderStride	= plotOrderStride{1};
								stnIdx					= 1;
							else 
								plotOrderStride = plotOrderStride{2};
								stnIdx					= 2;
							end
					else 

							% left foot stride
							footLabel 	 = 'L';
							controLatIdx = 1;
							if strcmp(mostAffSide,'L')
								plotOrderStride	= plotOrderStride{2};
								stnIdx					= 2;
							else 
								plotOrderStride = plotOrderStride{1};
								stnIdx					=1;
							end
					end

					if sProcess.options.normalizeOnStride.Value
							normFactor 	= repmat(mean(finalTF,2),[1 numel(timeWindow) 1]);
							A 					= finalTF./normFactor;
							swingData 	= A(controLatIdx,swingWnd,:);
							stanceData	= A(controLatIdx,stanceWnd,:);
					else
							swingData		= finalTF(controLatIdx,swingWnd,:);
							stanceData	= finalTF(controLatIdx,stanceWnd,:)

					end

					% for each subject and stn-(1)/+(2) we save the time-averaged swing and stance power
					% and save them in separate variables
					swing{subjectIdx,stnIdx} = cat(1,swing{subjectIdx,stnIdx},swingData);
					stance{subjectIdx,stnIdx} = cat(1,stance{subjectIdx,stnIdx},stanceData);
					plotIdx = plotIdx + 1;
					clear finalTF;

			end % for stride

			clear finalTF;

	end % for sInputs files


	f2 = figure('papertype','a4','paperposition',[0 0 1 1],'paperunits','normalized',...
										'paperorientation','portrait','Visible','on');
	
	folder = 'original';
	pvalue = ones(numel(f),2);

	for ii = 1:nSubjects
			
			% nStrides,time(161),f(84)
			swingData  = cat(4,swing{ii,:});
			stanceData = cat(4,stance{ii,:});


			% for a given stride/stance I randomly pick a swing phase
			% compute the zscore averaged across time
			% repeat this estimate with replacement k times
			% then compute the CI for p<0.05

			for surrIdx = 1:nSurrogates

					swingIndices = randperm(1:size(swingData,1));
					swingDataSurr= swingData(swingIndices,:,:,:);
					
					% f x STN - 84 x 2 
					zScore(surrIdx,:,:) = squeeze(mean(mean(bsxfun(@rdivide,bsxfun(@minus,stanceData,...
													mean(swingDataSurr,2)),std(swingDataSurr,[],2))),2));


			end

			% for each STN-/+ we compute a two tailed permutation test
			confidenceIntervals = prctile(zScore,[0.025 0.975]);


			% this is the STN- 
			subplot(nSubjects,4,4*(ii-1)+1)
			imagesc(tAxis,f,squeeze(stnMeans{ii,1})',zLimit);
			axis xy;
			set(gca,'XTickLabel',[]);
			xlim([min(tAxis) max(tAxis)]);
			ylim([6 80]);

			subplot(nSubjects,4,4*(ii-1)+2,'NextPlot','add')
			plot(squeeze( mean(stnMeans{ii,1},2)),f)
			plot([-3 3;-3 3],[0 0;90 90],'k--');	
			plot(significanceMask(:,1),f);
			ylim([6 80]);
			xlim([-5 5]);

			subplot(nSubjects,4,4*(ii-1)+3)
			imagesc(tAxis,f,squeeze(stnMeans{ii,2})',zLimit);
			axis xy;
			set(gca,'XTickLabel',[]);
			xlim([min(tAxis) max(tAxis)]);
			ylim([6 80]);

			subplot(nSubjects,4,4*(ii-1)+4,'NextPlot','add')
			plot(squeeze(mean(stnMeans{ii,2},2)),f)
			plot([-3 3;-3 3],[0 0;90 90],'k--');	
			plot(significanceMask(:,2),f,'ro','MarkerSize',2),2;
			ylim([6 80]);
			xlim([-5 5]);

%			fname = fullfile('/home','lgabri','Desktop','walking',folder,'wavelet',subjectNameOrdered{ii},...
%					'avgZScoreStancevSwing.ps');


	end
			annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
			annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');

			fname = fullfile('/home','lgabri','Desktop','walking',folder,'wavelet',...
					'avgZScoreStancevSwing.ps');

			print(f2,'-dpsc2',fname);
		
end % function
