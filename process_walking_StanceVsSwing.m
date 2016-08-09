function varargout = process_walking_StanceVsSwing( varargin )
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
    sProcess.Comment     = 'Stance vs Swing';
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


		sProcess.options.normalizeOnStride.Comment = 'Normalize On Stride';
		sProcess.options.normalizeOnStride.Type = 'checkbox';
		sProcess.options.normalizeOnStride.Value = true; 

		sProcess.options.method.Comment = 'Method:';
		sProcess.options.method.Type    = 'text';
		sProcess.options.method.Value   = 'zscore';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

	nFiles = numel(sInputs);

  trialRejectionFile = fullfile('/media/gabri/My Passport','trialRejection.csv');
  [patNames,trialStrings,stepIds] = textread(trialRejectionFile,...
																				'%s %*s %*s %*s%s%d%*s','delimiter',',');

	sideFile = fullfile('/media/gabri/My Passport','patientSides.csv');
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

	swingTcourse = cell(nSubjects,2);
	stanceTcourse = cell(nSubjects,2);

	for fileIdx = 1:nFiles; 


			walkingStruct = in_bst_timefreq(sInputs(fileIdx).FileName); 

			parentStruct 	= bst_process('GetInputStruct',walkingStruct.DataFile);
			parentData 		= in_bst(parentStruct.FileName);
			channelData		= in_bst_channel(sInputs(fileIdx).ChannelFile);
			iChannels			= channel_find(channelData.Channel,'SEEG');
			signals				= parentData.F(iChannels,:);

			mostAffSide		= mostAffSides(strcmpi(subjectNames,...
														(sInputs(fileIdx).SubjectName) ));

			if isempty(currentSubject) || ...
							~strcmp(currentSubject,sInputs(fileIdx).SubjectName)

				currentSubject = sInputs(fileIdx).SubjectName;
				subjectIdx = subjectIdx + 1;
				subjectNameOrdered{subjectIdx} = sInputs(fileIdx).SubjectName;
			end

			fprintf('Analyzing %s\n',sInputs(fileIdx).SubjectName);

			fs = round(1/mean(diff( parentData.Time )));

			% filter cardiac events from gait-related events
			evGroupNames = {parentData.Events.label};
			gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel|toe)'));

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
			nStrideLeft = sum(strcmp(evNames,'heelcontact_L'))-1; 
			nStrideRight= sum(strcmp(evNames,'heelcontact_R'))-1; 
			nStrides 		= nStrideLeft + nStrideRight;

			% prepare event mask to reject artefactual or incomplete step
			trialString = regexp(sInputs(fileIdx).Condition,'trial\d+','match');
			subjMask 		= ismember(lower(patNames),lower(sInputs(fileIdx).SubjectName)); 
			trialMask 	= (ismember(lower(trialStrings),lower(trialString)));
			stepRej  		= stepIds(and(subjMask,trialMask));

			% stride mask each stride is composed by two steps
			strideIndexes = [1:nStrides;2:nStrides+1]';
			strideRej 		= find(sum(ismember(strideIndexes,stepRej),2));

			% prepare plot order
			nRowsInPlot		= nStrideLeft+nStrideRight+1;
			plotOrder 		= mat2cell(reshape(1:nRowsInPlot*4,4,nRowsInPlot)',...
																		ones(nRowsInPlot,1),[2 2]);

			% we have to correct the event adding the offset
			% since they are referred to the 0 of the raw data 
			evOffset 			= round(walkingStruct.Time(1)*fs);
			strideStart 	= strideStart - evOffset;

			% we create the normalized stride time vector
			referenceTimeVector = -1:1/fs:(1-1/fs);
			doubleSuppDur 			= floor(0.19*fs);
			singleSuppDur 			= floor(0.4*fs);

			referenceStance 		= 400 +[ -doubleSuppDur-singleSuppDur -singleSuppDur 0 ...
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

					footLabel 		 	= regexp(evNames(strideIdx),'[L|R]','match');
					footLabel				= footLabel{:}{:}; % this is the label of the central hc 

%					plotOrderStride	= plotOrder(plotIdx,:);

					% we need to switch foot Label since a left foot stride
					% will have a hcR event as central (t0) event
					% stnIdx are ordered as STN- first and STN+ next
					if strcmp(footLabel,'L')
							% right foot stride
							footLabel 	 = 'R';
							controLatIdx = 2;
							if strcmp(mostAffSide,'L')
%								plotOrderStride	= plotOrderStride{1};
								stnIdx					= 1;
							else 
%								plotOrderStride = plotOrderStride{2};
								stnIdx					= 2;
							end
					else 

							% left foot stride
							footLabel 	 = 'L';
							controLatIdx = 1;
							if strcmp(mostAffSide,'L')
%								plotOrderStride	= plotOrderStride{2};
								stnIdx					= 2;
							else 
%								plotOrderStride = plotOrderStride{1};
								stnIdx					= 1;
							end
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
					originalVector = [1 (strideStart((-2:1:2) + strideIdx) ...
																	- timeWindow(1))  800];

					if sProcess.options.doWarping.Value

							% compute the mixing matrix that maps the single orignal
							% stride on the normalized stride
							mixingMatrix 	 = mytimewarp(referenceVector,originalVector,3);

							% apply warping at each channel separately for both TF 
							finalTF(1,:,:) = mixingMatrix * squeeze(dataTF(1,:,:));
							finalTF(2,:,:) = mixingMatrix * squeeze(dataTF(2,:,:));
					
							% and raw data
							finalRaw(1,:)	 = mixingMatrix * strideRaw(:,1);
							finalRaw(2,:)	 = mixingMatrix * strideRaw(:,2);

							zLimit 		= [-2 2];
							tEvAxis 	= repmat(referenceTimeVector(...
																	referenceStance([2 3 4])),2, 1);

							% tf_* -> hc_*
							stanceWnd	= referenceStance(2):referenceStance(3);
							% tf_* -> hc_*
							swingWnd	= referenceStance(4):referenceStance(5);
							tAxis			= referenceTimeVector;

					else
							finalTF 	= dataTF.*1e12;
							finalRaw	= strideRaw';
							zLimit 		= [min(finalTF(:)) max(finalTF(:))];
							% tf_* -> hc_*
							stanceWnd	= originalVector(3):originalVector(4);
							% tf_* -> hc_*
							swingWnd	= originalVector(5):originalVector(6);

							tAxis			= (-399:400)/fs;
							tEvAxis 	= repmat(originalVector([3 4 5])./400,2, 1);
					end
					% here depending on the options chosen we might 
					% end up having a time-warped data with fixed stride lengths
					% or individual stride lenghts with their original data
					% thus we need to account for different sizes

					if sProcess.options.normalizeOnStride.Value
							normFactor 	= repmat(mean(finalTF,2),[1 numel(timeWindow) 1]);
							A 					= finalTF./normFactor;
							swingData 	= A(controLatIdx,swingWnd,:);
							stanceData	= A(controLatIdx,stanceWnd,:);

					else
							swingData = finalTF(controLatIdx,swingWnd,:);
							stanceData = finalTF(controLatIdx,stanceWnd,:);

					end

					if strcmp(sProcess.options.method.Value,'zscore')
						% z score swing									
						swingZScored = bsxfun(@rdivide,bsxfun(@minus,finalTF(:,swingWnd,:),...
															mean(finalTF(:,swingWnd,:),2)),...
															std(finalTF(:,swingWnd,:),[],2));
						
						% z score stance
						stanceZScored = bsxfun(@rdivide,bsxfun(@minus,finalTF(:,stanceWnd,:),...
																			mean(finalTF(:,stanceWnd,:),2)),...
																			std(finalTF(:,stanceWnd,:),[],2));


						finalTFZScored = bsxfun(@rdivide,...
																bsxfun(@minus,finalTF(:,stanceWnd,:),...
																mean(finalTF(:,swingWnd,:),2)),...
																std(finalTF(:,swingWnd,:),[],2));

					elseif(strcmp(sProcess.options.method.Value,'rchange'))
						% relative change 
						finalTFZScored = bsxfun(@rdivide,...
																bsxfun(@minus,finalTF(:,stanceWnd,:),...
																mean(finalTF(:,swingWnd,:),2)),...
																mean(finalTF(:,swingWnd,:),2));

					else
						finalTFZScored = stanceZScored - swingZScored;

					end

					if sProcess.options.doWarping.Value
							stnResults{subjectIdx,stnIdx} = ...
																cat(1,stnResults{subjectIdx,stnIdx},...
																finalTFZScored(controLatIdx,:,:));

							% for each subject and stn-(1)/+(2) we save the time-averaged 
							% swing and stance power
							swingTcourse{subjectIdx,stnIdx} = ...
									cat(1,swingTcourse{subjectIdx,stnIdx},swingData);
							stanceTcourse{subjectIdx,stnIdx} = ...
									cat(1,stanceTcourse{subjectIdx,stnIdx},stanceData);
					else
							stnResults{subjectIdx,stnIdx} = ...
										cat(1,stnResults{subjectIdx,stnIdx},...
										squeeze(mean(finalTFZScored(controLatIdx,:,:),2))');
					end
%
%					finalTF = finalTF(:,referenceStance(1):referenceStance(end),:);
%					time 		= referenceTimeVector(referenceStance(1):referenceStance(end));
%					
%					if sProcess.options.saveOutput.Value
%							% save each stance in separate file
%							% get the output study
%							iStudy 						= sInputs(fileIdx).iStudy;
%							DataMat 					= walkingStruct;
%							DataMat.Freqs			= f;
%							DataMat.Time			= time;
%							DataMat.TF        = finalTF;
%							DataMat.DataType  = 'data';
%							DataMat.Comment		= sprintf('Stride %s (#%d)',footLabel,plotIdx);
%										
%							% Create a default output filename 
%							OutputFiles{fileIdx} = bst_process('GetNewFilename', ...
%									fileparts(sInputs(fileIdx).FileName), 'timefreq');
%
%							% Save on disk
%							save(OutputFiles{fileIdx}, '-struct', 'DataMat');
%
%							% Register in database
%							db_add_data(iStudy, OutputFiles{fileIdx}, DataMat);
%
%					end
					plotIdx = plotIdx + 1;
					clear finalTF;

			end % for stride

			clear finalTF;
			clear strideStart;

	end % for sInputs files

	% stnResults holds the zscored-power modulation (1,t,f) for stance {:,:,1} 
	% and swing {:,:,2} for each subjects {ii,:,:} and each strides and for 
	% both STN- {:,1,:} and STN+ {:,2,:}.
	% We first average each element of stnResults across time in order to have 
	% a frequency modulation
	stnMeans	 = cellfun(@mean,stnResults,'UniformOutput',false);
	stanceTAvg = cellfun(@mean,stanceTcourse,...
									repmat({2},size(stanceTcourse)),'UniformOutput',false);
	swingTAvg = cellfun(@mean,swingTcourse,...
									repmat({2},size(swingTcourse)),'UniformOutput',false);

	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
			'paperunits','normalized','paperorientation','portrait','Visible','on');
	
	folder = 'original';
	pvalue = zeros(numel(f),2);

	for ii = 1:nSubjects

			pvalue(:,1) = runPermutationTest(stnMeans{ii,1},...
					stanceTcourse{ii,1},swingTcourse{ii,1},100);
			pvalue(:,2) = runPermutationTest(stnMeans{ii,2},...
					stanceTcourse{ii,2},swingTcourse{ii,2},100);

			
			% need multiple comparison correction Bonferroni or FDR
			significanceMask = nan(size(pvalue));
			significanceMask(pvalue < 0.05) = 4;

			% this is the STN- 
			subplot(nSubjects,4,4*(ii-1)+1)
			imagesc(tAxis,f,squeeze(stnMeans{ii,1})',zLimit);
			axis xy;
			set(gca,'XTickLabel',[]);
			xlim([min(tAxis) max(tAxis)]);
			ylim([6 80]);

			subplot(nSubjects,4,4*(ii-1)+2,'NextPlot','add')
			plot(squeeze( mean(stnMeans{ii,1},2)),f)
%			plot(stnMeans{ii,1},f)
			plot([-3 3;-3 3],[0 0;90 90],'k--');	
			plot(significanceMask(:,1),f,'ro','MarkerSize',2);
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
%			plot(stnMeans{ii,2},f)
			plot([-3 3;-3 3],[0 0;90 90],'k--');	
			plot(significanceMask(:,2),f,'ro','MarkerSize',2),2;
			ylim([6 80]);
			xlim([-5 5]);

	end
	annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
	annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');

%			fname = fullfile('/home','lgabri','Desktop','walking',folder,'wavelet',...
%					'avgZScoreStancevSwing.ps');
%
%			print(f2,'-dpsc2',fname);
		
end % function


function pvalue = runPermutationTest(dataObs,stance,swing,nPermutation)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%
		pvalue  = zeros(84,1);
		nStance = size(stance,1);
		nSwing  = size(swing,1);
	
		dataObs = squeeze(mean(dataObs,2));

		% we perform a permutation test for each STN separatelly
		for permIdx = 1:nPermutation

			% concatenate each stride phase
			evGroup = [stance;swing];

			% riffle indices
			riffledIdx = randperm(nStance+nSwing);

			% separate in two groups with same probablilty as
			% observed data (ie. sample size)
			dataStance = evGroup(riffledIdx(1:nStance),:,:);
			dataSwing = evGroup(riffledIdx(nStance+1:end),:,:);

			% compute permutated statistics
			dataPerm = bsxfun(@rdivide,...
					bsxfun(@minus,dataStance,mean(dataSwing,2)),...
					std(dataSwing,[],2));

			dataPerm = squeeze(mean(mean(dataPerm),2));

			pvalue = pvalue + double(dataPerm > dataObs)./nPermutation;

		end

		pvalue = fdrCorrection(pvalue,0.05);

end

function pvalue = fdrCorrection(pvalue, alpha)
%FDRCORRECTION Description
%	PVALUE  = FDRCORRECTION() Long description
%

	tmpPvalue 	= sort(pvalue(:));
	N 					= numel(pvalue);
	FDR 				= alpha.*(1:N)./N;
	thr 				= FDR(find(tmpPvalue <= FDR',1,'last'));
	if isempty(thr)
			pvalue = ones(84,1);
	else
		pvalue(pvalue >= thr) = 1;
end
%	K 					= N*alpha;
%	ttmpPvalue	= sort(tmpPvalue(tmpPvalue < alpha));
%
%	if( numel(ttmpPvalue) > K)
%			b	= ttmpPvalue(end-K:end);
%
%			loc = arrayfun(@(x) find( tmpPvalue == x ), unique(b),'uni',false);
%			loc = cat(1,loc{:});
%			tmpPvalue(loc) = 1;
%			pvalue = reshape(tmpPvalue,size(pvalue));
%	else
%			pvalue = ones(size(pvalue));
%	end
%	
end
