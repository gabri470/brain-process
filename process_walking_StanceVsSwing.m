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

	INFO_DIR = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','info');
  trialRejectionFile = fullfile(INFO_DIR,'trialRejection.csv');
  [patNames,trialStrings,stepIds] = textread(trialRejectionFile,'%s %*s %*s %*s%s%d%*s','delimiter',',');

	sideFile = fullfile(INFO_DIR,'patientSides.csv');
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
				fprintf('Analyzing %s\n',sInputs(fileIdx).SubjectName);
			end


			fs = round(1/mean(diff( parentData.Time )));

			% filter cardiac and peakVelocity events from gait-related events
			evGroupNames = {parentData.Events.label};
			gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel|toe)'));

			% concat all heel contact events in order to have
			% a vector of latencies of this form: e.g.
			% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
			[strideStart,ord]		= sort([parentData.Events(gaitEventGroups).samples]);

			evLabels = cell(1,sum(gaitEventGroups));
			evIdx = find(gaitEventGroups);

			% extract event names
			for iidx = 1:numel(evIdx)
			
					
					evLabels{iidx} = repmat({parentData.Events(evIdx(iidx)).label},...
							[1 numel(parentData.Events(evIdx(iidx)).samples)]);		

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

			strideCheck = evNames(bsxfun(@plus,(-2:2),(3:2:numel(evNames)-2)'));
			nStrides = size(strideCheck,1);
			matchingString = {'heelcontact_[L|R]', 'toeoff_[R|L]',...
						'heelcontact_[R|L]','toeoff_[L|R]','heelcontact_[L|R]'};
			orderCheck = nan(1,nStrides);

			for el = 1:nStrides
					orderCheck(el) = sum(cellfun(@isempty,cellfun(@regexp,strideCheck(el,:),...
															matchingString,'uni',false)));
			end

			for strideIdx = 3:2:numel(evNames)-2
					
					% if the stride contains bad steps 
					% we skip it and continue to the next
					if ismember(plotIdx,strideRej) || orderCheck(plotIdx) > 0
							plotIdx = plotIdx + 1;
							continue;
					end


					try
							
							footLabel 		 	= regexp(evNames(strideIdx),'[L|R]','match');
							footLabel				= footLabel{:}{:}; % this is the label of the central hc 
					catch 
							evNames
							strideIdx
					end


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
								stnIdx					= 1;
							end
					end

					%								[ idx ] 	
					% (hc_L) *tof_R  hc_R   *tof_L (hc_L) 
					%    ^ start-2		|t0				      ^ start + 2 == end stride
					timeWindow = strideStart(strideIdx) + (-399:400);
					freqMask 	 = walkingStruct.Freqs > 6;
					try
						dataTF 		 = walkingStruct.TF(:,timeWindow,freqMask);
					catch
							p=0;
					end

					strideRaw	 = signals(:,timeWindow)'.*1e6;
					f 				 = walkingStruct.Freqs(freqMask);	

					% then create the time-warping vector
					try
						originalVector = [1 (strideStart((-2:1:2) + strideIdx) ...
																	- timeWindow(1))  800];
					catch
							p = 0;
					end

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
					end
					% here depending on the options chosen we might 
					% end up having a time-warped data with fixed stride lengths
					% or individual stride lenghts with their original data
					% thus we need to account for different sizes
					%
					
					stnRawData{subjectIdx,stnIdx} = cat(1,stnRawData{subjectIdx,stnIdx},finalTF);

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

	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
			'paperunits','normalized','paperorientation','portrait','Visible','on');
	
	pvalue = zeros(numel(f),2);

	patientsOrder = {'wue03','wue09','wue04','wue02','wue10','wue07','wue06','wue11'};
	[~,ord] = ismember(patientsOrder,subjectNameOrdered);

	plotIdx = 1;
	for ii = ord

			fprintf('Statistcs on %s\n',subjectNameOrdered{ii});

			[corrPvalue(:,1), pvalue(:,1)] = runPermutationTest(stnMeans{ii,1},...
					stanceTcourse{ii,1},swingTcourse{ii,1},100);
			[corrPvalue(:,2), pvalue(:,2)] = runPermutationTest(stnMeans{ii,2},...
					stanceTcourse{ii,2},swingTcourse{ii,2},100);

			
			% need multiple comparison correction  FDR
			corrSignificanceMask = nan(size(corrPvalue));
			corrSignificanceMask(corrPvalue < 0.05) = 4;

			unCorrSignificanceMask = nan(size(pvalue));
			unCorrSignificanceMask(pvalue < 0.05) = 4;


			% this is the STN- 
			annotation('textbox',[0.05, 0.85-(plotIdx-1)*0.1, 0.1, 0.05],...
								'String',subjectNameOrdered{ii},'LineStyle','None');

			subplot(nSubjects,2,2*(plotIdx-1)+1,'NextPlot','add')
			plot(f,squeeze( mean(stnMeans{ii,1},2)))
			plot([0 0;90 90],[-3 3;-3 3],'k--');	
			plot(f,unCorrSignificanceMask(:,1),'k.','MarkerSize',16);
			plot(f,corrSignificanceMask(:,1),'r.','MarkerSize',16);
			xlim([6 80]);
			ylim([-5 5]);


			
			subplot(nSubjects,2,2*(plotIdx-1)+2,'NextPlot','add')
			plot(f,squeeze(mean(stnMeans{ii,2},2)))
			plot([0 0;90 90],[-3 3;-3 3],'k--');	
			plot(f,unCorrSignificanceMask(:,2),'k.','MarkerSize',16);
			plot(f,corrSignificanceMask(:,2),'r.','MarkerSize',16);
			xlim([6 80]);
			ylim([-5 5]);
			plotIdx = plotIdx + 1;

	end
	annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
	annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');

	fname = fullfile('/','home','lgabri','Dropbox','Isaias_group','walking','figs',...
			'avgZScoreStancevSwing_Corr.ps');

	print(f2,'-dpsc2',fname);
		
end % function


function [corrPvalue,unCorrPvalue] = runPermutationTest(dataObs,stance,swing,nPermutation)
%RUNPERMUTATIONTEST 
%	[CORRECTEDPVALUES,UNCORRECTEDPVALUES] = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) 
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
		% output also uncorrected pvalues
		unCorrPvalue = pvalue;
		corrPvalue = fdrCorrection(pvalue,0.05);

end

function pvalue = fdrCorrection(pvalue, alpha)
%FDRCORRECTION 
%	PVALUE  = FDRCORRECTION(pvalue, alpha) 
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

end
