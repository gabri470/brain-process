function varargout = process_walking_SingleVsDouble( varargin )
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
varargout = {};
macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Single Vs Double Supp';
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
		sProcess.options.normOnStride.Comment = 'Normalize on Stride';
		sProcess.options.normOnStride.Type = 'checkbox';
		sProcess.options.normOnStride.Value= false';


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
	nSubjects = numel(unique({sInputs.SubjectName}));
	stnResults = cell(nSubjects,2); 
	stnRawResults = cell(nSubjects,2);
	% the above var holds for each subjects the STN-/+ power for each stride
	
	currentSubject=[];
	subjectIdx = 0;
	subjectNameOrdered = cell(nSubjects,1);
	OutputFiles = {};


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


			% we need to read the angular velocity
			% and in order to do that we need to understand what trial and 
			% which sujbect we are analysing
%			conditionString = regexp(sInputs(fileIdx).Condition,'trial\d+','match');
%		  angVelocityFile = fullfile(DATA_FOLDER,currentSubject,...
%																strcat(currentSubject,'_walking',...
%																				conditionString,'.txt'));


		  fprintf(' Trial %s\n',sInputs(fileIdx).Condition);
			% compute sampling frequency
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
			evNames = [evLabels{:}];
			evNames = evNames(ord);
	
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

			referenceStance 		= 400 + [ -doubleSuppDur-singleSuppDur ...
																				-singleSuppDur 0 +doubleSuppDur ...
																		+doubleSuppDur+singleSuppDur ];

			referenceVector			= [1 referenceStance 800];
			plotIdx = 1;

			strideCheck = evNames(bsxfun(@plus,(-2:2),(3:2:numel(evNames)-2)'));
			nStrides = size(strideCheck,1);
			matchingString = {'heelcontact_[L|R]', 'toeoff_[R|L]',...
						'heelcontact_[R|L]','toeoff_[L|R]','heelcontact_[L|R]'};

			orderCheck = ones(1,nStrides);

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

					%								[ idx ] 	
					% (hc_L) *tof_R  hc_R   *tof_L (hc_L) 
					%    ^ start-2		|t0				      ^ start + 2 == end stride
					timeWindow = strideStart(strideIdx) + (-399:400);
					freqMask 	 = walkingStruct.Freqs > 6;
					dataTF 		 = walkingStruct.TF(:,timeWindow,freqMask);

					if (sProcess.options.normOnStride.Value) 
							% normalize on stride
							normFactor = repmat(mean(dataTF,2),[1 numel(timeWindow) 1]);
							dataTF 		 = (dataTF-normFactor)./normFactor;
					else
							% normalize on swing only
							% this means that we are not normalizing here 
							% but do it later after warping
							% for simplicity in the code. NOT SURE IF THIS IS THE CORRECT WAY
							% TODO check whether we have an effect of normalization order
 
					end


					fprintf('[%d]',strideIdx);
					fprintf('%s ',evNames{(-2:2) + strideIdx});
					fprintf('\n');
					strideRaw	 = signals(:,timeWindow)';
					f 				 = walkingStruct.Freqs(freqMask);	

					% then create the time-warping vector
					originalVector = [1 (strideStart((-2:1:2) + strideIdx)...
															- timeWindow(1))  800];

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

							zLimit 	= [-1 1];
							tAxis		= referenceTimeVector;
							tEvAxis = repmat(referenceTimeVector(referenceStance([2 3 4])),2, 1);

					else
							finalTF 	= dataTF;
							finalRaw	= strideRaw';
							zLimit 		= [min(finalTF(:)) max(finalTF(:))];
							tAxis			= (-399:400)/fs;
							tEvAxis 	= repmat(originalVector([3 4 5])./fs,2, 1);
					end

					footLabel 		 	= regexp(evNames(strideIdx),'[L|R]','match');
					footLabel				= footLabel{:}{:}; 
					% this is the label of the central event 


					if strcmp(footLabel,'L')
							% left foot swing => central hc_L
							%  => stnContra is rightSTN == idx 1
							controLatIdx = 1;
							if strcmp(mostAffSide,'L')
								% if STN- is L => this swing is relative to 
								% STN+ => plot on right side of page
								stnIdx					= 2;
							else 
								% if STN- is R => this swing is relative
								% to STN- => plot on left side of page
								stnIdx					= 1;
							end
					else 
							% right foot swing
							controLatIdx = 2;
							if strcmp(mostAffSide,'L')
								% if STN- is R => this swing is relative
								% to STN- => plot on left side of page
								stnIdx					= 1;
							else 
								% if STN- is L => this swing is relative to 
								% STN+ => plot on right side of page
								stnIdx					= 2;
							end
					end

					
					% we save for each subject the time-warped ERSD normalized over***
					stnRawResults{subjectIdx,stnIdx} = ...
							cat(1,stnRawResults{subjectIdx,stnIdx},...
																								finalTF(controLatIdx,:,:));

					if(~sProcess.options.normOnStride.Value)
							finalTF = bsxfun(@rdivide,bsxfun(@minus,finalTF,...
									mean(finalTF(:,referenceStance(2):referenceStance(3),:),2)),...
									std(finalTF(:,referenceStance(2):referenceStance(3),:),[],2));
					end

					stnResults{subjectIdx,stnIdx} = cat(1,stnResults{subjectIdx,stnIdx},...
																								finalTF(controLatIdx,:,:));


					finalTF = finalTF(:,referenceStance(1):referenceStance(end),:);
%					time 		= referenceTimeVector(referenceStance(1):referenceStance(end));

					plotIdx = plotIdx + 1;
					clear finalTF;

			end % for stride

			clear finalTF;

	end % for sInputs files

	% stnResults holds ERSD for each stride divide as STN-/+
	stnMeans = cellfun(@mean,stnResults,'UniformOutput',false);
	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
							 'paperunits','normalized','paperorientation',...
								'portrait','Visible','on');

	highBetaMask = f >= 6 & f <= 19;
	lowBetaMask = f >= 20 & f <= 35;
	gammaMask  = f > 13 & f < 80;
	for ii = 1:nSubjects

			stnMostAff = squeeze(stnMeans{ii,1});
			stnLessAff = squeeze(stnMeans{ii,2});

			subplot(nSubjects*2,2,4*(ii-1)+1,'NextPlot','add')
			imagesc(tAxis,f,squeeze(stnMeans{ii,1})',zLimit);
			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
			axis xy;
			set(gca,'XTickLabel',[]);
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			ylim([6 80]);

			subplot(nSubjects*2,2,4*(ii-1)+2,'NextPlot','add')
		  imagesc(tAxis,f,squeeze(stnMeans{ii,2})',zLimit);
			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			axis xy;
			set(gca,'XTickLabel',[]);
			ylim([6 80]);


			subplot(nSubjects*2,2,4*(ii-1)+3,'NextPlot','add')
			plot(tAxis,mean(stnMostAff(:,lowBetaMask),2),'r');
			plot(tAxis,mean(stnMostAff(:,highBetaMask),2),'g');
			plot(tAxis,mean(stnMostAff(:,gammaMask),2),'b');
			plot(tEvAxis,repmat([-3;3],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			ylim([-3 3]);

			subplot(nSubjects*2,2,4*(ii-1)+4,'NextPlot','add')
			plot(tAxis,mean(stnLessAff(:,lowBetaMask),2),'r');
			plot(tAxis,mean(stnLessAff(:,highBetaMask),2),'g');
			plot(tAxis,mean(stnLessAff(:,gammaMask),2),'b');
			plot(tEvAxis,repmat([-3;3],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			ylim([-3 3]);



	end

	fname = fullfile('/','home','lgabri','Dropbox','Isaias_group','walking','figs',...
			'avgZScoreSinVsDouble.ps');

	print(f2,'-dpsc',fname);

	for ii = 1:nSubjects

			[pvalue(1,:,:), unCorrPvalue(1,:,:)] = runPermutationTest(stnMeans{ii,1},...
																					stnRawResults{ii,1},100,referenceStance);
			[pvalue(2,:,:), unCorrPvalue(2,:,:)] = runPermutationTest(stnMeans{ii,2},...
																					stnRawResults{ii,2},100,referenceStance);

			statSignificance = ones(size(pvalue)).*0.3;
			statSignificance(unCorrPvalue < 0.05) = 0.6;
			statSignificance(pvalue < 0.05) = 1;

			f2 = figure(2);
			hold on
			h = imagesc(tAxis,f,squeeze(stnMeans{ii,1})',zLimit);
			box off
			axis off

			set(h,'AlphaData',squeeze(statSignificance(1,:,:))')
			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			axis xy;
			set(gca,'XTickLabel',[]);
			ylim([6 80]);

			fname = fullfile('/','home','lgabri','Dropbox','Isaias_group','walking','figs',...
					strcat(subjectNameOrdered{ii},'_stn-_sinVsDouble.png'));

			print(f2,'-dpng','-r300',fname);

			f3 = figure(3);
			hold on
			h = imagesc(tAxis,f,squeeze(stnMeans{ii,2})',zLimit);

			box off
			axis off
			set(h,'AlphaData',squeeze(statSignificance(2,:,:))')
			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			axis xy;
			set(gca,'XTickLabel',[]);
			ylim([6 80]);

			fname = fullfile('/','home','lgabri','Dropbox','Isaias_group','walking','figs',...
					strcat(subjectNameOrdered{ii},'_stn+_sinVsDouble.png'));

			print(f3,'-dpng','-r300',fname);

			clearvars statSignificance pvalue unCorrPvalue
			close(f2);
			close(f3);


	end
end % function

function [pvalue, unCorrpvalue] = runPermutationTest(dataObs,dataRaw,nPermutation,referenceStance)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%
		pvalue  	= zeros(1,800,84);
		nSwing  	= size(dataRaw,1);

		dataPerm	= dataRaw;

		% we perform a permutation test for each STN separatelly
		for permIdx = 1:nPermutation
			
			% for each swing we randomly split the signal in two chunks 
			% and rotate them
			for swingIdx = 1:nSwing

				dataPerm(swingIdx,:,:) = randCircShift(dataRaw(swingIdx,:,:)); 

			end
			
			% compute permutated statistics
			dataPerm = bsxfun(@rdivide,bsxfun(@minus,dataPerm,...
									mean(dataPerm(:,referenceStance(2):referenceStance(3),:),2)),...
									 std(dataPerm(:,referenceStance(2):referenceStance(3),:),[],2));

			% average across swing
			dataPerm = mean(dataPerm);

			% compute pvalues for all frequencies and all time points.
			pvalue = pvalue + double(dataPerm > dataObs)./nPermutation;

		end
		unCorrpvalue = pvalue;
		pvalue = fdrCorrection(pvalue,0.05);

end

function A = randCircShift(A)

		idx 			= randi(size(A,2),1);
		A(1,:,:) 	= cat(2,A(1,idx:end,:),A(1,1:idx-1,:));

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
