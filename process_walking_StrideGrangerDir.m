function varargout = process_walking_StrideGrangerDir( varargin )
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
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Stride Granger';
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
%		sProcess.options.sensortypes.Comment = 'Sensor types or names: ';
%		sProcess.options.sensortypes.Type    = 'text';
%		sProcess.options.sensortypes.Value   = 'SEEG';
%		sProcess.options.doWarping.Comment   = 'time warping ';
%		sProcess.options.doWarping.Type      = 'checkbox';
%		sProcess.options.doWarping.Value     = true;
%		sProcess.options.saveOutput.Comment  = 'Save output to brainstormDB';
%		sProcess.options.saveOutput.Type     = 'checkbox';
%		sProcess.options.saveOutput.Value    = false;
%		sProcess.options.normOnStride.Comment= 'Normalize on Stride';
%		sProcess.options.normOnStride.Type   = 'checkbox';
%		sProcess.options.normOnStride.Value  = false';
%	

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
%	stnResults = cell(nSubjects,2); 
	stnData = cell(nSubjects,2); 
	stnRawResults = cell(nSubjects,2);
	stridePhaseDur = cell(nSubjects,2);
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

			% compute sampling frequency
			fs = round(1/mean(diff( parentData.Time )));

			% filter cardiac events from gait-related events
			evGroupNames = {parentData.Events.label};
			gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,...
					'(heel|toe|peakVeloc)'));

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

			end

			% re-order event names accordingly
			evNames 	= [evLabels{:}];
			evNames 	= evNames(ord);

			% we filter those events that 
			peakVelocIdx 	= find(~cellfun(@isempty,regexp(evNames,'peak')));
			peakVelocMask = peakVelocIdx - 2 > 0 & peakVelocIdx + 5 <= numel(evNames);

			% count how many strides we have recorded
			nStrideLeft = sum(strcmp(evNames,'peakVeloc_L'))-1; 
			nStrideRight= sum(strcmp(evNames,'peakVeloc_R'))-1; 
						
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
			referenceTimeVector = -1:1/fs:(1.5-1/fs);
			doubleSuppDur 			= floor(0.19*fs);
			singleSuppDur 			= floor(0.4*fs);
			acPhaseDur					= floor((.4*2/3)*fs);
			decPhaseDur					= floor((.4/3)*fs);


			referenceStance 		= 400 + [ -doubleSuppDur-acPhaseDur...
																				-acPhaseDur 0 +decPhaseDur ...
																		+doubleSuppDur+decPhaseDur ...
																		+doubleSuppDur+decPhaseDur+acPhaseDur ...
																		+doubleSuppDur+decPhaseDur+singleSuppDur ...
																		+doubleSuppDur+decPhaseDur+singleSuppDur+doubleSuppDur];

%			referenceStance 		= 400 + [ -doubleSuppDur-singleSuppDur ...
%																				-singleSuppDur 0 +doubleSuppDur ...
%																		+doubleSuppDur+singleSuppDur ];

			referenceVector			= [1 referenceStance 1000];
			plotIdx = 1;

%			strideCheck = evNames(bsxfun(@plus,(-2:4),(3:2:numel(evNames)-2)'));
 			strideCheck = evNames(bsxfun(@plus,(-2:5),(peakVelocIdx(peakVelocMask))'));
			% this is the label of the central event 
			nStrides = size(strideCheck,1);

			% check that data are order correctly for each stride
			matchingString = {'heelcontact_[L|R]', 'toeoff_[R|L]',...
																	'peakVeloc_[R|L]',...
													'heelcontact_[R|L]','toeoff_[L|R]','peakVeloc_[L|R]',...
													'heelcontact_[L|R]','toeoff_[R|L]'};
			orderCheck = nan(1,nStrides);

%			for el = 1:nStrides
%					orderCheck(el) = sum(cellfun(@isempty,cellfun(...
%														@regexp,strideCheck(el,:),...
%															matchingString,'uni',false)));
%			end

			for strideIdx = peakVelocIdx(peakVelocMask) 

	
					orderCheck = sum(cellfun(@isempty,cellfun(...
													@regexp,evNames((-2:5)+strideIdx),...
																	matchingString,'uni',false)));

					sizeCheck = strideStart(strideIdx) -499 > 0 & ...
											strideStart(strideIdx) + 500 <= size(walkingStruct.TF,2);

					% if the stride contains bad steps 
					% we skip it and continue to the next
					if ismember(plotIdx,strideRej) || orderCheck > 0 || ~sizeCheck
							plotIdx = plotIdx + 1;
							continue;
					end

					%								[ idx ] 	
					% (hc_L) *tof_R  hc_R   *tof_L (hc_L) 
					%    ^ start-2		|t0				      ^ start + 2 == end stride
					timeWindow = strideStart(strideIdx) + (-499:500);
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


%					fprintf('[%d]',strideIdx);
%					fprintf('%s ',evNames{(-2:4) + strideIdx});
%					fprintf('\n');
					strideRaw	 = signals(:,timeWindow)';
					f 				 = walkingStruct.Freqs(freqMask);	

					% then create the time-warping vector
					originalVector = [1 (strideStart((-2:5) + strideIdx)...
															- timeWindow(1))  1000];


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

%							zLimit 	= [-1 1];
							tAxis		= referenceTimeVector;
							tEvAxis = repmat(referenceTimeVector(referenceStance(2:7)),2, 1);

					else
%							finalTF 	= dataTF.*1e12;
%							finalRaw	= strideRaw';
%							zLimit 		= [min(finalTF(:)) max(finalTF(:))];
%							tAxis			= (-399:400)/fs;
%							tEvAxis 	= repmat(originalVector([3 4 5])./fs,2, 1);
					end

					footLabel 		 	= regexp(evNames(strideIdx),'[L|R]','match');
					footLabel				= footLabel{:}{:}; 
					% this is the label of the central event 

					fprintf('[%d]',strideIdx);
					fprintf('%s \n',evNames{strideIdx});


					if strcmp(footLabel,'L')
							% left foot swing => central hc_L
							%  => stnContra is rightSTN == idx 1
%							controLatIdx = 1;
							stnOrder = [1 2];
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
%							controLatIdx = 2;
							stnOrder = [2 1];
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

					stridePhaseDur{subjectIdx,stnIdx} = ...
							cat(1,stridePhaseDur{subjectIdx,stnIdx},...
																diff(originalVector));
					% we save for each subject the time-frequency data
					% in the way that stnIdx = 1 contains STN- controlateral strides
					% and stnIdx = 2 the strides controlateral to STN+
					% the resuling matrices for a given cell will be 
					% of 2 x T x f x Trials
					stnRawResults{subjectIdx,stnIdx} = ...
							cat(1,stnRawResults{subjectIdx,stnIdx},...
																finalTF(stnOrder,:,:));

					plotIdx = plotIdx + 1;
					clear finalTF;

			end % for stride

			clear finalTF;

	end % for sInputs files

	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
							 'paperunits','normalized','paperorientation',...
								'portrait','Visible','on');
			colormap([1 1 1; 1 1 1; jet(256)]);
%
	highBetaMask = f >= 6 & f <= 19;
	lowBetaMask = f >= 20 & f <= 35;
	gammaMask  = f > 35 & f < 80;

	patientsOrder = {'wue03','wue09','wue04','wue02','wue10','wue07','wue06','wue11'};
	[~,ord] = ismember(patientsOrder,subjectNameOrdered);

	groupData = ones(8,numel(tAxis),84,2);
	plotIdx = 1;

	for ii = ord

			% for a given subjects we should first compute the Wavelet cross coupling matrix
			dataStnMostContra = stnRawResults{ii,1};

			WS(1,2,:,:) = mean(dataStnMostContra(1,:,:,:)*conj(dataStnMostContra(2,:,:,:)),4);
			WS(1,1,:,:) = mean(dataStnMostContra(1,:,:,:)*conj(dataStnMostContra(1,:,:,:)),4);
			WS(2,2,:,:) = mean(dataStnMostContra(2,:,:,:)*conj(dataStnMostContra(2,:,:,:)),4);

			% these two loops would most probably take for ever to run
			for tIdx = 1:size(dataStnMostContra,2)
					[H,Z,S,psi] = sfactorization_wilson2x2(squeeze(WS(:,:,tIdx,:)),f);

%					I = log(WS


			end

			 
 



			subplot(nSubjects*2,2,4*(plotIdx-1)+1,'NextPlot','add')

			imagesc(tAxis,f,(stnMostAffERSD.*squeeze(statSignificance(1,:,:)))');
%			set(h,'AlphaData',squeeze(statSignificance(1,:,:))')
			plot(tEvAxis,repmat([min(f);max(f)],[1 6]),'k--');
			axis xy;
			axis xy;
			set(gca,'XTickLabel',[]);
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			ylim([6 80]);


			subplot(nSubjects*2,2,4*(plotIdx-1)+2,'NextPlot','add')
		  imagesc(tAxis,f,(stnLessAffERSD.*squeeze(statSignificance(2,:,:)))');
%			set(h,'AlphaData',squeeze(statSignificance(2,:,:))')
			plot(tEvAxis,repmat([min(f);max(f)],[1 6]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			axis xy;
			set(gca,'XTickLabel',[]);
			ylim([6 80]);

			% this is the STN- 
			annotation('textbox',[0.05, 0.85-(plotIdx-1)*0.1, 0.1, 0.05],...
								'String',subjectNameOrdered{ii},'LineStyle','None');

			subplot(nSubjects*2,2,4*(plotIdx-1)+3,'NextPlot','add')
			plot(tAxis,mean(stnMostAffERSD(:,highBetaMask),2),'r');
			plot(tAxis,mean(stnMostAffERSD(:,lowBetaMask),2),'g');
			plot(tAxis,mean(stnMostAffERSD(:,gammaMask),2),'b');
%			plot(tEvAxis,repmat([-3;3],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);

			subplot(nSubjects*2,2,4*(plotIdx-1)+4,'NextPlot','add')
			plot(tAxis,mean(stnLessAffERSD(:,highBetaMask),2),'r');
			plot(tAxis,mean(stnLessAffERSD(:,lowBetaMask),2),'g');
			plot(tAxis,mean(stnLessAffERSD(:,gammaMask),2),'b');
%			plot(tEvAxis,repmat([-3;3],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);

			plotIdx = plotIdx + 1;

			% groupData cointains single-subject mean time profiles
			% for 3 frequency bands of interest
%			groupData(ii,:,1,1) = mean(stnMostAff(:,lowBetaMask),2);
%			groupData(ii,:,1,2) = mean(stnLessAff(:,lowBetaMask),2);
%			groupData(ii,:,2,1) = mean(stnMostAff(:,highBetaMask),2);
%			groupData(ii,:,2,2) = mean(stnLessAff(:,highBetaMask),2);
%			groupData(ii,:,3,1) = mean(stnMostAff(:,gammaMask),2);
%			groupData(ii,:,3,2) = mean(stnLessAff(:,gammaMask),2);
%			groupData(ii,:,:,1) = stnMostAff;
%			groupData(ii,:,:,2) = stnLessAff;
		
%			fname = fullfile(getenv('home'),'dropbox','isaias_group',...
%					'walking','figs',strcat(subjectnameordered{ii},'_avgzscorestridemod.ps'));
%
%
%
%			print(f2, '-dpsc2',fname);


	end
	annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
	annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');

	fname = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','figs',...
			'avgZScoreStrideMod.png');

	print(f2,'-dpng',fname);
%
%	grPval = nan(size(groupData,2),3,2);
%	for fIdx= 1:3
%			for stnIdx = 1:2
%					for tIdx = 1:size(groupData,2)
%
%							[~,p] = ttest(squeeze(groupData(:,tIdx,fIdx,stnIdx)),0);
%							grPval(tIdx,fIdx,stnIdx) = fdrCorrection(p,0.05);
%
%					end
%			end
%	end

%	[h,p] = ttest(groupData);
%
%	grPval(grPval >= 0.05) = nan;
%	grPval(grPval < 0.05) = 1;
%
%	fn = figure;
%	subplot(2,1,1,'NextPlot','Add'),
%%			plot(tAxis,squeeze(mean(groupData(:,:,:,1))));
%%			plot(tAxis,squeeze(grPval(:,:,1)).*repmat([2 2.2 2.4],...
%%					[numel(tAxis),1]),'.','MarkerSize',2);
%
%			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
%	subplot(2,1,2,'NextPlot','Add'),
%%			plot(tAxis,squeeze(mean(groupData(:,:,:,2))));
%%			plot(tAxis,squeeze(grPval(:,:,2)).*repmat([2 2.2 2.4],...
%%					[numel(tAxis),1]),'.','MarkerSize',2);
%
%			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
%
%%
%%	fname = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','figs',...
%%			'grLvlAcVsDCPhases.ps');
%%
%%	print(fn,'-dpsc',fname);
%%
%%
%	for ii = 1:nSubjects
%
%			[pvalue(1,:,:), unCorrPvalue(1,:,:)] = runPermutationTest(stnMeans{ii,1},...
%																					stnRawResults{ii,1},100,referenceStance);
%			[pvalue(2,:,:), unCorrPvalue(2,:,:)] = runPermutationTest(stnMeans{ii,2},...
%																					stnRawResults{ii,2},100,referenceStance);
%
%
%			statSignificance = ones(size(pvalue)).*0.3;
%			statSignificance(unCorrPvalue < 0.05) = 0.6;
%			statSignificance(pvalue < 0.05) = 1;
%
%
%			f2 = figure(2);
%			hold on
%			h = imagesc(tAxis,f,squeeze(stnMeans{ii,1})',zLimit);
%			box off
%			axis off
%
%			set(h,'AlphaData',squeeze(statSignificance(1,:,:))')
%			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
%			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
%			axis xy;
%			set(gca,'XTickLabel',[]);
%			ylim([6 80]);
%
%			fname = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','figs',...
%					strcat(subjectNameOrdered{ii},'_stn-_stat.png'));
%
%%			print(f2,'-dpng','-r300',fname);
%
%			f3 = figure(3);
%			hold on
%			h = imagesc(tAxis,f,squeeze(stnMeans{ii,2})',zLimit);
%
%			box off
%			axis off
%			set(h,'AlphaData',squeeze(statSignificance(2,:,:))')
%			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
%%			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
%			axis xy;
%			set(gca,'XTickLabel',[]);
%			ylim([6 80]);
%
%			fname = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','figs',...
%					strcat(subjectNameOrdered{ii},'_stn+_stat.png'));
%
%%			print(f3,'-dpng','-r300',fname);
%
%			clearvars statSignificance pvalue unCorrPvalue
%			close(f2);
%			close(f3);
%
%
%	end
end % function

function [pvalue, unCorrpvalue] = runPermutationTest(obsERSD,stnData,nPermutation,referenceStance)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%
		pvalue  	= zeros(1000,84);
		nSwing  	= size(stnData,1);

		dataPerm  = stnData;

		% we perform a permutation test for each STN separatelly
		for permIdx = 1:nPermutation
			
			% for each swing we randomly split the signal in two chunks 
			% and rotate them
			for swingIdx = 1:nSwing

				dataPerm(swingIdx,:,:) = randCircShift(stnData(swingIdx,:,:)); 

			end
			
			% compute permutated statistics
			permERSD = computeERSD(dataPerm,referenceStance,1);
				
			% compute pvalues for all frequencies and all time points.
			pvalue = pvalue + double(permERSD > obsERSD)./nPermutation;

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
	if ~isempty(thr)
		pvalue(pvalue >= thr) = 1;
	else
		pvalue = ones(size(pvalue));
	end
		

end

function stnResult = computeERSD(stnData,referenceStance,method)
%	 COMPUTEERSD of a single STN for a single subject
%

% stnData contains each trial morphed in the 
% => stnData [ n x time x freq ]

% compute the normlization factor concatenating all baseline 
% and computing the mean across trials
	[n,~,f] = size(stnData);
	tBaseline = referenceStance(1):referenceStance(3);
	t = numel(tBaseline);
	
	numFactor = mean(mean(stnData(:,tBaseline,:),2));
	denFactor = std(reshape(stnData(:,tBaseline,:),[n*t,f]));

	if method
		% rel change
		stnResult = bsxfun(@rdivide,bsxfun(@minus,stnData,numFactor),numFactor);
	else
		% pseudo-zscore
		stnResult = bsxfun(@rdivide,bsxfun(@minus,stnData,numFactor),denFactor);
	end

	stnResult = squeeze(mean(stnResult));

end
