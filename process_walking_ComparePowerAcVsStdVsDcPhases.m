function varargout = process_walking_ComparePowerAcVsStdVsDcPhases( varargin )
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
    sProcess.Comment     = 'Compare Power AcVsStdVsDc';
    sProcess.FileTag     = '__';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Walking';
    sProcess.Index       = 801;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
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
	crossSpectrumOut = cell(nSubjects,1);

	f = 1:60;

	for fileIdx = 1:nFiles; 

			walkingStruct = in_bst_data(sInputs(fileIdx).FileName); 


			channelData		= in_bst_channel(sInputs(fileIdx).ChannelFile);
			iChannels			= channel_find(channelData.Channel,'SEEG');
			signals				= walkingStruct.F(iChannels,:);
			mostAffSide		= mostAffSides(strcmpi(subjectNames,...
																			(sInputs(fileIdx).SubjectName) ));

			if isempty(currentSubject) || ...
							~strcmp(currentSubject,sInputs(fileIdx).SubjectName)

				currentSubject = sInputs(fileIdx).SubjectName;
				subjectIdx = subjectIdx + 1;
				subjectNameOrdered{subjectIdx} = sInputs(fileIdx).SubjectName;
				fprintf('Analyzing %s\n',sInputs(fileIdx).SubjectName);
			end

		  fprintf(' Trial %s\n',sInputs(fileIdx).Condition);
			% compute sampling frequency
			Fs = round(1/mean(diff( walkingStruct.Time )));

			% filter cardiac and peakVelocity events from gait-related events
			evGroupNames = {walkingStruct.Events.label};
			gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel)'));

			% concat all heel contact events in order to have
			% a vector of latencies of this form: e.g.
			% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
			[strideStart,ord]		= sort([walkingStruct.Events(gaitEventGroups).samples]);

			evLabels = cell(1,sum(gaitEventGroups));
			evIdx = find(gaitEventGroups);

			% extract event names
			for iidx = 1:numel(evIdx)
					
					evLabels{iidx} = repmat({walkingStruct.Events(evIdx(iidx)).label},...
							[1 numel(walkingStruct.Events(evIdx(iidx)).samples)]);		

			end

			% re-order event names accordingly
			evNames = [evLabels{:}];
			evNames = evNames(ord);
	
			% count how many strides we have recorded
			nSteps 	= numel(evNames)-1;

			% we have to correct the event adding the offset
			% since they are referred to the 0 of the raw data 
			evOffset 			= round(walkingStruct.Time(1)*Fs);
			strideStart 	= strideStart - evOffset;

			[ftData, DataMat, ChannelMat] = out_fieldtrip_data( sInputs(fileIdx).FileName );
			chancomb = [{channelData.Channel(iChannels).Name}];

%			% quantify steady state duration
%			steadyDuration = (strideStart(end)-2*Fs)-(strideStart(1)+2*Fs);
%
			walkingDuration = strideStart(end)-strideStart(1) ;
%
%			trialString = regexp(sInputs(fileIdx).FileName,'trial\d+','match');
%			fprintf('%s %s %d %f %f\n',sInputs(fileIdx).SubjectName, trialString{:},nSteps,steadyDuration,walkingDuration);

%			if ~(steadyDuration >= 2*Fs && strideStart(end)-strideStart(1) >= 6*Fs)
%					warning('Steady Duration < 2 seconds');
%					continue
%			end
			
			% 3 x 3 rows represent event classes and columns start, end, and offset 
			% for each class
			eventWindow = [strideStart(1), walkingDuration - 1, 0;...
										strideStart(end)-walkingDuration, strideStart(end),0];

			cfg					= [];
			cfg.trl			= eventWindow;
			ftData 			= ft_redefinetrial(cfg,ftData);
			tapNW				= 2;
			
			cfg 				= [];
			cfg.output 	='powandcsd';
			cfg.taper 	= 'dpss';
			cfg.channel = {channelData.Channel(iChannels).Name};
			cfg.channelcmb = chancomb;
			cfg.method  = 'mtmfft';
			cfg.foi 		= f;
      cfg.keeptrials = 'yes';
			cfg.pad 		= 'nextpow2';
			cfg.tapsmofrq = tapNW*Fs/length(ftData.time{1});
			% CrossSpectrum.powspctr tr x 2 x f  
			% 						 .crssspctr tr x 1 x f 
			% complex values
			[CrossSpectrum] = ft_freqanalysis(cfg, ftData);
		
				
			
			crossSpectrumOut(subjectIdx) = {cat(4,crossSpectrumOut{subjectIdx},...
																			cat(2,CrossSpectrum.powspctrm,...
																				abs(CrossSpectrum.crsspctrm)))};

			a = isnan(cat(2,CrossSpectrum.powspctrm,abs(CrossSpectrum.crsspctrm)));
			if any(a(:))
					subjectIdx
			end


	end % for sInputs files

	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
					 'paperunits','normalized','paperorientation',...
						'portrait','Visible','on');

	patientsOrder = {'wue03','wue09','wue04','wue02','wue07','wue06','wue11'};
	[~,ord] = ismember(patientsOrder,subjectNameOrdered);

	normBand = f >= 7 & f<= 60;
	plotIdx = 1;
	for subjectIdx = ord

			data = cat(4,crossSpectrumOut{subjectIdx,:});
			data = mean(data,4);
			normFactor = mean( data(:,:,normBand),3);

			data = data ./ repmat(normFactor,[1 1 numel(f)]);

			subplot(2, nSubjects, plotIdx);
			plot(f, squeeze(data(1,:,:)),'LineWidth',2);
			title(subjectNameOrdered(subjectIdx));
			xlim([6,60])
			ylim([0, 12])
			subplot(2,nSubjects, nSubjects+plotIdx)
			plot(f, squeeze(data(2,:,:)),'LineWidth',2);
			xlim([6,60])
			ylim([0, 12])

			plotIdx = plotIdx + 1;

	end

end % function

function [pvalue, unCorrpvalue] = runPermutationTest(obsERSD,stnData,nPermutation,referenceStance)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%
		pvalue  	= zeros(800,84);
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
			permERSD = computeERSD(dataPerm,referenceStance,2);
				
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
	pvalue(pvalue >= thr) = 1;

end
function stnResult = computeERSD(stnData,referenceStance,method)
%	 COMPUTEERSD of a single STN for a single subject
%

% stnData contains each trial morphed in the 
% => stnData [ n x time x freq ]

% compute the normlization factor concatenating all baseline 
% and computing the mean across trials
	[n,~,f] = size(stnData);
	tBaseline = referenceStance(2):referenceStance(3);
	t = numel(tBaseline);
	
	numFactor = mean(mean(stnData(:,tBaseline,:),1));
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
