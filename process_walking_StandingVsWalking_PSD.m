function varargout = process_walking_StandingVsWalking_PSD( varargin )
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
varargout =  {};
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Standing vs Walking';
    sProcess.FileTag     = '__';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Walking';
    sProcess.Index       = 801;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

	

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

	DATA_FOLDER = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','info');

	sideFile = fullfile(DATA_FOLDER,'patientSides.csv');
	[nameSubjects, mostAffSides] = textread(sideFile,'%s %s\n','delimiter',',');

	subjectNames = unique({sInputs.SubjectName});
	nSubjects = numel(subjectNames);

	mostAffSides = mostAffSides(ismember(nameSubjects,subjectNames));

	% we need to find files that have to be processed and group 
	% them for subjects and conditions
	% First group all sInput.Comment together
	conditionStrings = {sInputs.Condition};
	
	standingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(s|S)tanding'));
	walkingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(w|W)alking'));
	restingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(r|R)esting'));

	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
							 'paperunits','normalized','paperorientation',...
								'portrait','Visible','on');
	f = 1:60;

	OutputFiles = {};

	for subjIdx = 1:nSubjects
			% for each subject separately we pick standing condition
			subjectMask 		= ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));

			standingFileIdx = find(subjectMask & standingConMask);
			walkingFileIdx 	= find(subjectMask & walkingConMask);
			restingFileIdx 	= find(subjectMask & restingConMask);

			% read standing time-freq data
			standingStruct 	= in_bst_data(sInputs(standingFileIdx).FileName);
			standChannels 	= in_bst_channel(sInputs(standingFileIdx).ChannelFile);
			standiChannels	= channel_find( standChannels.Channel,'SEEG');

			standSpectrum   = computeSpectrum( sInputs(standingFileIdx).FileName,standChannels,standiChannels);

			% read resting time-freq data
			restingStruct 	= in_bst_data(sInputs(restingFileIdx).FileName);
			restChannels 		= in_bst_channel(sInputs(restingFileIdx).ChannelFile);
			restiChannels		= channel_find( restChannels.Channel,'SEEG');

			restSpectrum 		= computeSpectrum( sInputs(restingFileIdx).FileName, restChannels, restiChannels);

			walkSpectrum 		= nan(numel(f),numel(walkingFileIdx));

			for walkIdx = 1:numel(walkingFileIdx)

					walkingStruct = in_bst_data(sInputs(walkingFileIdx(walkIdx)).FileName);
					walkChannels 	= in_bst_channel(sInputs(walkingFileIdx(walkIdx)).ChannelFile);
					walkiChannels	= channel_find(walkChannels.Channel,'SEEG');
		
					% filter cardiac and peakVelocity events from gait-related events
					evGroupNames = {walkingStruct.Events.label};
					gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel)'));

					Fs = 1/mean(diff(walkingStruct.Time));

					% concat all heel contact events in order to have
					% a vector of latencies of this form: e.g.
					% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
					eventSamples = sort([walkingStruct.Events(gaitEventGroups).samples]);

					% we have to correct the event adding the offset
					% since they are referred to the 0 of the raw data 
					evOffset 			= round(walkingStruct.Time(1)*Fs);
					eventSamples  = eventSamples - evOffset;

					% analysis window is from the first heel contact to the last heel contact
					timeWindow 		= eventSamples(1):eventSamples(end);

					% we then pack trails together in order to have a matrix 2 x walkingRefLength x f x trials
					% that we will rotate to match the form 2 x windows x walkingRefLength x f as 
					% standing condition data are represented in
					walkSpectrum(:,walkIdx) = computeSpectrum(sInputs(walkingFileIdx(walkIdx)).FileName,...
																										walkChannels,walkiChannels,timeWindow);


		  end % walking trial loop


			subplot(2,4,subjIdx,'NextPlot','add')
			plot(f,restSpectrum,'LineWidth',2);
			plot(f,standSpectrum,'LineWidth',2);
			plot(f,mean(walkSpectrum,2),'LineWidth',2);	
			xlim([6 60]);
			ylim([0 12]);
			title(subjectNames(subjIdx));

	end % subject loop
	clearvars walkData standData
end % function


function [pvalue, unCorrpvalue] = runPermutationTest(dataObs,dataA,dataB,nPermutation,alpha)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%
		pvalue  	 = zeros(1,90);
		nStanding  = size(dataB,2);
		pooledData = cat(2,dataA,dataB);

		% we perform a permutation test for each STN separatelly
		for permIdx = 1:nPermutation

			riffledIndices = randperm(size(pooledData,2));

			standPsd = pooledData(:,riffledIndices(1:nStanding),:);
			walkPsd	= pooledData(:,riffledIndices(nStanding+1:end),:);
			
			% compute permutated statistics
			dataPerm = (squeeze(mean(walkPsd,2)) - squeeze(mean(standPsd,2)))./squeeze(mean(standPsd,2));
		
			% compute pvalues for all frequencies and all time points.
			pvalue = pvalue + double(dataPerm' > dataObs)./nPermutation;

		end
		unCorrpvalue = pvalue;
		pvalue = fdrCorrection(pvalue,alpha);

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
function crossSpectrum = computeSpectrum(filename,channelData,iChannels,timeWindow)
% Description
%	CROSSSPECTRUM = () Long description
%
%

	[ftData, DataMat, ChannelMat] = out_fieldtrip_data( filename );

	if nargin < 4
			timeWindow = ftData.time;
	end
	tapNW						= 2;
	f 							= 1:60;
	chancomb 				= [{channelData.Channel(iChannels).Name}];
    
  Fs = 1/mean(diff(ftData.time{1}));

	cfg 						= [];
	cfg.output 			='powandcsd';
	cfg.taper 			= 'dpss';
	cfg.channel 		= {channelData.Channel(iChannels).Name};
	cfg.channelcmb 	= chancomb;
	cfg.method  		= 'mtmfft';
	cfg.foi 				= 1:60;
	cfg.pad 				= 'nextpow2';
	cfg.tapsmofrq 	= tapNW*Fs/length(ftData.time{1});
	% CrossSpectrum.powspctr tr x 2 x f  
	% 						 .crssspctr tr x 1 x f 
	% complex values
	[CrossSpectrum] = ft_freqanalysis(cfg, ftData);

	normBand 				= f >= 7 & f<= 60;
	normFactor			= mean(abs(CrossSpectrum.crsspctrm(1,normBand)),2);
	crossSpectrum   = abs(CrossSpectrum.crsspctrm)./repmat(normFactor,[1 numel(f)]);

end


