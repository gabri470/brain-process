function varargout = process_walking_phaseTE_test( varargin )

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
    sProcess.Comment     = 'phasTE';
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
	conditionStrings 	= {sInputs.Condition};
	
	standingConMask 	= ~cellfun(@isempty,regexp(conditionStrings,'(s|S)tanding'));
	walkingConMask 		= ~cellfun(@isempty,regexp(conditionStrings,'(w|W)alking'));
	restingConMask 		= ~cellfun(@isempty,regexp(conditionStrings,'(r|R)esting'));

	f1 = figure('papertype','a4','paperposition',[0 0 1 1],...
									 'paperunits','normalized','paperorientation',...
										'portrait','visible','on');
	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
									 'paperunits','normalized','paperorientation',...
										'portrait','Visible','on');


	OutputFiles = {};

	for subjIdx = 1:nSubjects
			% for each subject separately we pick standing condition
			subjectMask 		= ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));

			standingFileIdx = find(subjectMask & standingConMask);
			walkingFileIdx 	= find(subjectMask & walkingConMask);
			restingFileIdx 	= find(subjectMask & restingConMask);

			% get bad channels	
			standChFlag			= getfield(in_bst_data(sInputs(standingFileIdx).FileName,'ChannelFlag'),'ChannelFlag');

			[standingStruct,~,standChannels]  = out_fieldtrip_data(sInputs(standingFileIdx).FileName);
			standChannels 	= in_bst_channel(sInputs(standingFileIdx).ChannelFile);
			standiChannels	= channel_find( standChannels.Channel,'SEEG');


			standingStruct  = preproc(standingStruct,standiChannels,standChFlag,3);
			% we should here quantify the STN coherence
			standCoh				= computeCoherence(standingStruct,[{standChannels.Channel(standiChannels).Name}]);

			restingStruct   = struct('dimord',[],'trial',[],'time',[],'label',[],'elec',[]);

			% resting recordings can be more than 1
			for restIdx = 1:numel(restingFileIdx)
					% read resting time-freq data
					restChFlag 			= getfield(...
															in_bst_data(sInputs(restingFileIdx(restIdx)).FileName,...
																						'ChannelFlag'),'ChannelFlag');


					[restIdxStruct,~,restChannels]  = out_fieldtrip_data(sInputs(restingFileIdx(restIdx)).FileName);
					restChannels 		= in_bst_channel(sInputs(restingFileIdx(restIdx)).ChannelFile);
					restiChannels		= channel_find(restChannels.Channel,'SEEG');
					restIdxStruct  	= preproc(restIdxStruct,restiChannels,restChFlag,3);

					if isempty(restIdxStruct)
							continue
					end

					if restIdx == 1
							restingStruct = restIdxStruct;
					else
							restingStruct.trial = [restingStruct.trial restIdxStruct.trial];
							restingStruct.time	= [restingStruct.time restIdxStruct.time];
					end
			end

			% compute the STN coherence during resting
			restingStruct = rmfield(restingStruct,'sampleinfo');
			restingStruct.hdr.nTrials = numel(restingStruct.trial);
			restCoh	= computeCoherence(restingStruct,restingStruct.label);

			walkingStruct   = struct('dimord',[],'trial',[],'time',[],'label',[],'elec',[]);

			for walkIdx = 1:numel(walkingFileIdx)
					walkChFlag 			= getfield(...
															in_bst_data(sInputs(walkingFileIdx(walkIdx)).FileName,...
																						'ChannelFlag'),'ChannelFlag');

					[walkIdxStruct,~,walkChannels]  = out_fieldtrip_data(sInputs(walkingFileIdx(walkIdx)).FileName);
					walkChannels 	= in_bst_channel(sInputs(walkingFileIdx(walkIdx)).ChannelFile);
					walkiChannels	= channel_find(walkChannels.Channel,'SEEG');
				
					walkIdxStruct = preproc(walkIdxStruct,walkiChannels,walkChFlag,3);
					if isempty(walkIdxStruct)
							continue
					end

					if walkIdx == 1
							walkingStruct = walkIdxStruct;
					else
							walkingStruct.trial = [walkingStruct.trial walkIdxStruct.trial];
							walkingStruct.time	= [walkingStruct.time walkIdxStruct.time];
					end
		
		  end % walking trial loop

			walkingStruct = rmfield(walkingStruct,'sampleinfo');
			walkingStruct.hdr.nTrials = numel(walkingStruct.trial);
			walkCoh	= computeCoherence(walkingStruct,walkingStruct.label);

			ax1 = subplot(2,4,subjIdx,'NextPlot','add');
			plot(restCoh.freq,abs(restCoh.cohspctrm),'LineWidth',2);
			plot(standCoh.freq,abs(standCoh.cohspctrm),'LineWidth',2);
			plot(walkCoh.freq,abs(walkCoh.cohspctrm),'LineWidth',2);	
			legend({'rest','stand','walk'});		
			xlim([6 60]);
			title(subjectNames(subjIdx));
			set(ax1,'Parent',f1);


			ax2 = subplot(2,4,subjIdx,'NextPlot','add');
			plot(restCoh.freq,abs(imag(restCoh.cohspctrm)),'LineWidth',2);
			plot(standCoh.freq,abs(imag(standCoh.cohspctrm)),'LineWidth',2);
			plot(walkCoh.freq,abs(imag(walkCoh.cohspctrm)),'LineWidth',2);	
			legend({'rest','stand','walk'});
			xlim([6 60]);
			title(subjectNames(subjIdx));
			set(ax2,'Parent',f2);

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

function coh = computeCoherence(data,channelNames)


	cfg = [];

	% check if any trial contains NaN values and discard it
	trlMask					= cellfun(@(x) sum(isnan(x),2),data.trial,'uni',false);
	trlMask 				= reshape(cat(1,trlMask{:}),2,numel(data.trial));
	trlMask					= sum(trlMask,1) == 0;
	cfg.trials			= trlMask;
	
	fs 							= 1/mean(diff(data.time{1}));
	cfg.channel 		= channelNames;
	cfg.channelcmb 	= channelNames';
	cfg.output 			='powandcsd';
	cfg.taper 			= 'dpss';
	tapNW						= 2;
	cfg.keeptrials 	= 'yes';
	cfg.method  		= 'mtmfft';
	cfg.foi 				= 1:3:60;
	cfg.pad 				= 'nextpow2';
	cfg.tapsmofrq 	= tapNW*fs/length(data.time{1});

	freq						= ft_freqanalysis(cfg, data);

	cfg.method 			= 'coh';
	cfg.complex 		= 'complex';
    
	coh 						= ft_connectivityanalysis(cfg,freq);

end


function pTE=phaseTE(Xf,lag)
% Inputs:
% Xf is 2D matrix, i.e. channels x samples
% lag is the expected delay, expressed in samples
% Outputs:
% pTE(1,2) is phase TE from 1 to 2
% pTE(2,1) is phase TE from 2 to 1
%
% Function written by Nitin Williams, with help from Felix Siebenhuehner and Muriel Lobier 
%
data 						= ft_preprocessing(cfg,data);
		%% Inst. phase
		Xfh=angle(hilbert(Xf'))';

		%% Number of bins

		for idx=1:size(Xfh,1),
				[~,dirsd(idx,1)]=circ_std(Xfh(idx,:),[],[],2);
		end

		N=size(Xf,2);
		bw=(((3.49).*dirsd)./(N^(1/3)));
		numbins=round((2*pi)./bw);

		%% Computing normalised histograms

		tX=Xfh(1,1+lag:end);
		tY=Xfh(2,1+lag:end);
		tXd=Xfh(1,1:end-lag);
		tYd=Xfh(2,1:end-lag);

		tXh=discretize(tX,linspace(-pi,pi,numbins(1,1)));
		tYh=discretize(tY,linspace(-pi,pi,numbins(2,1)));
		tXdh=discretize(tXd,linspace(-pi,pi,numbins(1,1)));
		tYdh=discretize(tYd,linspace(-pi,pi,numbins(2,1)));

		[P_YYd,~,~]=crosstab(tYh,tYdh); P_YYd=P_YYd./sum(P_YYd(:));
		[P_YdXd,~,~]=crosstab(tYdh,tXdh); P_YdXd=P_YdXd./sum(P_YdXd(:));
		[P_Yd,~,~]=crosstab(tYdh); P_Yd=P_Yd./sum(P_Yd(:));
		[P_YYdXd,~,~]=crosstab(tYh,tYdh,tXdh); P_YYdXd=P_YYdXd./sum(P_YYdXd(:));

		[P_XXd,~,~]=crosstab(tXh,tXdh); P_XXd=P_XXd./sum(P_XXd(:));
		[P_XdYd,~,~]=crosstab(t101Xdh,tYdh); P_XdYd=P_XdYd./sum(P_XdYd(:));
		[P_Xd,~,~]=crosstab(tXdh); P_Xd=P_Xd./sum(P_Xd(:));
		[P_XXdYd,~,~]=crosstab(tXh,tXdh,tYdh); P_XXdYd=P_XXdYd./sum(P_XXdYd(:));

		%% Entropy & pTE calculation

		H_YYd = sum(-(P_YYd(P_YYd>0).*(log2(P_YYd(P_YYd>0)))));
		H_YdXd = sum(-(P_YdXd(P_YdXd>0).*(log2(P_YdXd(P_YdXd>0)))));
		H_Yd = sum(-(P_Yd(P_Yd>0).*(log2(P_Yd(P_Yd>0)))));
		H_YYdXd = sum(-(P_YYdXd(P_YYdXd>0).*(log2(P_YYdXd(P_YYdXd>0)))));

		pTE(1,1)=NaN;
		pTE(2,2)=NaN;
		pTE(1,2)=H_YYd+H_YdXd-H_Yd-H_YYdXd;

		H_XXd = sum(-(P_XXd(P_XXd>0).*(log2(P_XXd(P_XXd>0)))));
		H_XdYd = sum(-(P_XdYd(P_XdYd>0).*(log2(P_XdYd(P_XdYd>0)))));
		H_Xd = sum(-(P_Xd(P_Xd>0).*(log2(P_Xd(P_Xd>0)))));
		H_XXdYd = sum(-(P_XXdYd(P_XXdYd>0).*(log2(P_XXdYd(P_XXdYd>0)))));

		pTE(2,1)=H_XXd+H_XdYd-H_Xd-H_XXdYd;

end

function data = preproc(data,iChannels,chFlag,trialLength)
%	preproc  
%	data = preproc() split data in trials
%

	if any(ismember(find(chFlag==-1),iChannels))
			data = [];
			return
	end

	cfg = [];
	% at this point the recordings are just a single 
	% continuos stream of samples. 
	begRecording 		= min(data.time{1});
	endRecording 		= max(data.time{1});

	fs 							= 1/mean(diff(data.time{1}));

	% we should split this in ntrials of 3s
	totLength 			= endRecording-begRecording;
	nTrials   			= floor(totLength/trialLength);
	offset 					= totLength/2;
	analysisWindow	= nTrials*trialLength;
	startTime 			= offset-analysisWindow/2;
	endTime   			= analysisWindow/2+offset;

	trials					= cat(2,(startTime:trialLength:endTime-trialLength)',...
													(startTime+trialLength:trialLength:endTime)',zeros(nTrials,1));

	trials					= floor(trials*fs)+1;
	cfg.trl					= trials;

	data						= ft_redefinetrial(cfg,data);

	cfg 						= [];
	cfg.continuous  = 'yes';
	cfg.channel			=	data.label(iChannels); 
	data 						= ft_preprocessing(cfg,data);

end

