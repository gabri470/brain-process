function varargout = process_walking_StandingVsWalking_PSD( varargin )

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
    sProcess.Comment     = 'Compare PSDs across Conditions';
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
function OutputFiles = Run(~, sInputs) %#ok<DEFNU>

    DATA_FOLDER = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','info');

    sideFile = fullfile(DATA_FOLDER,'patientSides.csv');
    [nameSubjects, mostAffSides] = textread(sideFile,'%s %s\n','delimiter',',');

    subjectNames = unique({sInputs.SubjectName});
    nSubjects = numel(subjectNames);

    mostAffSides = mostAffSides(ismember(nameSubjects,subjectNames));
    
    patientsOrder = {'wue03','wue09','wue04','wue02','wue10','wue07','wue06','wue11'};
	[~,subjOrder] = ismember(patientsOrder,subjectNames);

    % we need to find files that have to be processed and group
    % them for subjects and conditions
    % First group all sInput.Comment together
    conditionStrings = {sInputs.Condition};

    standingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(s|S)tanding'));
    walkingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(w|W)alking'));
    restingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(r|R)esting'));

    figure('papertype','a4','paperorientation','portrait','Visible','on');
    f = 1:60;

    OutputFiles = {};

    standingData = cell(nSubjects,2);
    walkingData = cell(nSubjects,2);
    restingData = cell(nSubjects,2);
    plotIdx = 1;

    for subjIdx = subjOrder
        
        % for each subject separately we pick standing condition
        subjectMask 		= ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));

        standingFileIdx = find(subjectMask & standingConMask);
        walkingFileIdx 	= find(subjectMask & walkingConMask);
        restingFileIdx 	= find(subjectMask & restingConMask);

        % STN- is the left most ( first index ) element of the arrays
        if strcmp(mostAffSides(subjIdx), 'R')
            chOrder = [1 2];
        else
            chOrder = [2 1];
        end

        % read standing time-freq data
        %standingStruct 	= in_bst_data(sInputs(standingFileIdx).FileName);
        standChannels 	= in_bst_channel(sInputs(standingFileIdx).ChannelFile);
        standiChannels	= channel_find( standChannels.Channel,'SEEG');

        standSpectrum   = computeSpectrum( sInputs(standingFileIdx).FileName,standChannels,standiChannels,chOrder);
        standingData(subjIdx,:) = [{standSpectrum(1,:,:)},{standSpectrum(2,:,:)}];

        restSpectrum    = zeros(2,numel(f),numel(restingFileIdx));
        for restIdx = 1:numel(restingFileIdx)
            % read resting time-freq data
            %restingStruct 	= in_bst_data(sInputs(restingFileIdx(restIdx)).FileName);
            restChannels 	= in_bst_channel(sInputs(restingFileIdx(restIdx)).ChannelFile);
            restiChannels	= channel_find( restChannels.Channel,'SEEG');

            fprintf('Analysing: %s with %d NaNs\n', sInputs(restingFileIdx(restIdx)).FileName,sum(sum(isnan(restSpectrum(:,restIdx)))));
            restSpectrum(:,:,restIdx)	= computeSpectrum( sInputs(restingFileIdx(restIdx)).FileName, restChannels, restiChannels,chOrder);

        end
        restingData{subjIdx} = restSpectrum;
        restingData(subjIdx,:) = [{restSpectrum(1,:,:)},{restSpectrum(2,:,:)}];

        walkSpectrum = nan(2,numel(f),numel(walkingFileIdx));

        for walkIdx = 1:numel(walkingFileIdx)

            %walkingStruct   = in_bst_data(sInputs(walkingFileIdx(walkIdx)).FileName);
            walkChannels 	= in_bst_channel(sInputs(walkingFileIdx(walkIdx)).ChannelFile);
            walkiChannels	= channel_find(walkChannels.Channel,'SEEG');


            % filter cardiac and peakVelocity events from gait-related events
            %evGroupNames = {walkingStruct.Events.label};
            %gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel)'));

            %&Fs = 1/mean(diff(walkingStruct.Time));

            % concat all heel contact events in order to have
            % a vector of latencies of this form: e.g.
            % hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
            %eventSamples = sort([walkingStruct.Events(gaitEventGroups).samples]);

            % we have to correct the event adding the offset
            % since they are referred to the 0 of the raw data
            %evOffset 	   = round(walkingStruct.Time(1)*Fs);
            %eventSamples  = eventSamples - evOffset;

            % analysis window is from the first heel contact to the last heel contact
            %timeWindow 		= eventSamples(1):eventSamples(end);

            % we then pack trails together in order to have a matrix 2 x walkingRefLength x f x trials
            % that we will rotate1 to match the form 2 x windows x walkingRefLength x f as
            % standing condition data are represented in
            walkSpectrum(:,:,walkIdx) = computeSpectrum(sInputs(walkingFileIdx(walkIdx)).FileName,...
                walkChannels,walkiChannels,chOrder);
            fprintf('Analysing: %s with %d NaNs\n', sInputs(walkingFileIdx(walkIdx)).FileName,...  
                           sum(sum(isnan(walkSpectrum(:,walkIdx)))));


        end % walking trial loop
        walkingData(subjIdx,:) = [{walkSpectrum(1,:,:)},{walkSpectrum(2,:,:)}];
        
        

  
        subplot(nSubjects,2,2*(plotIdx-1)+1,'NextPlot','add')
        plot(f,standSpectrum(chOrder(1),:,:),'LineWidth',2,'Color',[255 109 182]./255);
        plot(f,squeeze(mean(restSpectrum(chOrder(1),:,:),3)),'LineWidth',2,'Color',[0 109 219]./255);
        plot(f,mean(walkSpectrum(chOrder(1),:,:),3),'LineWidth',2,'Color',[76 255 36]./255);
        xlim([6 60]);
        ylim([0 12]);
        xlabel('Freq. Hz');
        ylabel('norm. pow');
        
        subplot(nSubjects,2,2*(plotIdx-1)+2,'NextPlot','add')
        plot(f,standSpectrum(chOrder(2),:,:),'LineWidth',2,'Color',[255 109 182]./255);
        plot(f,squeeze(mean(restSpectrum(chOrder(2),:,:),3)),'LineWidth',2,'Color',[0 109 219]./255);
        plot(f,mean(walkSpectrum(chOrder(2),:,:),3),'LineWidth',2,'Color',[76 255 36]./255);

        xlim([6 60]);
        ylim([0 12]);
        
        xlabel('Freq. Hz');
        ylabel('norm. pow');
        
        plotIdx = plotIdx + 1;



    end % subject loop
    clearvars walkData standData

    % contains PSD data across subjects for ordered STN along the second
    % dimension ( nSubj x STN ). Each element of cell array contains the
    % PSD 1 x f x trial
    
    % take avg across trial
    standingData = cellfun(@(x) mean(x,3),standingData,'UniformOutput',false);
    restingData  = cellfun(@(x) mean(x,3),restingData,'UniformOutput',false);
    walkingData  = cellfun(@(x) mean(x,3),walkingData,'UniformOutput',false);    

    
    standingData = reshape(cat(1,standingData{:}),nSubjects,2,numel(f));
    restingData = reshape(cat(1,restingData{:}),nSubjects,2,numel(f));
    walkingData = reshape(cat(1,walkingData{:}),nSubjects,2,numel(f));
    
    [standingConfLim,standingMeans] = myBootstrap(standingData,nSubjects,10);
    [restingConfLim,restingMeans] = myBootstrap(restingData,nSubjects,10);
    [walkingConfLim,walkingMeans] = myBootstrap(walkingData,nSubjects,10);
    
    figure
    subplot(2,1,1,'NextPlot','add')
    title('STN-')
    plot(f,squeeze(standingMeans(1,1,:)),'LineWidth',2,'Color',[255 109 182]./255);
    plot(f,squeeze(restingMeans(1,1,:)),'LineWidth',2,'Color',[0 109 219]./255);
    plot(f,squeeze(walkingMeans(1,1,:)),'LineWidth',2,'Color',[76 255 36]./255);
    
    legend({'stand','rest','walk'})
    fill_between(f,squeeze(standingConfLim(1,1,:)),squeeze(standingConfLim(2,1,:)),f,'EdgeColor','none','FaceAlpha',0.2,'FaceColor',[255 109 182]./255);  
    fill_between(f,squeeze(restingConfLim(1,1,:)),squeeze(restingConfLim(2,1,:)),f,'EdgeColor','none','FaceAlpha',0.2,'FaceColor',[0 109 219]./255);
    fill_between(f,squeeze(walkingConfLim(1,1,:)),squeeze(walkingConfLim(2,1,:)),f,'EdgeColor','none','FaceAlpha',0.2,'FaceColor',[76 255 36]./255);
    xlabel('Freq. (Hz)');
    ylabel('Norm. Power');
    ylim([0 12])
    xlim([5 60])
    
    subplot(2,1,2,'NextPlot','add')
    title('STN+')
    plot(f,squeeze(standingMeans(1,2,:)),'LineWidth',2,'Color',[255 109 182]./255);
    plot(f,squeeze(restingMeans(1,2,:)),'LineWidth',2,'Color',[0 109 219]./255);
    plot(f,squeeze(walkingMeans(1,2,:)),'LineWidth',2,'Color',[76 255 36]./255);
    legend({'stand','rest','walk'})
    fill_between(f,squeeze(standingConfLim(1,2,:)),squeeze(standingConfLim(2,2,:)),f,'EdgeColor','none','FaceAlpha',0.2,'FaceColor',[255 109 182]./255);  
    fill_between(f,squeeze(restingConfLim(1,2,:)),squeeze(restingConfLim(2,2,:)),f,'EdgeColor','none','FaceAlpha',0.2,'FaceColor',[0 109 219]./255);
    fill_between(f,squeeze(walkingConfLim(1,2,:)),squeeze(walkingConfLim(2,2,:)),f,'EdgeColor','none','FaceAlpha',0.2,'FaceColor',[76 255 36]./255);
    xlabel('Freq. (Hz)');
    ylabel('Norm. Power');
    
    ylim([0 12])
    xlim([5 60])

    
end % function

function [confLimit,dataMean] = myBootstrap(data,nSubject,nBootstraps)

    % data matrix has nSubjects x nStn (ordered (stn-/stn+) ) x f
    % bootstrap the C.L. for mean
    bootstrapIndexes = randi(nSubject,nBootstraps);
    
    dataMean = mean(data);
    % currBootstraps will contain nBootstraps x nStn x f
    currBootstraps = nan(nBootstraps,2,60);
    for idx = 1:nBootstraps
        currBootstraps(idx,:,:) = mean(data(bootstrapIndexes(idx),:,:),1);
    end
    
    % confLimit will contain [UB LB] x nStn x f
    % dataMean is 1 x 2 x f
    % thus we replicate on the first dim
    A = repmat(dataMean,2,1,1);
    % currBootstrap after mean should be 1 x 2 x f
    % thus we replicate on the first dim
    B = prctile(currBootstraps,[5 95]);
    % coeff correction is 2 x 2 x f
    % C = repmat([1;-1],[1,2,60]);
    confLimit = B;
end


% function [pvalue, unCorrpvalue] = runPermutationTest(dataObs,dataA,dataB,nPermutation,alpha)
% %RUNPERMUTATIONTEST Description
% %	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
% %
% 		pvalue     = zeros(1,90);
% 		nStanding  = size(dataB,2);
% 		pooledData = cat(2,dataA,dataB);
% 
% 		% we perform a permutation test for each STN separatelly
% 		for permIdx = 1:nPermutation
% 
% 			riffledIndices = randperm(size(pooledData,2));
% 
% 			standPsd = pooledData(:,riffledIndices(1:nStanding),:);
% 			walkPsd	= pooledData(:,riffledIndices(nStanding+1:end),:);
% 			
% 			% compute permutated statistics
% 			dataPerm = (squeeze(mean(walkPsd,2)) - squeeze(mean(standPsd,2)))./squeeze(mean(standPsd,2));
% 		
% 			% compute pvalues for all frequencies and all time points.
% 			pvalue = pvalue + double(dataPerm' > dataObs)./nPermutation;
% 
% 		end
% 		unCorrpvalue = pvalue;
% 		pvalue = fdrCorrection(pvalue,alpha);
% 
% end


% function pvalue = fdrCorrection(pvalue, alpha)
% %FDRCORRECTION Description
% %	PVALUE  = FDRCORRECTION() Long description
% %
% 
% 	tmpPvalue 	= sort(pvalue(:));
% 	N 			= numel(pvalue);
% 	FDR 		= alpha.*(1:N)./N;
% 	thr 		= FDR(find(tmpPvalue <= FDR',1,'last'));
% 	pvalue(pvalue >= thr) = 1;
% 
% end
function crossSpectrum = computeSpectrum(filename,channelData,iChannels,chOrder)
% Description
%	CROSSSPECTRUM = () Long description
%
%
    dataStruct   = in_bst_data(filename);
    channelFlags = dataStruct.ChannelFlag==1;
	[ftData, ~, ~] = out_fieldtrip_data(filename);
% 
% 	if nargin < 4
% 			timeWindow = ftData.time;
%     end
    
    channelIndexes  = zeros(numel(channelData.Channel),1);
    channelIndexes(iChannels) = 1;
    goodChannelMask = channelIndexes & channelFlags;
    
    chLabels        = {channelData.Channel(goodChannelMask).Name};
	tapNW			= 2;
	f 				= 1:60;
	chancomb 		= {channelData.Channel(goodChannelMask).Name};
    
    Fs 				= 1/mean(diff(ftData.time{1}));

	cfg 			= [];
	cfg.output 		='powandcsd';
	cfg.taper 		= 'dpss';
	cfg.channel 	= chLabels;
	cfg.channelcmb 	= chancomb;
	cfg.method  	= 'mtmfft';
	cfg.foi 		= f;
	cfg.pad 		= 'nextpow2';
	cfg.tapsmofrq 	= tapNW*Fs/length(ftData.time{1});
	% CrossSpectrum.powspctr tr x 2 x f  
	% 			   .crssspctr tr x 1 x f 
	% complex values
	[CrossSpectrum] = ft_freqanalysis(cfg, ftData);

	normBand 		= f >= 7 & f<= 60;
% 	normFactor		= mean(abs(CrossSpectrum.crsspctrm(1,normBand)),2);
% 	crossSpectrum   = abs(CrossSpectrum.crsspctrm)./repmat(normFactor,[1 numel(f)]);
 	normFactor		= mean(abs(CrossSpectrum.powspctrm(1,normBand)),2);
 	crossSpectrum   = CrossSpectrum.powspctrm./repmat(normFactor,[1 numel(f)]);
    crossSpectrum   = crossSpectrum(chOrder,:);

end


