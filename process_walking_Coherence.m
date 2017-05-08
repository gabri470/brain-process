function varargout = process_walking_Coherence( varargin )

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
sProcess.Comment     = 'Coherence';
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

%	DATA_FOLDER = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','info');

%	sideFile = fullfile(DATA_FOLDER,'patientSides.csv');
%	[nameSubjects, mostAffSides] = textread(sideFile,'%s %s\n','delimiter',',');

subjectNames = unique({sInputs.SubjectName});
nSubjects = numel(subjectNames);

%	mostAffSides = mostAffSides(ismember(nameSubjects,subjectNames));

% we need to find files that have to be processed and group
% them for subjects and conditions
% First group all sInput.Comment together
conditionStrings 	= {sInputs.Condition};

standingConMask 	= ~cellfun(@isempty,regexp(conditionStrings,'(s|S)tanding'));
walkingConMask 		= ~cellfun(@isempty,regexp(conditionStrings,'(w|W)alking'));
restingConMask 		= ~cellfun(@isempty,regexp(conditionStrings,'(r|R)esting'));

f1 = figure('papertype','a4','paperorientation',...
    'portrait','visible','on');

f2 = figure('papertype','a4','paperorientation',...
    'portrait','Visible','on');

f3 = figure('papertype','a4','paperorientation',...
    'portrait','Visible','on');

f4 = figure('papertype','a4','paperorientation',...
    'portrait','Visible','on');

OutputFiles = {};

nPermutation = 10;
%nFilters     = 12;
alpha		 = 0.05;

standingData = zeros(nSubjects,4);
restingData = zeros(nSubjects,4);
walkingData = zeros(nSubjects,4);

grpRestingData = zeros(nSubjects,119);
grpStandingData= zeros(nSubjects,119);
grpWalkingData = zeros(nSubjects,119);
    

for subjIdx = 1:nSubjects
    % for each subject separately we pick standing condition
    subjectMask 		= ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));
    
    standingFileIdx = find(subjectMask & standingConMask);
    walkingFileIdx 	= find(subjectMask & walkingConMask);
    restingFileIdx 	= find(subjectMask & restingConMask);
    
    % get bad channels
    standChFlag		= getfield(in_bst_data(sInputs(standingFileIdx).FileName,'ChannelFlag'),'ChannelFlag');
    
    [standingStruct,~,~]  = out_fieldtrip_data(sInputs(standingFileIdx).FileName);
    standChannels 	= in_bst_channel(sInputs(standingFileIdx).ChannelFile);
    standiChannels	= channel_find( standChannels.Channel,'SEEG');
    
    
    standingStruct  = preproc(standingStruct,standiChannels,standChFlag,3);
    % we should here quantify the STN coherence
    standCoh		= computeCoherence(standingStruct,{standChannels.Channel(standiChannels).Name});
    [~,standPvalue,standAvgSurr,~] = runPermutationTest(standCoh,standingStruct,nPermutation,alpha);
    
    %			filteredstandingStruct = narrowBandFiltering(standingStruct, nFilters,standiChannels,standChFlag,3);
    %			[standPlv, f] = computePhaseMetric(filteredstandingStruct,standingStruct.label,nPermutation);
    restingStruct = struct('dimord',[],'trial',[],'time',[],'label',[],'elec',[]);
    
    % resting recordings can be more than 1
    for restIdx = 1:numel(restingFileIdx)
        % read resting time-freq data
        restChFlag 	= getfield(...
            in_bst_data(sInputs(restingFileIdx(restIdx)).FileName,...
            'ChannelFlag'),'ChannelFlag');
        
        
        [restIdxStruct,~,~] = out_fieldtrip_data(sInputs(restingFileIdx(restIdx)).FileName);
        restChannels 		= in_bst_channel(sInputs(restingFileIdx(restIdx)).ChannelFile);
        restiChannels		= channel_find(restChannels.Channel,'SEEG');
        restIdxStruct       = preproc(restIdxStruct,restiChannels,restChFlag,3);
        
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

    
    [~,restPvalue, restAvgSurr,~] = runPermutationTest(restCoh,restingStruct,nPermutation,alpha);
    
    walkingStruct   = struct('dimord',[],'trial',[],'time',[],'label',[],'elec',[]);
    
    for walkIdx = 1:numel(walkingFileIdx)
        walkChFlag 	= getfield(...
            in_bst_data(sInputs(walkingFileIdx(walkIdx)).FileName,...
            'ChannelFlag'),'ChannelFlag');
        
        [walkIdxStruct,~,~]  = out_fieldtrip_data(sInputs(walkingFileIdx(walkIdx)).FileName);
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
    [~,walkPvalue, walkAvgSurr,~] = runPermutationTest(walkCoh,walkingStruct,nPermutation,alpha);
    
    f = 1:.5:60;
    
    thetaBand = f >= 4  & f < 8;
    alphaBand = f >= 8  & f < 13;
    betaBand  = f >= 13 & f < 30;
    gammaBand = f >= 30 & f < 60;
    
    standingData(subjIdx,1) = mean(abs(standCoh.cohspctrm(thetaBand)));
    standingData(subjIdx,2) = mean(abs(standCoh.cohspctrm(alphaBand)));
    standingData(subjIdx,3) = mean(abs(standCoh.cohspctrm(betaBand)));
    standingData(subjIdx,4) = mean(abs(standCoh.cohspctrm(gammaBand)));
    
    restingData(subjIdx,1) = mean(abs(restCoh.cohspctrm(thetaBand)));
    restingData(subjIdx,2) = mean(abs(restCoh.cohspctrm(alphaBand)));
    restingData(subjIdx,3) = mean(abs(restCoh.cohspctrm(betaBand)));
    restingData(subjIdx,4) = mean(abs(restCoh.cohspctrm(gammaBand)));
    
    walkingData(subjIdx,1) = mean(abs(walkCoh.cohspctrm(thetaBand)));
    walkingData(subjIdx,2) = mean(abs(walkCoh.cohspctrm(alphaBand)));
    walkingData(subjIdx,3) = mean(abs(walkCoh.cohspctrm(betaBand)));
    walkingData(subjIdx,4) = mean(abs(walkCoh.cohspctrm(gammaBand)));
    
    grpRestingData(subjIdx,:) = abs(restCoh.cohspctrm);
    grpStandingData(subjIdx,:)= abs(standCoh.cohspctrm);
    grpWalkingData(subjIdx,:) = abs(walkCoh.cohspctrm);
   
    ax1 = subplot(2,4,subjIdx,'NextPlot','add');
    plot(restCoh.freq,abs(restCoh.cohspctrm),'LineWidth',1);
    plot(standCoh.freq,abs(standCoh.cohspctrm),'LineWidth',1);
    plot(walkCoh.freq,abs(walkCoh.cohspctrm),'LineWidth',1);
    legend({'rest','stand','walk'});
    xlim([6 60]);
    ylim([0 1]);
    title(subjectNames(subjIdx));
    xlabel('Hz');
    ylabel('Coh');
    set(ax1,'Parent',f1);
    
    ax2 = subplot(2,4,subjIdx,'NextPlot','add');
    plot(standCoh.freq,abs(standCoh.cohspctrm),'LineWidth',1);
    standCoh.cohspctrm( standPvalue >= 0.05 ) = NaN;
    plot(standCoh.freq,abs(standCoh.cohspctrm),'LineWidth',2);
    plot(standCoh.freq,abs(standAvgSurr),'--');

    xlim([6 60]);
    ylim([0 1]);
    title(subjectNames(subjIdx));
    xlabel('Hz');
    ylabel('Coh');
    set(ax2,'Parent',f2);
    
    ax3 = subplot(2,4,subjIdx,'NextPlot','add');

    plot(walkCoh.freq,abs(walkCoh.cohspctrm),'LineWidth',1);
    walkCoh.cohspctrm( walkPvalue >= 0.05) = NaN;
    plot(walkCoh.freq,abs(walkCoh.cohspctrm),'LineWidth',2);
    plot(walkCoh.freq,abs(walkAvgSurr),'--');

    xlim([6 60]);
    ylim([0 1]);
    title(subjectNames(subjIdx));
    xlabel('Hz');
    ylabel('Coh');
    set(ax3,'Parent',f3);
    
    ax4 = subplot(2,4,subjIdx,'NextPlot','add');
    plot(restCoh.freq,abs(restCoh.cohspctrm),'LineWidth',1);
    restCoh.cohspctrm( restPvalue >= 0.05 ) = NaN;
    plot(restCoh.freq,abs(restCoh.cohspctrm),'LineWidth',2);
    plot(restCoh.freq,abs(restAvgSurr),'--');

    xlim([6 60]);
    ylim([0 1]);
    title(subjectNames(subjIdx));
    xlabel('Hz');
    ylabel('Coh');
    set(ax4,'Parent',f4);
    
end % subject loop

% mean in bands as PLV and CC to have similar freq resolution
fn = 2;
% band-flat top and band pass width
wb = 0.5;
% stop band
ws = 2;
% multiplier
m = sqrt(2);
passBands = [];
fBands = zeros(11,1);
f = 1:.5:60;
for fIdx = 1:11
    
    passBands(fIdx,:) = [fn-(wb*fn)/2  fn+(wb * fn)/2];
    fBands(fIdx) = fn;
    fn = fn * m;
end
[~, hpMask] = meshgrid(f,passBands(:,1));
[fMask, lpMask] = meshgrid(f,passBands(:,2));
passBandsMask = (fMask >= hpMask & fMask <= lpMask);

grpAvgRestingData = avgInfBands(grpRestingData,passBandsMask);
grpAvgStandingData = avgInfBands(grpStandingData,passBandsMask);
grpAvgWalkingData = avgInfBands(grpWalkingData,passBandsMask);

[grpRestingDataCL, grpRestingDataMean] = myBootstrap(grpAvgRestingData,nSubjects,10,passBandsMask);
[grpWalkingDataCL, grpWalkingDataMean] = myBootstrap(grpAvgWalkingData,nSubjects,10,passBandsMask);
[grpStandingDataCL, grpStandingDataMean] = myBootstrap(grpAvgStandingData,nSubjects,10,passBandsMask);

figure
f = fBands;

pRestvsWalk  = mySignRank(grpAvgRestingData-grpAvgWalkingData);
pRestvsStand = mySignRank(grpAvgRestingData-grpAvgStandingData);

subplot(1,1,1,'NextPlot','add')
plot(f,grpRestingDataMean,'LineWidth',1,'Color',[0 109 219]./255);
plot(f,grpStandingDataMean,'LineWidth',1,'Color',[255 109 182]./255);
plot(f,grpWalkingDataMean,'LineWidth',1,'Color',[76 255 36]./255);

legend({'rest' 'stand' 'walk'});
plot(f,(pRestvsWalk)*0.5,'r*','MarkerSize',1);
plot(f,(pRestvsStand)*0.5,'k*','MarkerSize',1);
fill_between(f,grpRestingDataCL(1,:),grpRestingDataCL(2,:),f,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,grpStandingDataCL(1,:),grpStandingDataCL(2,:),f,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,grpWalkingDataCL(1,:),grpWalkingDataCL(2,:),f,'FaceColor',[76 255 36]./255,'FaceAlpha',0.2,'EdgeColor','None');
xlabel('Freq. (Hz)');
ylabel('Coh');
xlim([0 60]);

deltaCohMean(1,:) = mean( (restingData - standingData)./standingData );
deltaCohMean(2,:) = mean( (restingData - walkingData)./walkingData );
deltaCohMean(3,:) = mean( (standingData - walkingData)./walkingData);

deltaCohStd(1,:) = std( (restingData - standingData)./standingData  )./sqrt(nSubjects);
deltaCohStd(2,:) = std( (restingData - walkingData)./walkingData )./sqrt(nSubjects);
deltaCohStd(3,:) = std( (standingData - walkingData)./walkingData )./sqrt(nSubjects);

[nGroups, nBars] = size(deltaCohMean); 
groupWidth = min(0.8, nBars/(nBars+1.5));
figure, bar(deltaCohMean);
hold on
for i = 1:nBars
      % magic numbers ... 
      x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);  % Aligning error bar with individual bar
      errorbar(x, deltaCohMean(:,i), deltaCohStd(:,i), 'k', 'linestyle', 'none');

end

legend({'theta','alpha','beta','gamma'})
set(gca,'XTick',[1 2 3],'XTickLabel',{'rest-stand','rest-walk','stand-walk'});


end % function

function avg = avgInfBands(data,bandMasks)
    nSubjects = size(data,1);
    nFreq = size(bandMasks,1);
    avg = nan(nSubjects,nFreq);
    for fIdx = 1:nFreq
        avg(:,fIdx) = squeeze(mean(data(:,bandMasks(fIdx,:)),2));
    end
end

function h = mySignRank(data)

    % data nSubjects x f
    [~,nFreq] = size(data);
    p = ones(1,nFreq);
    h = ones(1,nFreq);
    
    for fIdx = 1:nFreq
        p(fIdx) = signrank(squeeze(data(:,fIdx)));
    end
    h(p >= 0.05) = NaN;
    
end

function [confLimit,dataMean] = myBootstrap(data,nSubject,nBootstraps,passBands)

    % data matrix has  nSubjects x f
    % bootstrap the C.L. for mean
    bootstrapIndexes = randi(nSubject,nBootstraps,nSubject);
    [nFreq,~] = size(passBands);
    
    %dataMean will be a 1 x 13
    dataMean = mean(data);
    
    
    % currBootstraps will contain nBootstraps x f
    currBootstraps = nan(nBootstraps,nFreq);
    for idx = 1:nBootstraps
       
        currBootstraps(idx,:) = squeeze(mean(data(bootstrapIndexes(idx,:),:)));
        
    end
    
    % confLimit will contain [UB LB] x nStn x f
    confLimit = prctile(currBootstraps,[5 95]);
        
end

function [pvalue, unCorrpvalue, avgSurrogate,stdSurrogate] = runPermutationTest(dataObs, data, nPermutation,alpha)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%

dataObs 	 = abs(dataObs.cohspctrm);

pvalue  	 = zeros(size(dataObs));
avgSurrogate = zeros(size(dataObs));
avgSurrSquare = zeros(size(dataObs));

% we perform a permutation test
for permIdx = 1:nPermutation
    for trialIdx = 1:numel(data.trial)
        
        dat = data.trial{trialIdx};
        splitOffset = randi(size(dat,2),1);
        % probably this can be done also using circshift?
        dat(2,:) = [dat(2,splitOffset:end) dat(2,1:splitOffset-1)];
        data.trial{trialIdx} = dat;
        
    end
    
    % compute permutated statistics
    dataPerm 			= computeCoherence(data,data.label);
    dataPerm 			= abs(dataPerm.cohspctrm);
    % x
    avgSurrogate 	= avgSurrogate + dataPerm;
    % x^2
    avgSurrSquare = avgSurrSquare + dataPerm.^2;
    
    % compute pvalues for all frequencies and all time points.
    pvalue = pvalue + double((abs(dataPerm) >= abs(dataObs)))./nPermutation;
    
end

avgSurrogate = avgSurrogate ./ nPermutation;
stdSurrogate = sqrt((avgSurrSquare./nPermutation) - avgSurrogate);
unCorrpvalue = pvalue;
pvalue 		 = fdrCorrection(pvalue,alpha);

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
    pvalue  = ones(size(thr));
else
    pvalue(pvalue >= thr) = 1;
end

end

function coh = computeCoherence(data,channelNames)

cfg = [];

% check if any trial contains NaN values and discard it
trlMask     = cellfun(@(x) sum(isnan(x),2),data.trial,'uni',false);

trlMask     = reshape(cat(1,trlMask{:}),2,numel(data.trial));

trlMask     = sum(trlMask,1) == 0;
cfg.trials	= trlMask;

fs 			= 1/mean(diff(data.time{1}));
cfg.channel = channelNames;
cfg.channelcmb 	= channelNames';
cfg.output 	='powandcsd';
cfg.taper 	= 'dpss';
tapNW		= 2;
cfg.keeptrials 	= 'yes';
cfg.method  = 'mtmfft';
cfg.foi 	= 1:.5:60;
cfg.pad 	= 'nextpow2';
cfg.tapsmofrq 	= tapNW*fs/length(data.time{1});

freq		= ft_freqanalysis(cfg, data);

cfg.method 	= 'coh';
cfg.complex = 'complex';

coh 		= ft_connectivityanalysis(cfg,freq);

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
begRecording = min(data.time{1});
endRecording = max(data.time{1});

fs 			 = 1/mean(diff(data.time{1}));

% we should split this in ntrials of 3s
totLength       = endRecording-begRecording;
nTrials         = floor(totLength/trialLength);
offset          = totLength/2;
analysisWindow	= nTrials*trialLength;
startTime 		= offset - analysisWindow/2;
endTime   		= analysisWindow/2+offset;;


trials			= cat(2,(startTime:trialLength:endTime-trialLength)',...
    (startTime+trialLength:trialLength:endTime)',zeros(nTrials,1));

trials			= floor(trials*fs)+1;
cfg.trl			= trials;

data            = ft_redefinetrial(cfg,data);

cfg 						= [];
cfg.continuous  = 'yes';
cfg.channel		=	data.label(iChannels);
cfg.detrend		= 'yes';
data 			= ft_preprocessing(cfg,data);

end

