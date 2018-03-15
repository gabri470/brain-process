function varargout = process_walking_phaseSynch( varargin )

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
sProcess.Comment     = 'Phase Synch';
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

standingConMask = ~cellfun(@isempty,regexp(conditionStrings,'(s|S)tanding'));
walkingConMask 	= ~cellfun(@isempty,regexp(conditionStrings,'(w|W)alking'));
restingConMask 	= ~cellfun(@isempty,regexp(conditionStrings,'(r|R)esting'));
offConMask 		= ~cellfun(@isempty,regexpi(conditionStrings,'off'));
montageConMask 	= ~cellfun(@isempty,regexpi(conditionStrings,'visite'));

% order patients in descending order based on dopamine-deplition
% patientsOrder = {'wue03','wue09','wue04','wue02','wue10','wue07','wue06','wue11'};
%[~,subjOrder] = ismember(patientsOrder,subjectNames);

f1 = figure('papertype','a4','paperorientation',...
    'portrait','visible','on');


OutputFiles = {};

nPermutation = 3;
nFilters	 = 13;
% alpha		 = 0.05;
trialLength  = 3;

standingData = cell(nSubjects,1);
walkingData = cell(nSubjects,1);
restingData = cell(nSubjects,1);
plotIdx = 1;

bst_progress('start','Estimating Phase Sync','Estimating Phase Sync',0,nSubjects);
for subjIdx = 1:nSubjects
    
    % for each subject separately we pick standing condition
    subjectMask 	= ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));
    
    standingFileIdx = find(subjectMask & standingConMask & offConMask);% & ~montageConMask);
    walkingFileIdx 	= find(subjectMask & walkingConMask & offConMask);% & ~montageConMask);
    restingFileIdx 	= find(subjectMask & restingConMask & offConMask);% & ~montageConMask);
    
    if ~isempty(standingFileIdx)
        % get bad channels
        standChFlag		= getfield(in_bst_data(sInputs(standingFileIdx).FileName,'ChannelFlag'),'ChannelFlag');

        [standingStruct,~,~]  = out_fieldtrip_data(sInputs(standingFileIdx).FileName);
        standChannels 	= in_bst_channel(sInputs(standingFileIdx).ChannelFile);
        standiChannels	= channel_find( standChannels.Channel,'SEEG');    

        [filteredstandingStruct,f] = narrowBandFiltering(standingStruct, nFilters,standiChannels,standChFlag,trialLength);
        % [plv, f, nPLV] = computePhaseMetric
        [standPlv,standnPlv, standAmp] = computePhaseMetric(filteredstandingStruct,standingStruct.label,nPermutation);
        
    else
        standPlv = zeros(nFilters,1);
        standnPlv = zeros(nFilters,1);
        standAmp = zeros(nFilters,1);
    end
    
    standingData{plotIdx,1} = standPlv;
    standingData{plotIdx,2} = standAmp;
    standingData{plotIdx,3} = standnPlv;
    standingData{plotIdx,4} = imag(standPlv);
        
    if ~isempty(restingFileIdx)
        restingStruct = struct('dimord',[],'trial',[],'time',[],'label',[],'elec',[]);

        % resting recordings can be more than 1
        for restIdx = 1:numel(restingFileIdx)
            % read resting time-freq data
            restChFlag 	= getfield(...
                in_bst_data(sInputs(restingFileIdx(restIdx)).FileName,...
                'ChannelFlag'),'ChannelFlag');

            % this contains a single file with multiple trials
            [restIdxStruct,~,~]         = out_fieldtrip_data(sInputs(restingFileIdx(restIdx)).FileName);
            restChannels                = in_bst_channel(sInputs(restingFileIdx(restIdx)).ChannelFile);
            restiChannels               = channel_find(restChannels.Channel,'SEEG');
            [filteredrestingStruct,f]   = narrowBandFiltering(restIdxStruct, nFilters,restiChannels,restChFlag,trialLength);

            if isempty(filteredrestingStruct)
                continue
            end

            if restIdx == 1
                restingStruct = filteredrestingStruct;
            else
                for fIdx = 1:nFilters
                    restingStruct(fIdx).trial = [restingStruct(fIdx).trial filteredrestingStruct(fIdx).trial];
                    restingStruct(fIdx).time  = [restingStruct(fIdx).time filteredrestingStruct(fIdx).time];                   
                    restingStruct(fIdx).hdr.nTrials = numel(restingStruct(fIdx).trial);
                end
            end
        end

        [restPlv, restnPlv, restAmp] = computePhaseMetric(restingStruct,restIdxStruct.label,nPermutation);

    else
        restPlv = zeros(nFilters,1);
        restnPlv = zeros(nFilters,1);
        restAmp = zeros(nFilters,1);
    end
    
    restingData{plotIdx,1} = restPlv;
    restingData{plotIdx,2} = restAmp;
    restingData{plotIdx,3} = restnPlv;
    restingData{plotIdx,4} = imag(restPlv);
    
    if ~isempty(walkingFileIdx)
        walkingStruct   = struct('dimord',[],'trial',[],'time',[],'label',[],'elec',[]);

        for walkIdx = 1:numel(walkingFileIdx)
            % read walking time-freq data
            walkChFlag 	= getfield(...
                in_bst_data(sInputs(walkingFileIdx(walkIdx)).FileName,...
                'ChannelFlag'),'ChannelFlag');

            % this contains a single file with multiple trials
            [walkIdxStruct,~,~]         = out_fieldtrip_data(sInputs(walkingFileIdx(walkIdx)).FileName);
            walkChannels                = in_bst_channel(sInputs(walkingFileIdx(walkIdx)).ChannelFile);
            walkiChannels               = channel_find(walkChannels.Channel,'SEEG');
            [filteredwalkingStruct,f]   = narrowBandFiltering(walkIdxStruct, nFilters,walkiChannels,walkChFlag,trialLength);

            if isempty(filteredwalkingStruct)
                continue
            end

            if walkIdx == 1
                walkingStruct = filteredwalkingStruct;
            else
                for fIdx = 1:nFilters
                    walkingStruct(fIdx).trial = [walkingStruct(fIdx).trial filteredwalkingStruct(fIdx).trial];
                    walkingStruct(fIdx).time  = [walkingStruct(fIdx).time filteredwalkingStruct(fIdx).time];
                    % compute the STN coherence during walking
                    % walkingStruct(fIdx) = rmfield(walkingStruct(fIdx),'sampleinfo');
                    walkingStruct(fIdx).hdr.nTrials = numel(walkingStruct(fIdx).trial);
                end
            end
        end


        [walkPlv, walknPlv,walkAmp] = computePhaseMetric(walkingStruct,walkIdxStruct.label,nPermutation);

        
    else
        walkPlv = zeros(nFilters,1);
        walknPlv = zeros(nFilters,1);
        walkAmp = zeros(nFilters,1);
    end
    
    walkingData{plotIdx,1} = walkPlv;
    walkingData{plotIdx,2} = walkAmp;
    walkingData{plotIdx,3} = walknPlv;
    walkingData{plotIdx,4} = imag(walkPlv);
    
    ax1 = subplot(nSubjects,4,4*(plotIdx-1)+1,'NextPlot','add');
    plot(f,abs(restPlv),'LineWidth',1,'Color',[0 109 219]./255);
    plot(f,abs(standPlv),'LineWidth',1,'Color',[255 109 182]./255);
    plot(f,abs(walkPlv),'LineWidth',1,'Color',[76 255 36]./255);
    %legend({'rest','stand','walk','C.I'});
    xlim([6 60]);
    ylim([0 0.5]);

    xlabel('Hz');
    ylabel('PLV');
    set(ax1,'Parent',f1);
    
    ax2 = subplot(nSubjects,4,4*(plotIdx-1)+2,'NextPlot','add');
    plot(f,abs(imag(restPlv)),'LineWidth',1,'Color',[0 109 219]./255);
    plot(f,abs(imag(standPlv)),'LineWidth',1,'Color',[255 109 182]./255);
    plot(f,abs(imag(walkPlv)),'LineWidth',1,'Color',[76 255 36]./255);
    %legend({'rest','stand','walk'});
    xlim([6 60]);
    ylim([0 0.5]);

    xlabel('Hz');
    ylabel('iPLV');
    set(ax2,'Parent',f1);
    
    ax3 = subplot(nSubjects,4,4*(plotIdx-1)+3,'NextPlot','add');
    plot(f,abs(restAmp),'LineWidth',1,'Color',[0 109 219]./255);
    plot(f,abs(standAmp),'LineWidth',1,'Color',[255 109 182]./255);
    plot(f,abs(walkAmp),'LineWidth',1,'Color',[76 255 36]./255);
    xlim([6 60]);
    ylim([0 0.5]);

    xlabel('Hz');
    ylabel('AAc');
    set(ax3,'Parent',f1);
    
    ax4 = subplot(nSubjects,4,4*(plotIdx-1)+4,'NextPlot','add');
    plot(f,abs(restnPlv),'LineWidth',1,'Color',[0 109 219]./255);
    plot(f,abs(standnPlv),'LineWidth',1,'Color',[255 109 182]./255);
    plot(f,abs(walknPlv),'LineWidth',1,'Color',[76 255 36]./255);
    plot(f,repmat(2.42,size(f)),'--');
    %legend({'rest','stand','walk','C.I'});
    xlim([6 60]);
    %ylim([0 1]);
    
    xlabel('Hz');
    ylabel('nPLV');
    set(ax4,'Parent',f1);
    drawnow
    plotIdx = plotIdx + 1;
    bst_progress('inc',1);
    
end % subject loop

% restingData, walkingData and standingData contain cPLV CC and nPLV for all subjects
% in the form of nSubjects x 3 cells of nFreq x 1 elements 

% reformat data in a 3D matrix with nSubjects x Freq x metric
% first concat all elements in nFreq x ( nSubjects x metric )
restData = reshape(cat(2,restingData{:}),[nFilters nSubjects 4]);
walkData = reshape(cat(2,walkingData{:}),[nFilters nSubjects 4]);
standData = reshape(cat(2,standingData{:}),[nFilters nSubjects 4]);

[confLimRest, avgRest]   = myBootstrap(abs(restData),nSubjects,10);
[confLimWalk, avgWalk]   = myBootstrap(abs(walkData),nSubjects,10);
[confLimStand, avgStand] = myBootstrap(abs(standData),nSubjects,10);

figure('papertype','a4','paperorientation','portrait','Visible','on');

%% compare PLV across conditions
pRestvsWalk  = mySignRank(abs(restData)-abs(walkData),1);
pRestvsStand = mySignRank(abs(restData)-abs(standData),1);

subplot(1,4,1,'NextPlot','add')
plot(f,avgRest(:,1,1),'LineWidth',1,'Color',[0 109 219]./255);
plot(f,avgStand(:,1,1),'LineWidth',1,'Color',[255 109 182]./255);
plot(f,avgWalk(:,1,1),'LineWidth',1,'Color',[76 255 36]./255);
legend({'rest' 'stand' 'walk'});

plot(f,(pRestvsWalk)*0.3,'r','MarkerSize',1);
plot(f,(pRestvsStand)*0.3,'k','MarkerSize',1);
fill_between(f,confLimRest(1,:,1),confLimRest(2,:,1),f,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimStand(1,:,1),confLimStand(2,:,1),f,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimWalk(1,:,1),confLimWalk(2,:,1),f,'FaceColor',[76 255 36]./255,'FaceAlpha',0.2,'EdgeColor','None');
xlabel('Freq. (Hz)');
ylabel('PLV');
xlim([0 60]);


%% compare CC across conditions
pRestvsWalk = mySignRank(abs(restData)-abs(walkData),2);
pRestvsStand= mySignRank(abs(restData)-abs(standData),2);
subplot(1,4,2,'NextPlot','add')
plot(f,avgRest(:,1,2),'LineWidth',1,'Color',[0 109 219]./255);
plot(f,avgStand(:,1,2),'LineWidth',1,'Color',[255 109 182]./255);
plot(f,avgWalk(:,1,2),'LineWidth',1,'Color',[76 255 36]./255);
legend({'rest' 'stand' 'walk'});
plot(f,(pRestvsWalk)*0,'r','MarkerSize',1);
plot(f,(pRestvsStand)*0,'k','MarkerSize',1);
fill_between(f,confLimRest(1,:,2),confLimRest(2,:,2),f,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimStand(1,:,2),confLimStand(2,:,2),f,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimWalk(1,:,2),confLimWalk(2,:,2),f,'FaceColor',[76 255 36]./255,'FaceAlpha',0.2,'EdgeColor','None');
xlabel('Freq. (Hz)');
ylabel('CC');
xlim([0 60]);

%% compare nPLV across conditions
pRestvsWalk = mySignRank(abs((restData))-abs((walkData)),3);
pRestvsStand= mySignRank(abs((restData))-abs((standData)),3);
subplot(1,4,3,'NextPlot','add')
plot(f,avgRest(:,1,3),'LineWidth',1,'Color',[0 109 219]./255);
plot(f,avgStand(:,1,3),'LineWidth',1,'Color',[255 109 182]./255);
plot(f,avgWalk(:,1,3),'LineWidth',1,'Color',[76 255 36]./255);
plot(f,repmat(2.42,size(f)),'--');
legend({'rest' 'stand' 'walk'});
plot(f,(pRestvsWalk)*120,'r','MarkerSize',1);
plot(f,(pRestvsStand)*120,'k','MarkerSize',1);
fill_between(f,confLimRest(1,:,3),confLimRest(2,:,3),f,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimStand(1,:,3),confLimStand(2,:,3),f,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimWalk(1,:,3),confLimWalk(2,:,3),f,'FaceColor',[76 255 36]./255,'FaceAlpha',0.2,'EdgeColor','None');
xlabel('Freq. (Hz)');
ylabel('nPLV');
xlim([0 60]);

%% compare iPLV acorss conditions
pRestvsWalk = mySignRank(abs((restData))-abs((walkData)),4);
pRestvsStand= mySignRank(abs((restData))-abs((standData)),4);
subplot(1,4,4,'NextPlot','add')
plot(f,avgRest(:,1,4),'LineWidth',1,'Color',[0 109 219]./255);
plot(f,avgStand(:,1,4),'LineWidth',1,'Color',[255 109 182]./255);
plot(f,avgWalk(:,1,4),'LineWidth',1,'Color',[76 255 36]./255);

legend({'rest' 'stand' 'walk'});
plot(f,(pRestvsWalk)*0.3,'r','MarkerSize',1);
plot(f,(pRestvsStand)*0.3,'k','MarkerSize',1);
fill_between(f,confLimRest(1,:,4),confLimRest(2,:,4),f,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimStand(1,:,4),confLimStand(2,:,4),f,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(f,confLimWalk(1,:,4),confLimWalk(2,:,4),f,'FaceColor',[76 255 36]./255,'FaceAlpha',0.2,'EdgeColor','None');
xlabel('Freq. (Hz)');
ylabel('iPLV');
xlim([0 60]);

%% Compute delta plv in conditions for selected F bands
bandMasks = [f >= 4  & f < 8, f >= 8  & f < 13, f >= 13 & f < 30, f >= 30 & f < 60];

restingData = abs(avgInFBands(restData,bandMasks));
standingData = abs(avgInFBands(standData,bandMasks));
walkingData = abs(avgInFBands(walkData,bandMasks));

deltaPlvMean(1,:) = mean( (restingData - standingData)./standingData );
deltaPlvMean(2,:) = mean( (restingData - walkingData)./walkingData );
deltaPlvMean(3,:) = mean( (standingData - walkingData)./walkingData);

deltaCohStd(1,:) = std( (restingData - standingData)./standingData  )./sqrt(nSubjects);
deltaCohStd(2,:) = std( (restingData - walkingData)./walkingData )./sqrt(nSubjects);
deltaCohStd(3,:) = std( (standingData - walkingData)./walkingData )./sqrt(nSubjects);

[nGroups, nBars] = size(deltaPlvMean); 
groupWidth = min(0.8, nBars/(nBars+1.5));
figure, bar(deltaPlvMean);
hold on
for i = 1:nBars
      % magic numbers ... 
      x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);  % Aligning error bar with individual bar
      errorbar(x, deltaPlvMean(:,i), deltaCohStd(:,i), 'k', 'linestyle', 'none');

end

legend({'theta','alpha','beta','gamma'})
set(gca,'XTick',[1 2 3],'XTickLabel',{'rest-stand','rest-walk','stand-walk'});

end % function

function avg = avgInFBands(data,mask)
    
    %% data are stored in a f x sub x metric
    [~,nSubjects,~]=size(data);
    nFreq = size(mask,2);
    avg = nan(nSubjects,nFreq);
    for fIdx = 1:nFreq
        avg(:,fIdx) = mean(data(mask(:,fIdx),:,1),1);
    end

end
function [confLimit,dataMean] = myBootstrap(data,nSubject,nBootstraps)

    % data matrix has  f x nSubjects x metric
    % bootstrap the C.L. for mean
    bootstrapIndexes = randi(nSubject,nBootstraps,nSubject);
    [nFilters,~,nMetric] = size(data);
    
    %dataMean will be a 13 x 1 x 3
    dataMean = mean(data,2);
    
    % currBootstraps will contain nBootstraps x nStn x f
    currBootstraps = nan(nBootstraps,nFilters,nMetric);
    
    for idx = 1:nBootstraps
        currBootstraps(idx,:,:) = squeeze(mean(data(:,bootstrapIndexes(idx,:),:),2));
    end
    
    % confLimit will contain [UB LB] x nStn x f
    confLimit = prctile(currBootstraps,[5 95]);
        
end

function h = mySignRank(data,metricIdx)

    % data f x nSubjects x metric
    [nFreq,~,~] = size(data);
    p = ones(1,nFreq);
    h = ones(1,nFreq);
    
    for fIdx = 1:nFreq
        p(fIdx) = signrank(squeeze(data(fIdx,:,metricIdx)));
    end
    h(p >= 0.05) = NaN;
    
end


% function pvalue = fdrCorrection(pvalue, alpha)
% %FDRCORRECTION Description
% %	PVALUE  = FDRCORRECTION() Long description
% %
% 
% tmpPvalue 	= sort(pvalue(:));
% N 			= numel(pvalue);
% FDR 		= alpha.*(1:N)./N;
% thr 		= FDR(find(tmpPvalue <= FDR',1,'last'));
% 
% if isempty(thr)
%     pvalue  = ones(size(thr));
% else
%     pvalue(pvalue >= thr) = 1;
% end
% 
% end



function [plv, nPLV, amp] = computePhaseMetric(data,~,nPermutation)
% Description
%	[PTE,PLV, F] = computePhaseMetric(data,channelNames)


%chMask	 = ~cellfun(@isempty,regexp(channelNames,'E[0-9]*'));
nFilters = numel(data);
nTrials  = numel(data(1).trial);
plv 	 = complex(zeros(nFilters,nTrials),zeros(nFilters,nTrials));
amp      = zeros(nFilters,nTrials);

plvSurr	 = complex(zeros(nFilters,nTrials,nPermutation),zeros(nFilters,nTrials,nPermutation));

nChans   = numel(unique([data.label]));
if nChans > 2
    error('casino');
end

for fIdx = 1:nFilters
    
    % check if any trial contains NaN values and discard it
    trlMask    = cellfun(@(x) sum(isnan(x),2),data(fIdx).trial,'uni',false);
    trlMask    =  reshape(cat(1,trlMask{:}),2,numel(data(fIdx).trial));
    trlMask    = sum(trlMask,1) == 0;
    nTrials    = numel(trlMask);
    trlIndices = 1:nTrials;
    
    for trialIdx = trlIndices(trlMask)
        
        % function [pTE, plv, plvSurr] = phaseTE(Xf,lag,nPermutation)        
        [tmp,tmpSurr,tmpAmp] = phaseAmplitudeCorrelation(data(fIdx).trial{trialIdx},nPermutation);


        amp(fIdx,trialIdx) = tmpAmp;
        plv(fIdx,trialIdx) = tmp;
        plvSurr(fIdx,trialIdx,:) = tmpSurr;
        
    end
    
end
% get the mean across trials for each filter
plv 	= squeeze(nanmean(plv,2));
amp     = squeeze(nanmean(amp,2));

% get the mean across trials and across permutations for each filter
plvSurr = squeeze(mean(mean(plvSurr,2),3));

% this normalized PLV
nPLV 	= abs(plv)./abs(plvSurr);

end

function [plv,plvSurr,ampCorr,ampCorrSurr] = phaseAmplitudeCorrelation(X,nPermutation)
% Description
%	[PLV,PLVSURR,AMPCORR,AMPCORRSURR] = phaseAmplitudeCorrelation(X,nPermutation)
%
%%% Inst. phase
Xfc = hilbert(X')';
Xfh = angle(Xfc);
amp = abs(Xfc);

ampCorr = corrcoef(amp');
ampCorr = ampCorr(1,2);

phi=diff(Xfh);
plv=mean(exp(1i*phi),2);
plvSurr = zeros(1,nPermutation);
ampCorrSurr = zeros(1,nPermutation);

for pIdx = 1:nPermutation
    % compute surrogate
    offset              = randi(size(Xfc,2),1);
    Xfc(2,:)            = circshift(Xfc(2,:),offset);
    
    phi                 = diff(angle(Xfc));
    plvSurr(pIdx)       = mean(exp(1i*phi),2);
    tmp                 = corrcoef(abs(Xfc)');
    ampCorrSurr(pIdx)   = tmp(1,2);
end

end



function data = preproc(data,~,chFlag,trialLength)
%	preproc
%	data = preproc() split data in trials
%

channelNames = data.label;
iChannels    = find(~cellfun(@isempty,regexp(channelNames,'E[0-9]+')));

if numel(iChannels) > 2 && numel(iChannels) == numel(chFlag)
    iChannels = iChannels(chFlag==1);
elseif any(chFlag(iChannels)==-1)
    error('bad channels included');
end
cfg = [];

% at this point the recordings are just a single
% continuos stream of samples.
begRecording 	= min(data.time{1});
endRecording 	= max(data.time{1});

fs 				= 1/mean(diff(data.time{1}));

% we should split this in ntrials of 3s
totLength 		= endRecording-begRecording;
nTrials   		= floor(totLength/trialLength);
offset 			= totLength/2;
analysisWindow  = nTrials*trialLength;
startTime 		= offset - analysisWindow/2;
endTime   		= analysisWindow/2+offset;


trials			= cat(2,(startTime:trialLength:endTime-trialLength)',...
    (startTime+trialLength:trialLength:endTime)',zeros(nTrials,1));

trials			= floor(trials*fs)+1;
cfg.trl			= trials;

data			= ft_redefinetrial(cfg,data);

cfg 			= [];
cfg.continuous  = 'yes';
cfg.channel		= data.label(iChannels);
cfg.detrend		= 'yes';
data 			= ft_preprocessing(cfg,data);

end


function [dataFiltered, f] = narrowBandFiltering(data,nFilters,iChannels,chFlag,trialLength)
% Description
%	DATAFILTERED = () Long description
%

% nominal frequency and central frequency
fn = 2;
% band-flat top and band pass width
wb = 0.5;
% stop band
ws = 2;
% multiplier
m = sqrt(2);

f = zeros(nFilters,1);
dataOrig = data;
fs = 1/mean(diff(dataOrig.time{1}));  

for fIdx = 1:nFilters
    
    passBandLp = fn + (wb * fn)/2;
    passBandHp = fn - (wb * fn)/2;
    stopBandLp = min([0.95*(fs/2) fn * ws]);
    stopBandHp = fn / ws;
    
    f(fIdx) = fn;
    fn = fn * m;
    
    cfg.bpfilter = 'yes';
    
    cfg.bpfilttype = 'fir';
    cfg.bpfiltdir = 'twopass';
    
    
    fcaz = [stopBandHp,passBandHp,passBandLp,stopBandLp]./(fs/2);
    mags = [0 1 0];
    devs = [1e-3 1e-3 1e-3];
    
    [n,Wn,beta,ftype] = kaiserord(fcaz, mags, devs, 2);
    n = n + rem(n,2);  % ensure even order
    cfg.bpfreq   = [ passBandHp passBandLp];
    cfg.bpfiltord = n;
    cfg.bpfiltwintype = kaiser(n+1,beta);
    %cfg.bpfiltdev = devs;
    
    cfg.padding  = trialLength;
    cfg.continuous = 'yes';
    data = ft_preprocessing(cfg,dataOrig);
    
  
%     
%     [fftmp,h,ff] = ft_preproc_bandpassfilter(dataOrig.trial{1},Fs,[passBandHp passBandLp],4,'but','twopass','reduce',3,'hamming',0.001,'yes');
%     plot(ff,abs(h)); hold on
    dataFiltered(fIdx) = preproc(data,iChannels,chFlag,trialLength);
%     
end

end
