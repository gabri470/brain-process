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
nFilters	 = 12;
alpha		 = 0.05;

for subjIdx = 1:nSubjects
    % for each subject separately we pick standing condition
    subjectMask 		= ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));
    
    standingFileIdx = find(subjectMask & standingConMask);
    walkingFileIdx 	= find(subjectMask & walkingConMask);
    restingFileIdx 	= find(subjectMask & restingConMask);
    
    % get bad channels
    standChFlag			= getfield(in_bst_data(sInputs(standingFileIdx).FileName,'ChannelFlag'),'ChannelFlag');
    
    [standingStruct,~,~]  = out_fieldtrip_data(sInputs(standingFileIdx).FileName);
    standChannels 	= in_bst_channel(sInputs(standingFileIdx).ChannelFile);
    standiChannels	= channel_find( standChannels.Channel,'SEEG');
    
    
   
    % we should here quantify the STN coherence    
    [filteredstandingStruct,f] = narrowBandFiltering(standingStruct, nFilters,standiChannels,standChFlag,3);
    % [plv, f, nPLV] = computePhaseMetric
    [standPlv,standnPLV] = computePhaseMetric(filteredstandingStruct,standingStruct.label,nPermutation);
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
        [filteredrestingStruct,~]   = narrowBandFiltering(restIdxStruct, nFilters,restiChannels,restChFlag,3);
        
        if isempty(filteredrestingStruct)
            continue
        end
        
        if restIdx == 1
            restingStruct = filteredrestingStruct;
        else
            for fIdx = 1:nFilters
                restingStruct(fIdx).trial = [restingStruct(fIdx).trial filteredrestingStruct(fIdx).trial];
                restingStruct(fIdx).time  = [restingStruct(fIdx).time filteredrestingStruct(fIdx).time];
                % compute the STN coherence during resting
                % restingStruct(fIdx) = rmfield(restingStruct(fIdx),'sampleinfo');
                restingStruct(fIdx).hdr.nTrials = numel(restingStruct(fIdx).trial);
            end
        end
    end
   
    [restPLV, restnPLV] = computePhaseMetric(restingStruct,restIdxStruct.label,nPermutation);
    
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
        [filteredwalkingStruct,~]   = narrowBandFiltering(walkIdxStruct, nFilters,walkiChannels,walkChFlag,3);
        
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
    

    [walkPlv, walknPLV] = computePhaseMetric(walkingStruct,walkIdxStruct.label,nPermutation);
    
    ax1 = subplot(2,4,subjIdx,'NextPlot','add');
    plot(f,abs(restPLV),'LineWidth',1);
    plot(f,abs(standPlv),'LineWidth',1);
    plot(f,abs(walkPlv),'LineWidth',1);
    legend({'rest','stand','walk'});
    xlim([6 60]);
    ylim([0 1]);
    title(subjectNames(subjIdx));
    xlabel('Hz');
    ylabel('Coh');
    set(ax1,'Parent',f1);
    
    ax2 = subplot(2,4,subjIdx,'NextPlot','add');
    plot(f,abs(imag(restPLV)),'LineWidth',1);
    plot(f,abs(imag(standPlv)),'LineWidth',1);
    plot(f,abs(imag(walkPlv)),'LineWidth',1);
    legend({'rest','stand','walk'});
    xlim([6 60]);
    ylim([0 1]);
    title(subjectNames(subjIdx));
    xlabel('Hz');
    ylabel('Coh');
    set(ax2,'Parent',f2);
%     
%     ax3 = subplot(2,4,subjIdx,'NextPlot','add');
%     %		  plot(restCoh.freq,abs(imag(restCoh.cohspctrm)),'LineWidth',2);
%     %		  plot(standCoh.freq,abs(imag(standCoh.cohspctrm)),'LineWidth',2);
%     plot(walkCoh.freq,abs(walkCoh.cohspctrm),'LineWidth',1);
%     walkCoh.cohspctrm( walkPvalue >= 0.05) = NaN;
%     plot(walkCoh.freq,abs(walkCoh.cohspctrm),'LineWidth',2);
%     plot(walkCoh.freq,abs(walkAvgSurr),'--');
%     %		  legend({'rest','stand','walk'});
%     xlim([6 60]);
%     ylim([0 1]);
%     title(subjectNames(subjIdx));
%     xlabel('Hz');
%     ylabel('Coh');
%     set(ax3,'Parent',f3);
%     
%     ax4 = subplot(2,4,subjIdx,'NextPlot','add');
%     plot(restCoh.freq,abs(restCoh.cohspctrm),'LineWidth',1);
%     restCoh.cohspctrm( restPvalue >= 0.05 ) = NaN;
%     plot(restCoh.freq,abs(restCoh.cohspctrm),'LineWidth',2);
%     plot(restCoh.freq,abs(restAvgSurr),'--');
%     %		  plot(standCoh.freq,abs(standCoh.cohspctrm),'LineWidth',2);
%     %		  plot(walkCoh.freq,abs(walkCoh.cohspctrm),'LineWidth',2);
%     %     legend({'rest','stand','walk'});
%     xlim([6 60]);
%     ylim([0 1]);
%     title(subjectNames(subjIdx));
%     xlabel('Hz');
%     ylabel('Coh');
%     set(ax4,'Parent',f4);
    
end % subject loop
%
% print(f1,'/home/gabri/Dropbox/Isaias_group/walking/all.ps','-dpsc2');
% print(f2,'/home/gabri/Dropbox/Isaias_group/walking/standing.ps','-dpsc2 ');
% print(f3,'/home/gabri/Dropbox/Isaias_group/walking/walking.ps','-dpsc2 ');
% print(f4,'/home/gabri/Dropbox/Isaias_group/walking/resting.ps','-dpsc2 ');

clearvars walkData standData
end % function


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
    dataPerm 	= computePhaseMetric(data,data.label);
    dataPerm 	= abs(dataPerm.cohspctrm);
    % x
    avgSurrogate = avgSurrogate + dataPerm;
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
N 			= numel(pvalue);
FDR 		= alpha.*(1:N)./N;
thr 		= FDR(find(tmpPvalue <= FDR',1,'last'));

if isempty(thr)
    pvalue  = ones(size(thr));
else
    pvalue(pvalue >= thr) = 1;
end

end



function [plv, nPLV] = computePhaseMetric(data,channelNames,nPermutation)
% Description
%	[PTE,PLV, F] = computePhaseMetric(data,channelNames)


%chMask	 = ~cellfun(@isempty,regexp(channelNames,'E[0-9]*'));
nFilters = numel(data);
nTrials  = numel(data(1).trial);
plv 	 = complex(zeros(nFilters,nTrials),zeros(nFilters,nTrials));

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
        
        try
            [tmp,tmpSurr] = phaseAmplitudeCorrelation(data(fIdx).trial{trialIdx},nPermutation);
        catch
            nTrials
            trialIdx
        end
        plv(fIdx,trialIdx) = tmp;
        plvSurr(fIdx,trialIdx,:) = tmpSurr;
        
    end
    
end
% get the mean across trials for each filter
plv 	= squeeze(mean(plv,2));

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



function data = preproc(data,iChannels,chFlag,trialLength)
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
wb = 0.25;
% stop band
ws = 2;
% multiplier
m = sqrt(2);

f = zeros(nFilters,1);
dataOrig = data;
for fIdx = 1:nFilters
    
    passBandLp = fn + (wb * fn)/2;
    passBandHp = fn - (wb * fn)/2;
    stopBandLp = fn * ws;
    stopBandHp = fn / ws;
    f(fIdx) = fn;
    fn = fn * m;
    
    cfg.lpfreq = passBandLp;
    cfg.hpfreq = passBandHp;
    cfg.lpfilter = 'yes';
    cfg.hpfilter = 'yes';
    cfg.lpfiltertype = 'but';
    cfg.hpfiltertype = 'but';
    cfg.lpfiltdir = 'twopass';
    cfg.hpfiltdir = 'twopass';
    cfg.padding  = trialLength;
    cfg.continuous = 'yes';
    data = ft_preprocessing(cfg,dataOrig);
    dataFiltered(fIdx) = preproc(data,iChannels,chFlag,trialLength);
    
end

end
