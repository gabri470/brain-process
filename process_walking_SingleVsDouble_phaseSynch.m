function varargout = process_walking_SingleVsDouble_phaseSynch(varargin )
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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Single vs Double ';
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

singleSuppTcourse = cell(nSubjects,2);
doubleSuppTcourse = cell(nSubjects,2);
singleSuppPlv = cell(nSubjects,1);
doubleSuppPlv = cell(nSubjects,1);

% we need to find files that have to be processed and group
% them for subjects and conditions
% First group all sInput.Comment together
conditionStrings 		= {sInputs.Condition};


walkingConMask      = ~cellfun(@isempty,regexp(conditionStrings,'(w|W)alking'));

offConMask          = ~cellfun(@isempty,regexpi(conditionStrings,'off'));
notMontageConMask 	= cellfun(@isempty,regexpi(conditionStrings,'visite'));

fileIndices = find( walkingConMask & offConMask & notMontageConMask );

nFilters = 13;

% nominal frequency and central frequency
fn = 2;
% band-flat top and band pass width
wb = 0.5;
% stop band
% ws = 2;
% multiplier
m = sqrt(2);

f = zeros(nFilters,1);

for fIdx = 1:nFilters

    f(fIdx) = fn;
    fn = fn * m;
end

options.MorletFc = 1;
options.MorletFwhmTc = 2;
options.Freqs = f;

for fileIdx = fileIndices
    
    
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
    subjMask 	= ismember(lower(patNames),lower(sInputs(fileIdx).SubjectName));
    trialMask = (ismember(lower(trialStrings),lower(trialString)));
    stepRej  	= stepIds(and(subjMask,trialMask));
    
    % stride mask each stride is composed by two steps
    strideIndexes = [1:nStrides;2:nStrides+1]';
    strideRej 	  = find(sum(ismember(strideIndexes,stepRej),2));
    
    % we have to correct the event adding the offset
    % since they are referred to the 0 of the raw data
    evOffset      = round(walkingStruct.Time(1)*fs);
    strideStart 	= strideStart - evOffset;
    
    % we create the normalized stride time vector
    referenceTimeVector = -1:1/fs:(1-1/fs);
    doubleSuppDur 		= floor(0.19*fs);
    singleSuppDur 		= floor(0.4*fs);
    
    referenceStance 	= 400 +[ -doubleSuppDur-singleSuppDur -singleSuppDur 0 ...
        +doubleSuppDur +doubleSuppDur+singleSuppDur ];
    
    referenceVector		= [1 referenceStance 800];
    
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
        
             
        footLabel = regexp(evNames(strideIdx),'[L|R]','match');
        footLabel = footLabel{:}{:}; % this is the label of the central hc

        
        
        % we need to switch foot Label since a left foot stride
        % will have a hcR event as central (t0) event
        % stnIdx are ordered as STN- first and STN+ next
        if strcmp(footLabel,'L')
            % right foot stride
            controLatIdx = 2;
            if strcmp(mostAffSide,'L')
                stnIdx	= 1;
            else
                stnIdx	= 2;
            end
        else
            
            % left foot stride
            controLatIdx = 1;
            if strcmp(mostAffSide,'L')
                stnIdx	= 2;
            else
                stnIdx	= 1;
            end
        end
        
        %								[ idx ]
        % (hc_L) *tof_R  hc_R   *tof_L (hc_L)
        %    ^ start-2		|t0				      ^ start + 2 == end stride
        timeWindow = strideStart(strideIdx) + (-399:400);
        freqMask   = f > 0;
        
        dataTF = morlet_transform(signals, parentData.Time, options.Freqs,...
												options.MorletFc, options.MorletFwhmTc, 'n');  
        dataTF = dataTF(:,timeWindow,freqMask);
        
        %dataTF 		 = walkingStruct.TF(:,timeWindow,freqMask);
        
        strideRaw	 = signals(:,timeWindow)'.*1e6;
        %f 			 = walkingStruct.Freqs(freqMask);
        
        % then create the time-warping vector
     
        originalVector = [1 (strideStart((-2:1:2) + strideIdx) ...
            - timeWindow(1))  800];

        
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
           
            
            % tf_* -> hc_* => single support controlateral singleSupp
            doubleSuppWnd = referenceStance(2):referenceStance(3);
            % hc_* -> tf_* => double support subsequent to above singleSupp phase
            singleSuppWnd	= referenceStance(3):referenceStance(4);
            tAxis			= referenceTimeVector;
            
        else
            finalTF 	= [];
            finalRaw	= strideRaw';
            zLimit 		= [min(finalTF(:)) max(finalTF(:))];
            % tf_* -> hc_*
            doubleSuppWnd	= originalVector(2):originalVector(3);
            % tf_* -> hc_*
            singleSuppWnd	= originalVector(3):originalVector(4);
            
            tAxis			= (-399:400)/fs;
        end
        % here depending on the options chosen we might
        % end up having a time-warped data with fixed stride lengths
        % or individual stride lenghts with their original data
        % thus we need to account for different sizes
        %
        
        %stnRawData{subjectIdx,stnIdx} = cat(1,stnRawData{subjectIdx,stnIdx},finalTF);
        
        if sProcess.options.normalizeOnStride.Value
            normFactor 	= repmat(mean(finalTF,2),[1 numel(timeWindow) 1]);
            A 					= finalTF./normFactor;
            singleSuppData 	= A(controLatIdx,singleSuppWnd,:);
            doubleSuppData	= A(controLatIdx,doubleSuppWnd,:);
            
        else
            singleSuppData = finalTF(controLatIdx,singleSuppWnd,:);
            doubleSuppData = finalTF(controLatIdx,doubleSuppWnd,:);
            
        end
        
        % ch x time x freq
        singleSuppPhaseData  = angle(finalTF(:,singleSuppWnd,:));
        doubleSuppPhaseData = angle(finalTF(:,doubleSuppWnd,:));
        
        singleSuppPlvSingleTrial = abs(mean(exp(1i.*squeeze(diff(singleSuppPhaseData,[],1)))));
        doubleSuppPlvSingleTrial = abs(mean(exp(1i.*squeeze(diff(doubleSuppPhaseData,[],1)))));
        
        singleSuppPlv{subjectIdx}  = cat(1,singleSuppPlv{subjectIdx},singleSuppPlvSingleTrial);
        doubleSuppPlv{subjectIdx}  = cat(1,doubleSuppPlv{subjectIdx},doubleSuppPlvSingleTrial);
        
        
        if strcmp(sProcess.options.method.Value,'zscore')
            % z score singleSupp
            singleSuppZScored = bsxfun(@rdivide,bsxfun(@minus,finalTF(:,singleSuppWnd,:),...
                mean(finalTF(:,singleSuppWnd,:),2)),...
                std(finalTF(:,singleSuppWnd,:),[],2));
            
            % z score doubleSupp
            doubleSuppZScored = bsxfun(@rdivide,bsxfun(@minus,finalTF(:,doubleSuppWnd,:),...
                mean(finalTF(:,doubleSuppWnd,:),2)),...
                std(finalTF(:,doubleSuppWnd,:),[],2));
            
            finalTFZScored = bsxfun(@rdivide,...
                bsxfun(@minus,finalTF(:,doubleSuppWnd,:),...
                mean(finalTF(:,singleSuppWnd,:),2)),...
                std(finalTF(:,singleSuppWnd,:),[],2));
            
        elseif(strcmp(sProcess.options.method.Value,'rchange'))
            % relative change
            finalTFZScored = bsxfun(@rdivide,...
                bsxfun(@minus,finalTF(:,doubleSuppWnd,:),...
                mean(finalTF(:,singleSuppWnd,:),2)),...
                mean(finalTF(:,singleSuppWnd,:),2));
            
        elseif(strcmp(sProcess.options.method.Value,'difference'))
            finalTFZScored = doubleSuppZScored - singleSuppZScored;
           
        end
        
        if sProcess.options.doWarping.Value
            stnResults{subjectIdx,stnIdx} = ...
                cat(1,stnResults{subjectIdx,stnIdx},...
                finalTFZScored(controLatIdx,:,:));
            
            % for each subject and stn-(1)/+(2) we save the time-averaged
            % singleSupp and doubleSupp power
            singleSuppTcourse{subjectIdx,stnIdx} = ...
                cat(1,singleSuppTcourse{subjectIdx,stnIdx},singleSuppData);
            doubleSuppTcourse{subjectIdx,stnIdx} = ...
                cat(1,doubleSuppTcourse{subjectIdx,stnIdx},doubleSuppData);
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

% stnResults holds the zscored-power modulation (1,t,f) for doubleSupp {:,:,1}
% and singleSupp {:,:,2} for each subjects {ii,:,:} and each strides and for
% both STN- {:,1,:} and STN+ {:,2,:}.
% We first average each element of stnResults across time in order to have
% a frequency modulation
stnMeans	= cellfun(@mean,stnResults,'UniformOutput',false);

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

for fIdx = 1:11
    
    passBands(fIdx,:) = [fn-(wb*fn)/2  fn+(wb * fn)/2];
    fBands(fIdx) = fn;
    fn = fn * m;
end

[~, hpMask] = meshgrid(f,passBands(:,1));
[fMask, lpMask] = meshgrid(f,passBands(:,2));
passBandsMask = (fMask >= hpMask & fMask < lpMask);

f2 = figure('papertype','a4','paperorientation','portrait','Visible','on');

pvalue = zeros(numel(f),2);
corrPvalue = zeros(numel(f),2);

patientsOrder = {'wue03','wue09','wue04','wue02','wue10','wue07','wue06','wue11'};
[~,ord] = ismember(patientsOrder,subjectNameOrdered);

plotIdx = 1;

singleSuppPlv = cellfun(@mean,singleSuppPlv,'UniformOutput',false);
doubleSuppPlv = cellfun(@mean,doubleSuppPlv,'UniformOutput',false);

singlePlv = cat(1,singleSuppPlv{:});
doublePlv = cat(1,doubleSuppPlv{:});

% swngPlv = avgInfBands(swngPlv,passBandsMask);
% stncPlv = avgInfBands(stncPlv,passBandsMask);
fBands = f;

[doubleSuppCL, doubleSuppMean] = myBootstrap(doublePlv,nSubjects,20);
[singleSuppCL, singleSuppMean] = myBootstrap(singlePlv,nSubjects,20);

p = mySignrank(doublePlv-singlePlv);

h1 = plot(fBands,doubleSuppMean,'Color',[0 109 219]./255);
hold on
h2 = plot(fBands,singleSuppMean,'Color',[255 109 182]./255);
legend([h1 h2], {'doubleSupp','singleSupp'});

fill_between(fBands,doubleSuppCL(1,:),doubleSuppCL(2,:),fBands,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(fBands,singleSuppCL(1,:),singleSuppCL(2,:),fBands,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');

xlim([min(fBands) 60]);

%for ii = ord
%     
%     fprintf('Statistcs on %s\n',subjectNameOrdered{ii});
%     
%     %% these below contain complex wavelet coeff
%     phaseStanceStnMost = angle(doubleSuppTcourse{ii,1});
%     phaseStanceStnLess = angle(doubleSuppTcourse{ii,2});
%     
%     phaseSwingStnMost  = angle(singleSuppTcourse{ii,1});
%     phaseSwingStnLess  = angle(singleSuppTcourse{ii,2});
%     
%     %% phaseStanceStnMost and other are trial x time x freq matrices 
%     plvStance = squeeze(mean(abs(sum(exp(phaseStanceStnMost-phaseStanceStnLess),2))));
%     plvSwing  = squeeze(mean(abs(sum(exp(phaseSwingStnMost-phaseSwingStnLess),2)))); 
%     
%          
%     subplot(nSubjects,2,(plotIdx-1),'NextPlot','add')
%     plot(f,plvStance,f,plvSwing)
%    
%    
%    
%     [corrPvalue(:,1), pvalue(:,1)] = runPermutationTest(stnMeans{ii,1},...
%         doubleSuppTcourse{ii,1},singleSuppTcourse{ii,1},100);
%     [corrPvalue(:,2), pvalue(:,2)] = runPermutationTest(stnMeans{ii,2},...
%         doubleSuppTcourse{ii,2},singleSuppTcourse{ii,2},100);
%     
%     
%     % need multiple comparison correction  FDR
%     cororrSignificanceMask = nan(size(pvalue));
%     unCorrSignificanceMask(pvalue < 0.05) = 4;
%     
%     
%     % this is the STN-
%     annotation('textbox',[0.05, 0.85-(plotIdx-1)*0.1, 0.1, 0.05],...
%         'String',subjectNameOrdered{ii},'LineStyle','None');
%     
%     subplot(nSubjects,2,2*(plotIdx-1)+1,'NextPlot','add')
%     plot(f,squeeze(mean(stnMeans{ii,1},2)))
%     plot([0 0;90 90],[-3 3;-3 3],'k--');
%     plot(f,unCorrSignificanceMask(:,1),'k.','MarkerSize',16);
%     plot(f,corrSignificanceMask(:,1),'r.','MarkerSize',16);
%     xlim([6 80]);
%     ylim([-5 5]);
%     
%     
%     
%     subplot(nSubjects,2,2*(plotIdx-1)+2,'NextPlot','add')
%     plot(f,squeeze(mean(stnMeans{ii,2},2)))
%     plot([0 0;90 90],[-3 3;-3 3],'k--');
%     plot(f,unCorrSignificanceMask(:,2),'k.','MarkerSize',16);
%     plot(f,corrSignificanceMask(:,2),'r.','MarkerSize',16);
%     xlim([6 60]);
%     ylim([-5 5]);rSignificanceMask = nan(size(corrPvalue));
%     corrSignificanceMask(corrPvalue < 0.05) = 4;
%     
%     unCorrSignificanceMask = nan(size(pvalue));
%     unCorrSignificanceMask(pvalue < 0.05) = 4;
%     
%     
%     % this is the STN-
%     annotation('textbox',[0.05, 0.85-(plotIdx-1)*0.1, 0.1, 0.05],...
%         'String',subjectNameOrdered{ii},'LineStyle','None');
%     
%     subplot(nSubjects,2,2*(plotIdx-1)+1,'NextPlot','add')
%     plot(f,squeeze(mean(stnMeans{ii,1},2)))
%     plot([0 0;90 90],[-3 3;-3 3],'k--');
%     plot(f,unCorrSignificanceMask(:,1),'k.','MarkerSize',16);
%     plot(f,corrSignificanceMask(:,1),'r.','MarkerSize',16);
%     xlim([6 80]);
%     ylim([-5 5]);
%     
%     
%     
%     subplot(nSubjects,2,2*(plotIdx-1)+2,'NextPlot','add')
%     plot(f,squeeze(mean(stnMeans{ii,2},2)))
%     plot([0 0;90 90],[-3 3;-3 3],'k--');
%     plot(f,unCorrSignificanceMask(:,2),'k.','MarkerSize',16);
%     plot(f,corrSignificanceMask(:,2),'r.','MarkerSize',16);
%     xlim([6 60]);
%     ylim([-5 5]);
%    plotIdx = plotIdx + 1;
%    
%end
% annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
% annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');
% 
% fname = fullfile('/','home','lgabri','Dropbox','Isaias_group','walking','figs',...
%     'avgZScoreStancevSwing_Corr.ps');
% 
% print(f2,'-dpsc2',fname);

end % function


function avg = avgInfBands(data,bandMasks)
    nSubjects = size(data,1);
    nFreq = size(bandMasks,1);
    avg = nan(nSubjects,nFreq);
    for fIdx = 1:nFreq
        avg(:,fIdx) = squeeze(mean(data(:,bandMasks(fIdx,:)),2));
    end
end

function [p] = mySignrank(data)
    % data subject x freq
    [~,nTimes] = size(data);
    
    p = ones(1,nTimes);
    for ii=1:nTimes
        %p(ii) = signrank(data(:,ii));
        p(ii) = ttest(data(:,ii));
    end

end

function [confLimit,dataMean] = myBootstrap(data,nSubject,nBootstraps)

    % data matrix has  nSubjects x f
    % bootstrap the C.L. for mean
    bootstrapIndexes = randi(nSubject,nBootstraps,nSubject);
    [~,nFreq] = size(data);
    
    %dataMean will be a 1 x f
    dataMean = mean(data);
    
    % currBootstraps will contain nBootstraps x f
    currBootstraps = nan(nBootstraps,nFreq);
    
    for idx = 1:nBootstraps
        currBootstraps(idx,:,:) = squeeze(mean(data(bootstrapIndexes(idx,:),:)));
    end
    
    % confLimit will contain [UB LB] x nStn x f
    confLimit = prctile(currBootstraps,[5 95]);
        
end

function [corrPvalue,unCorrPvalue] = runPermutationTest(dataObs,doubleSupp,singleSupp,nPermutation)
%RUNPERMUTATIONTEST
%	[CORRECTEDPVALUES,UNCORRECTEDPVALUES] = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION)
%
pvalue  = zeros(size(doubleSupp,3),1);
nDoubleSupp = size(doubleSupp,1);
nSingleSupp  = size(singleSupp,1);

dataObs = squeeze(mean(dataObs,2));

% we perform a permutation test for each STN separatelly
for permIdx = 1:nPermutation
    
    % concatenate each stride phase
    evGroup = [doubleSupp;singleSupp];
    
    % riffle indices
    riffledIdx = randperm(nDoubleSupp+nSingleSupp);
    
    % separate in two groups with same probablilty as
    % observed data (ie. sample size)
    dataStance = evGroup(riffledIdx(1:nDoubleSupp),:,:);
    dataSwing = evGroup(riffledIdx(nDoubleSupp+1:end),:,:);
    
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
f           = numel(pvalue);
tmpPvalue 	= sort(pvalue(:));
N 			= numel(pvalue);
FDR 		= alpha.*(1:N)./N;
thr 		= FDR(find(tmpPvalue <= FDR',1,'last'));
if isempty(thr)
    pvalue = ones(f,1);
else
    pvalue(pvalue >= thr) = 1;
end

end
