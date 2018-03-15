function varargout = process_walking_StrideAnalyses( varargin )
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
sProcess.Comment     = 'Stride Analyses';
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
sProcess.options.doWarping.Comment   = 'time warping ';
sProcess.options.doWarping.Type      = 'checkbox';
sProcess.options.doWarping.Value     = true;
sProcess.options.saveOutput.Comment  = 'Save output to brainstormDB';
sProcess.options.saveOutput.Type     = 'checkbox';
sProcess.options.saveOutput.Value    = false;
sProcess.options.normOnStride.Comment= 'Normalize on Stride';
sProcess.options.normOnStride.Type   = 'checkbox';
sProcess.options.normOnStride.Value  = false';


end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
%% TODO remove %%
close all;

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

% we need to find files that have to be processed and group
% them for subjects and conditions
% First group all sInput.Comment together
conditionStrings 	= {sInputs.Condition};


walkingConMask      = ~cellfun(@isempty,regexp(conditionStrings,'(w|W)alking'));

offConMask          = ~cellfun(@isempty,regexpi(conditionStrings,'off'));
notMontageConMask 	= cellfun(@isempty,regexpi(conditionStrings,'visite'));

fileIndices = find( walkingConMask & offConMask & notMontageConMask );

f1 = figure('papertype','a4','paperorientation','portrait','Visible','on');

for fileIdx = fileIndices
    
    walkingStruct = in_bst_timefreq(sInputs(fileIdx).FileName);
    
    parentStruct 	= bst_process('GetInputStruct',walkingStruct.DataFile);
    parentData 		= in_bst(parentStruct.FileName);
    channelData		= in_bst_channel(sInputs(fileIdx).ChannelFile);
    iChannels		= channel_find(channelData.Channel,'SEEG');
    signals			= parentData.F(iChannels,:);
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
    peakVelocIdx    = find(~cellfun(@isempty,regexp(evNames,'peak')));
    peakVelocMask   = peakVelocIdx - 2 > 0 & peakVelocIdx + 5 <= numel(evNames);
    
    % count how many strides we have recorded
    nStrideLeft     = sum(strcmp(evNames,'peakVeloc_L'))-1;
    nStrideRight    = sum(strcmp(evNames,'peakVeloc_R'))-1;
    
    nStrides 		= nStrideLeft + nStrideRight;
    
    % prepare event mask to reject artefactual or incomplete step
    trialString = regexp(sInputs(fileIdx).Condition,'trial\d+','match');
    subjMask 	= ismember(lower(patNames),lower(sInputs(fileIdx).SubjectName));
    trialMask 	= (ismember(lower(trialStrings),lower(trialString)));
    stepRej  	= stepIds(and(subjMask,trialMask));
    
    % stride mask each stride is composed by two steps
    strideIndexes = [1:nStrides;2:nStrides+1]';
    strideRej 	  = find(sum(ismember(strideIndexes,stepRej),2));
    
    % we have to correct the event adding the offset
    % since they are referred to the 0 of the raw data
    evOffset 		= round(walkingStruct.Time(1)*fs);
    strideStart 	= strideStart - evOffset;
    
    % we create the normalized stride time vector
%     referenceTimeVector = -1:1/fs:(1.5-1/fs);
    
    doubleSuppDur 		= floor(0.19*fs);
    singleSuppDur 		= floor(0.4*fs);
    acPhaseDur			= floor((.4*2/3)*fs);
    decPhaseDur			= floor((.4/3)*fs);
    
    % create the reference vector for morphing
    referenceStance 	= 400 + [ -doubleSuppDur-acPhaseDur...
        -acPhaseDur 0 +decPhaseDur ...
        +doubleSuppDur+decPhaseDur ...
        +doubleSuppDur+decPhaseDur+acPhaseDur ...
        +doubleSuppDur+decPhaseDur+singleSuppDur ...
        +doubleSuppDur+decPhaseDur+singleSuppDur+doubleSuppDur];
        
    referenceVector		= [1 referenceStance 1000];
    plotIdx = 1;
    
    % strideCheck = evNames(bsxfun(@plus,(-2:4),(3:2:numel(evNames)-2)'));
    % strideCheck = evNames(bsxfun(@plus,(-2:5),(peakVelocIdx(peakVelocMask))'));
    % this is the label of the central event
    % nStrides = size(strideCheck,1);
    
    % check that data are ordered correctly for each stride
    matchingString = {'heelcontact_[L|R]', 'toeoff_[R|L]',...
        'peakVeloc_[R|L]','heelcontact_[R|L]',...
        'toeoff_[L|R]','peakVeloc_[L|R]',...
        'heelcontact_[L|R]','toeoff_[R|L]'};
        
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
        
        % we take into account 2.5 seconds long time window centered at
        % peakVelocity
        freqMask 	 = walkingStruct.Freqs > 6;
        dataTF 		 = abs(walkingStruct.TF(:,timeWindow,freqMask)).^2;
        
        if (sProcess.options.normOnStride.Value)
            % normalize on the mean of the entire recording
            normFactor = mean(walkingStruct.TF(:,:,freqMask),2);
            normDataTF = bsxfun(@rdivide,dataTF,bsxfun(@minus,dataTF,normFactor),normFactor);
            
        else
            % normalize on swing only
            % this means that we are not normalizing here
            % but do it later after warping
            % for simplicity in the code. NOT SURE IF THIS IS THE CORRECT WAY
            % TODO check whether we have an effect of normalization order
            
        end
        
        
        %	fprintf('[%d]',strideIdx);
        %	fprintf('%s ',evNames{(-2:4) + strideIdx});
        %	fprintf('\n');
        strideRaw = signals(:,timeWindow)';
        f 		  = walkingStruct.Freqs(freqMask);
        
        % then create the time-warping vector in samples
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
            finalRaw(1,:) = mixingMatrix * strideRaw(:,1);
            finalRaw(2,:) = mixingMatrix * strideRaw(:,2);
            
            tAxis	= [-499:500]./fs;
%             tEvAxis = repmat(referenceTimeVector(referenceStance(2:end-1)),2, 1);
            tEvAxis = repmat((referenceStance(1:end-1)-400)./400,2,1);
            
        else
            %	finalTF 	= dataTF.*1e12;
            %   finalRaw	= strideRaw';
            %	zLimit 		= [min(finalTF(:)) max(finalTF(:))];
            %	tAxis		= (-399:400)/fs;
            %	tEvAxis 	= repmat(originalVector([3 4 5])./fs,2, 1);
        end
        
        footLabel = regexp(evNames(strideIdx),'[L|R]','match');
        footLabel = footLabel{:}{:};
        % this is the label of the central event
        
        fprintf('[%d]',strideIdx);
        fprintf('%s \n',evNames{strideIdx}); 
        
        if strcmp(footLabel,'L')
            % left foot swing => central hc_L
            %  => stnContra is rightSTN == idx 1
            controLatIdx = 1;
            if strcmp(mostAffSide,'L')
                % if STN- is L => this swing is relative to
                % STN+ => plot on right side of page
                stnIdx = 2;
            else
                % if STN- is R => this swing is relative
                % to STN- => plot on left side of page
                stnIdx = 1;
            end
        else
            % right foot swing
            controLatIdx = 2;
            if strcmp(mostAffSide,'L')
                % if STN- is R => this swing is relative
                % to STN- => plot on left side of page
                stnIdx = 1;
            else
                % if STN- is L => this swing is relative to
                % STN+ => plot on right side of page
                stnIdx = 2;
            end
        end
        
        stridePhaseDur{subjectIdx,stnIdx} = ...
            cat(1,stridePhaseDur{subjectIdx,stnIdx},...
            diff(originalVector));
        
        % we save for each subject the time-warped ERSD normalized over***
        stnRawResults{subjectIdx,stnIdx} = ...
            cat(1,stnRawResults{subjectIdx,stnIdx},...
            finalTF(controLatIdx,:,:));
        
        stnData{subjectIdx,stnIdx} = ...
            cat(1,stnData{subjectIdx,stnIdx},...
            finalRaw(controLatIdx,:));
        
        plotIdx = plotIdx + 1;
        clear finalTF;
        
    end % for stride
    
    clear finalTF;
    
end % for sInputs files

f2 = figure('papertype','a4','paperorientation','portrait','Visible','on');
colormap(jet(256));

f3 = figure('papertype','a4','paperorientation','portrait','Visible','on');
colormap(jet(256));
%

lowBetaMask = f >= 8 & f <= 16;
highBetaMask = f > 16 & f <= 30;
gammaMask  = f > 30 & f < 80;

% patientsOrder = {'wue03','wue09','wue04','wue02','wue10','wue07','wue06','wue11'};
patientsOrder = {'wue09','wue04','wue02','wue07','wue06'};
[~,ord] = ismember(patientsOrder,subjectNameOrdered);

plotIdx = 1;
% method 3 == ERSD - EEGLAB gain model (Grandchamp Delorme Front. 
method = 3;

tEvAxis = tEvAxis(:,1:4);

% groupStnMostAffERSD
% groupStnLessAffERSD

% loop across subjects
for ii = ord
    % for a given subject get data in the form
    % trial x time x freq. 
    stnMostAff = stnRawResults{ii,1};
    stnLessAff = stnRawResults{ii,2};
    
    stnMostAffERSD = computeERSD(stnMostAff,referenceStance,method);
    stnLessAffERSD = computeERSD(stnLessAff,referenceStance,method);
    
    groupStnMostAffERSD(ii,:,:,:) = stnMostAffERSD;
    groupStnLessAffERSD(ii,:,:,:) = stnLessAffERSD;
%     [pvalue(1,:,:), unCorrPvalue(1,:,:)] = ...
%         runPermutationTest(stnMostAffERSD,...
%         stnMostAff,100,referenceStance,method);
%     
%     [pvalue(2,:,:), unCorrPvalue(2,:,:)] = ...
%         runPermutationTest(stnLessAffERSD,...
%         stnLessAff,100,referenceStance,method);
      

    [pvalue(1,:,:), unCorrPvalue(1,:,:)] = ...
        runStatisticalValidationERSP(stnMostAffERSD,...
        stnMostAff,100,referenceStance,method);
    
    [pvalue(2,:,:), unCorrPvalue(2,:,:)] = ...
        runStatisticalValidationERSP(stnLessAffERSD,...
        stnLessAff,100,referenceStance,method);
    
    statSignificance = zeros(size(pvalue));
    statSignificance(unCorrPvalue < 0.05) = 1;
    
    ax1=subplot(nSubjects,2,2*(plotIdx-1)+1,'NextPlot','add');
    h=imagesc(tAxis,f(lowBetaMask | highBetaMask | gammaMask), stnMostAffERSD(:,(lowBetaMask | highBetaMask | gammaMask))',[-3 3]);
    h.AlphaData = squeeze(statSignificance(1,:,(lowBetaMask | highBetaMask | gammaMask)))';
    plot(tEvAxis,repmat([min(f);60],[1 numel(tEvAxis(1,:))]),'k--');
    axis xy;
    
    set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc'});
    xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
    ylim([min(f(lowBetaMask | highBetaMask | gammaMask)), max(f(lowBetaMask | highBetaMask | gammaMask))]);
 
    set(ax1,'Parent',f2);  
    ax2=subplot(nSubjects,2,2*(plotIdx-1)+2,'NextPlot','add');
    h=imagesc(tAxis,f(lowBetaMask | highBetaMask| gammaMask),stnLessAffERSD(:,(lowBetaMask | highBetaMask | gammaMask))',[-3 3]);
    h.AlphaData = squeeze(statSignificance(2,:,(lowBetaMask | highBetaMask| gammaMask)))';
    plot(tEvAxis,repmat([min(f);60],[1 numel(tEvAxis(1,:))]),'k--');
    xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
    axis xy;
    set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc'});
    
    ylim([min(f(lowBetaMask | highBetaMask| gammaMask)), max(f(lowBetaMask | highBetaMask| gammaMask))]);
    annotation('textbox',[0.05, 0.85-(plotIdx-1)*0.1, 0.1, 0.05],...
        'String',subjectNameOrdered{ii},'LineStyle','None');
    set(gca,'Parent',f2);
   
    ax3=subplot(nSubjects,2,2*(plotIdx-1)+1,'NextPlot','add');
    plot(tAxis,mean(stnMostAffERSD(:,highBetaMask),2),'r');
    plot(tAxis,mean(stnMostAffERSD(:,lowBetaMask),2),'g');
    plot(tAxis,mean(stnMostAffERSD(:,gammaMask),2),'b');
    plot(tEvAxis,repmat([-2;2],[1 numel(tEvAxis(1,:))]),'k--');

    xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
    set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc'});
    ylim([-1 1]);
    set(ax3,'Parent',f3);
    
    ax4=subplot(nSubjects,2,2*(plotIdx-1)+2,'NextPlot','add');
    plot(tAxis,mean(stnLessAffERSD(:,highBetaMask),2),'r');
    plot(tAxis,mean(stnLessAffERSD(:,lowBetaMask),2),'g');
    plot(tAxis,mean(stnLessAffERSD(:,gammaMask),2),'b');
    plot(tEvAxis,repmat([-2;2],[1 numel(tEvAxis(1,:))]),'k--')

    xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
    ylim([-1 1]);
    set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc'});
    annotation('textbox',[0.05, 0.85-(plotIdx-1)*0.1, 0.1, 0.05],...
        'String',subjectNameOrdered{ii},'LineStyle','None');
    set(ax4,'Parent',f3);

    plotIdx = plotIdx + 1;  
    
end

annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');

%% group level test
grpStnMostAff = squeeze(mean(groupStnMostAffERSD));
grpStnLessAff = squeeze(mean(groupStnLessAffERSD));
figure,
colormap(jet(256));
subplot(1,2,1)
imagesc(tAxis,f,grpStnMostAff',[-2 2]);
xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
ylim([min(f) 30]);
axis xy;
set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc'});

subplot(1,2,2)
imagesc(tAxis,f,grpStnLessAff',[-2 2]);
xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
ylim([min(f) 30]);
axis xy;
set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc'});

figure,
grpStnMostHighBeta = squeeze(mean(groupStnMostAffERSD(:,:,highBetaMask),3));
grpStnLessHighBeta = squeeze(mean(groupStnLessAffERSD(:,:,highBetaMask),3));
grpStnMostLowBeta  = squeeze(mean(groupStnMostAffERSD(:,:,lowBetaMask),3));
grpStnLessLowBeta  = squeeze(mean(groupStnLessAffERSD(:,:,lowBetaMask),3));
grpStnMostGamma    = squeeze(mean(groupStnMostAffERSD(:,:,gammaMask),3));
grpStnLessGamma    = squeeze(mean(groupStnLessAffERSD(:,:,gammaMask),3));

[grpStnMostHighBetaCL,grpStnMostHighBetaMean] = myBootstrap(grpStnMostHighBeta,nSubjects,10);
[grpStnLessHighBetaCL,grpStnLessHighBetaMean] = myBootstrap(grpStnLessHighBeta,nSubjects,10);
[grpStnMostLowBetaCL, grpStnMostLowBetaMean]  = myBootstrap(grpStnMostLowBeta,nSubjects,10);
[grpStnLessLowBetaCL, grpStnLessLowBetaMean]  = myBootstrap(grpStnLessLowBeta,nSubjects,10);
[grpStnMostGammaCL, grpStnMostGammaMean]      = myBootstrap(grpStnMostGamma,nSubjects,10);
[grpStnLessGammaCL, grpStnLessGammaMean]      = myBootstrap(grpStnLessGamma,nSubjects,10);


%% compare STN- and STN+ for LowBeta
p = mySignrank(grpStnMostLowBeta-grpStnLessLowBeta);

subplot(3,1,1,'NextPlot','add')
plot(tAxis,grpStnMostLowBetaMean,'Color',[0 109 219]./255);
plot(tAxis,grpStnLessLowBetaMean,'Color',[255 109 182]./255);
fill_between(tAxis,grpStnMostLowBetaCL(1,:),grpStnMostLowBetaCL(2,:),tAxis,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(tAxis,grpStnLessLowBetaCL(1,:),grpStnLessLowBetaCL(2,:),tAxis,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');

plot(tAxis(p<(0.05/numel(p))),grpStnMostLowBetaMean(p<(0.05/numel(p))),'s','MarkerFaceColor',[0 109 219]./255,'MarkerEdgeColor','none');
plot(tAxis(p<(0.05/numel(p))),grpStnLessLowBetaMean(p<(0.05/numel(p))),'s','MarkerFaceColor',[255 109 182]./255,'MarkerEdgeColor','none');
plot(tEvAxis,repmat([-2;2],[1 numel(tEvAxis(1,:))]),'k--')
set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc','To'});
xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
% ylim([-.5 .5]);

%% compare STN- and STN+ for highBeta
p = mySignrank(grpStnMostHighBeta-grpStnLessHighBeta);

subplot(3,1,2,'NextPlot','add')
plot(tAxis,grpStnMostHighBetaMean,'Color',[0 109 219]./255)
plot(tAxis,grpStnLessHighBetaMean,'Color',[255 109 182]./255);
legend({'STN-','STN+'});
plot(tAxis(p<(0.05/numel(p))),grpStnMostHighBetaMean(p<(0.05/numel(p))),'s','MarkerFaceColor',[0 109 219]./255,'MarkerEdgeColor','none')
plot(tAxis(p<(0.05/numel(p))),grpStnLessHighBetaMean(p<(0.05/numel(p))),'s','MarkerFaceColor',[255 109 182]./255,'MarkerEdgeColor','none');
fill_between(tAxis,grpStnMostHighBetaCL(1,:),grpStnMostHighBetaCL(2,:),tAxis,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(tAxis,grpStnLessHighBetaCL(1,:),grpStnLessHighBetaCL(2,:),tAxis,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');
set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc','To'});
plot(tEvAxis,repmat([-2;2],[1 numel(tEvAxis(1,:))]),'k--')
xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
% ylim([-.5 .5]);

%% compare STN- and STN+ for gamma
p = mySignrank(grpStnMostGamma-grpStnLessGamma);

subplot(3,1,3,'NextPlot','add')
plot(tAxis,grpStnMostGammaMean,'Color',[0 109 219]./255)
plot(tAxis,grpStnLessGammaMean,'Color',[255 109 182]./255);
fill_between(tAxis,grpStnMostGammaCL(1,:),grpStnMostGammaCL(2,:),tAxis,'FaceColor',[0 109 219]./255,'FaceAlpha',0.2,'EdgeColor','None');
fill_between(tAxis,grpStnLessGammaCL(1,:),grpStnLessGammaCL(2,:),tAxis,'FaceColor',[255 109 182]./255,'FaceAlpha',0.2,'EdgeColor','None');
plot(tEvAxis,repmat([-2;2],[1 numel(tEvAxis(1,:))]),'k--')
plot(tAxis(p<(0.05/numel(p))),grpStnMostGammaMean(p<(0.05/numel(p))),'s','MarkerFaceColor',[0 109 219]./255,'MarkerEdgeColor','none')
plot(tAxis(p<(0.05/numel(p))),grpStnLessGammaMean(p<(0.05/numel(p))),'s','MarkerFaceColor',[255 109 182]./255,'MarkerEdgeColor','none');
set(gca,'XTick',tEvAxis(1,:),'XTickLabel',{'Hc','To','Vp','Hc','To','Vp','Hc','To'});
xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
% ylim([-.5 .5]);
% fname = fullfile(getenv('HOME'),'Dropbox','Isaias_group','walking','figs',...
%     'avgZScoreStrideMod.png');

% print(f2,'-dpng',fname);
end % function

function [p] = mySignrank(data)
    % data subject x time
    [~,nTimes] = size(data);
    
    p = ones(1,nTimes);
    for ii=1:nTimes
        p(ii) = signrank(data(:,ii));
    end

end

function [confLimit,dataMean] = myBootstrap(data,nSubject,nBootstraps)

    % data matrix has nSubjects x t 
    % bootstrap the C.L. for mean
    bootstrapIndexes = randi(nSubject,nBootstraps,nSubject);
    [~,nTimes] = size(data);
    
    % compute mean across subjects
    dataMean = mean(data);
        
    % currBootstraps will contain nBootstraps x nStn x f
    currBootstraps = nan(nBootstraps,nTimes);
    for idx = 1:nBootstraps
        currBootstraps(idx,:) = squeeze(mean(data(bootstrapIndexes(idx,:),:)));
    end
    
    % confLimit will contain [UB LB] x nStn x f
    % dataMean is 1 x 2 x f
    % thus we replicate on the first dim
%     A = repmat(dataMean,1,2,1);
    % currBootstrap after mean should be 1 x 2 x f
    % thus we replicate on the first dim
    confLimit = prctile(currBootstraps,[2.5 97.5],1);

end

function [pvalue, unCorrpvalue] = runPermutationTest(obsERSD,stnData,nPermutation,referenceStance,method)
%RUNPERMUTATIONTEST Description
%	PVALUE = RUNPERMUTATIONTEST(STANCE,SWING,NPERMUTATION) Long description
%
pvalue    = zeros(size(obsERSD));
nSwing    = size(stnData,1);

dataPerm  = stnData;

% we perform a permutation test for each STN separatelly
for permIdx = 1:nPermutation
    
    % for each swing we randomly split the signal in two chunks
    % and rotate them
    for swingIdx = 1:nSwing
        
        dataPerm(swingIdx,:,:) = randCircShift(stnData(swingIdx,:,:));
        
    end
    
    % compute permutated statistics
    permERSD = computeERSD(dataPerm,referenceStance,method);
    
    % compute pvalues for all frequencies and all time points.
    pvalue = pvalue + double(permERSD > obsERSD)./nPermutation;
    
end
unCorrpvalue = pvalue;
pvalue = fdrCorrection(pvalue,0.05);

end

function [pvalue, unCorrpvalue] = runStatisticalValidationERSP(obsERSD,stnData,nPermutation,referenceStance,method)

pvalue = zeros(size(obsERSD));

[nTrials, ~, ~] = size(stnData);

dataPerm  = stnData;

baselineIndexes = referenceStance(1):referenceStance(2);
nTimeBaseline = numel(baselineIndexes);

% we perform a permutation test for each STN separatelly
for permIdx = 1:nPermutation
    
    % extract baseline samples within trial
    baselineSamples = stnData(:,baselineIndexes,:);
    
    % permute the values in time and trials for each frequency
    baselineSamples = baselineSamples(randperm(nTrials),randperm(nTimeBaseline),:);
    
    % save permuted baseline samples into Permuted Data array
    dataPerm(:,nTimeBaseline,:) = baselineSamples;
    
    % compute permutated statistics
    permERSD = computeERSD(dataPerm,referenceStance,method);
    
    % compute pvalues for all frequencies and all time points.
    pvalue = pvalue + double(permERSD > obsERSD)./nPermutation;
    
end

unCorrpvalue = pvalue;
pvalue = fdrCorrection(pvalue,0.05);

end

function A = randCircShift(A)

idx 	= randi(size(A,2),1);
A(1,:,:)= cat(2,A(1,idx:end,:),A(1,1:idx-1,:));

end

function pvalue = fdrCorrection(pvalue, alpha)
%FDRCORRECTION Description
%	PVALUE  = FDRCORRECTION() Long description
%

tmpPvalue 	= sort(pvalue(:));
N 			= numel(pvalue);
FDR 		= alpha.*(1:N)./N;
thr 		= FDR(find(tmpPvalue <= FDR',1,'last'));
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

% compute the normalization factor concatenating all baseline
% and computing the mean across trials
[n,Time,f] = size(stnData);
tBaseline = referenceStance(1):referenceStance(end);
t = numel(tBaseline);

numFactor = mean(stnData(:,tBaseline,:),2);
denFactor = std(reshape(stnData(:,tBaseline,:),[n*t,f]));

% baseline correction
switch(method)
    case 1
        % rel change
        stnResult = bsxfun(@rdivide,bsxfun(@minus,stnData,numFactor),numFactor);
    case 2
        % pseudo-zscore
        stnResult = bsxfun(@rdivide,bsxfun(@minus,stnData,numFactor),denFactor);
    case 3
        % full-trial baseline correction
        stnResult  = stnData ./ repmat(numFactor,[1 Time 1]);
%         stnResult  = real(10.*log10(stnResult));
end

% mean across trials
stnResult = squeeze(mean(stnResult));

% recompute ERSD
baseline   = referenceStance(1):referenceStance(2);
normFactor = mean(stnResult(baseline,:));
stnResult  = bsxfun(@rdivide, stnResult,normFactor);
stnResult  = real(10.*log10(stnResult));

end
