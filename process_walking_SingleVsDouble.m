function varargout = process_walking_singlevsDouble( varargin )
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

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Swing Modulation';
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

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

	nFiles = numel(sInputs);

	DATA_FOLDER='/media/lgabri/My Passport/';

  trialRejectionFile = fullfile(DATA_FOLDER,'trialRejection.csv');
  [patNames,trialStrings,stepIds] = textread(trialRejectionFile,...
																				'%s %*s %*s %*s%s%d%*s','delimiter',',');
	sideFile = fullfile(DATA_FOLDER,'patientSides.csv');
	[subjectNames, mostAffSides] = textread(sideFile,'%s %s\n','delimiter',',');
	nSubjects = numel(unique([{sInputs.SubjectName}]));
	stnResults = cell(nSubjects,2); 
	% the above var holds for each subjects the STN-/+ power for each stride
	
	currentSubject=[];
	subjectIdx = 0;
	subjectNameOrdered = cell(nSubjects,1);	OutputFiles = {};


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
			end


			% we need to read the angular velocity
			% and in order to do that we need to understand what trial and 
			% which sujbect we are analysing
			conditionString = regexp(sInputs(fileIdx).Condition,'trial\d+','match');
		  angVelocityFile = fullfile(DATA_FOLDER,currentSubject,...
																strcat(currentSubject,'_walking',...
																				conditionString,'.txt'));


			
			% compute sampling frequency
			fs = round(1/mean(diff( parentData.Time )));

			% filter cardiac events from gait-related events
			evGroupNames = {parentData.Events.label};
			gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'[heel|toe]'));

			% concat all heel contact events in order to have
			% a vector of latencies of this form: e.g.
			% hc_L *tof_R hc_R *tof_L hc_L *tof_R hc_R *tof_L hc_L *tof_R
			[strideStart,ord]		= sort([parentData.Events(gaitEventGroups).samples]);

			% extract event names
			for evIdx = find(gaitEventGroups)
			
					evLabels{:,evIdx} = repmat({parentData.Events(evIdx).label},...
							[1 numel(parentData.Events(evIdx).samples)]);		

			end

			% re-order event names accordingly
			evNames 	= [evLabels{:}];
			evNames 	= evNames(ord);

			% we filter those events that 
			peakVelocIdx 	= find(~cellfun(@isempty,regexp(evNames,'peak')));
			peakVelocMask = peakVelocIdx - 2 > 0 & peakVelocIdx +2 <= numel(evNames);


	
			% count how many strides we have recorded
%			nStrideLeft = sum(strcmp(evNames,'heelcontact_L'))-1; 
%			nStrideRight= sum(strcmp(evNames,'heelcontact_R'))-1; 

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

			% prepare plot order
			nRowsInPlot		= nStrideLeft+nStrideRight+1;
			plotOrder 		= mat2cell(reshape(1:nRowsInPlot*4,4,nRowsInPlot)',...
																ones(nRowsInPlot,1),[2 2]);

			% we have to correct the event adding the offset
			% since they are referred to the 0 of the raw data 
			evOffset 			= round(walkingStruct.Time(1)*fs);
			strideStart 	= strideStart - evOffset;

			% we create the normalized stride time vector
			referenceTimeVector = -1:1/fs:(1-1/fs);
			doubleSuppDur 			= floor(0.19*fs);
			singleSuppDur 			= floor(0.4*fs);

			referenceStance 		= 400 + [ -doubleSuppDur-singleSuppDur ...
																				-singleSuppDur 0 +doubleSuppDur ...
																		+doubleSuppDur+singleSuppDur ];

			referenceVector			= [1 referenceStance 800];
			plotIdx = 1;

			for strideIdx = peakVelocIdx(peakVelocMask) 

					% if the stride contains bad steps 
					% we skip it and continue to the next
					if ismember(plotIdx,strideRej)
							plotIdx = plotIdx + 1;
							continue;
					end

					%								[ idx ] 	
					% (hc_L) *tof_R  hc_R   *tof_L (hc_L) 
					%    ^ start-2		|t0				      ^ start + 2 == end stride
					timeWindow = strideStart(strideIdx) + (-399:400);
					freqMask 	 = walkingStruct.Freqs > 6;
					dataTF 		 = walkingStruct.TF(:,timeWindow,freqMask).*1e12;
					normFactor = repmat(mean(dataTF,2),[1 numel(timeWindow) 1]);
					dataTF 		 = (dataTF-normFactor)./normFactor;

					fprintf('[%d]',strideIdx);
					fprintf('%s ',evNames{(-2:2) + strideIdx});
					fprintf('\n');
					strideRaw	 = signals(:,timeWindow)'.*1e6;
					f 				 = walkingStruct.Freqs(freqMask);	

					% then create the time-warping vector
					originalVector = [1 (strideStart((-2:1:2) + strideIdx)...
															- timeWindow(1))  800];

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

							zLimit 	= [-1 1];
							tAxis		= referenceTimeVector;
							tEvAxis = repmat(referenceTimeVector(referenceStance([2 3 4])),2, 1);

					else
							finalTF 	= dataTF.*1e12;
							finalRaw	= strideRaw';
							zLimit 		= [min(finalTF(:)) max(finalTF(:))];
							tAxis			= (-399:400)/fs;
							tEvAxis 	= repmat(originalVector([3 4 5])./400,2, 1);
					end

					footLabel 		 	= regexp(evNames(strideIdx),'[L|R]','match');
					footLabel				= footLabel{:}{:}; % this is the label of the central hc 

					plotOrderStride	= plotOrder(plotIdx,:);

					if strcmp(footLabel,'L')
							% left foot swing => central hc_L
							%  => stnContra is rightSTN == idx 1
							controLatIdx = 1;
							if strcmp(mostAffSide,'L')
								% if STN- is L => this swing is relative to 
								% STN+ => plot on right side of page
								plotOrderStride	= plotOrderStride{2};
								stnIdx					= 2;
							else 
								% if STN- is R => this swing is relative
								% to STN- => plot on left side of page
								plotOrderStride = plotOrderStride{1};
								stnIdx					= 1;
							end
					else 
							% right foot swing
							controLatIdx = 2;
							if strcmp(mostAffSide,'L')
								% if STN- is R => this swing is relative
								% to STN- => plot on left side of page
								plotOrderStride	= plotOrderStride{1};
								stnIdx					= 1;
							else 
								% if STN- is L => this swing is relative to 
								% STN+ => plot on right side of page
								plotOrderStride = plotOrderStride{2};
								stnIdx					= 2;
							end
					end
					
					% we save for each subject the time-warped ERSD normalized over***
					stnResults{subjectIdx,stnIdx} = cat(1,stnResults{subjectIdx,stnIdx},...
																								finalTF(controLatIdx,:,:));

					finalTF = finalTF(:,referenceStance(1):referenceStance(end),:);
					time 		= referenceTimeVector(referenceStance(1):referenceStance(end));
					
					if sProcess.options.saveOutput.Value
							% save each stance in separate file
							% get the output study
							iStudy 						= sInputs(fileIdx).iStudy;
							DataMat 					= walkingStruct;
							DataMat.Freqs			= f;
							DataMat.Time			= time;
							DataMat.TF        = finalTF;
							DataMat.DataType  = 'data';
							DataMat.Comment		= sprintf('Stride %s (#%d)',footLabel,plotIdx);
										
							% Create a default output filename 
							OutputFiles{fileIdx} = bst_process('GetNewFilename', ...
									fileparts(sInputs(fileIdx).FileName), 'timefreq');

							% Save on disk
							save(OutputFiles{fileIdx}, '-struct', 'DataMat');

							% Register in database
							db_add_data(iStudy, OutputFiles{fileIdx}, DataMat);

					end
					plotIdx = plotIdx + 1;
					clear finalTF;

			end % for stride

%			conditionString = sInputs(fileIdx).Condition;
%			trialIdx = str2double(cell2mat(regexp(conditionString(strfind(conditionString,'trial'):...
%					strfind(conditionString,'trial')+7),'\d+','match')));
%
%			folder = cell2mat(lower(regexp(sInputs(fileIdx).Comment,...
%					'(Original|Normalized|RestingState)','match')));
%
%			fname = fullfile('~/Desktop','walking',folder,'wavelet',sInputs(fileIdx).SubjectName,...
%					sprintf('trial%d.png',trialIdx));
%			print(f1,'-dpng',fname);

			clear finalTF;

	end % for sInputs files

	% stnResults holds ERSD for each stride divide as STN-/+
	stnMeans = cellfun(@mean,stnResults,'UniformOutput',false);
	f2 = figure('papertype','a4','paperposition',[0 0 1 1],...
								'paperunits','normalized','paperorientation',...
								'portrait','Visible','on');
	
	folder = 'original';

	for ii = 1:nSubjects

			% this is the STN- 
			subplot(nSubjects,2,2*(ii-1)+1,'NextPlot','add')
			imagesc(tAxis,f,squeeze(stnMeans{ii,1})',zLimit);
			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
			axis xy;
			set(gca,'XTickLabel',[]);
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			ylim([6 80]);
			
			subplot(nSubjects,2,2*(ii-1)+2,'NextPlot','add')
			imagesc(tAxis,f,squeeze(stnMeans{ii,2})',zLimit);
			plot(tEvAxis,repmat([min(f);max(f)],[1 3]),'k--');
			xlim([min(tEvAxis(:)) max(tEvAxis(:))]);
			axis xy;
			set(gca,'XTickLabel',[]);
			ylim([6 80]);

	end
	annotation('textbox',[0.30,0.950,0.1,0.05],'String','STN-','LineStyle','None');
	annotation('textbox',[0.70,0.950,0.1,0.05],'String','STN+','LineStyle','None');

	fname = fullfile('/home','lgabri','Desktop',...
										'walking',folder,'wavelet','avgERSD.ps');

	print(f2,'-dpsc2',fname);

end % function
