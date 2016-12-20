function varargout = process_plot_gait_figures( varargin )
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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Plot Gait Figures';
    sProcess.FileTag     = '| figure';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Custom';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
   
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
		fid = fopen('~/Desktop/gaitini.log','w');

    % Initialize returned list of files
    OutputFiles = {};
    
		conditions  = unique( [{sInputs.Condition}]);
		subjects		= unique( [{sInputs.SubjectName}]);

		possConditionNames = { 'Standing','Reaction','GaitIni','Walking'};
		possMedLevel = {'ON','OFF'};

		panelTitles = strcat( repmat(possConditionNames',[1,2]),...
				repmat(possMedLevel,[4,1]) );

		% i have to read the foot from csv
		% WUE, MEDS, TRIAL, FOOT, MOVONSET(ignore)
		[subjectNames, meds, trials, foots] = textread('/home/gabri/Gait.csv',...
				'%s %s %d %s %*f','delimiter',',','headerlines',1);

		chControlatIndices = ones(size(trials));
		chIpsilatIndices = ones(size(trials));

		chControlatIndices(cellfun(@strcmp,foots,repmat({ 'right' },size(foots)))) = 3;
		chIpsilatIndices(cellfun(@strcmp,foots,repmat({ 'left' },size(foots)))) = 3;

		% reshape as trial, med, subj
		chControlatIndices = reshape(chControlatIndices,[4,2,3]);
		chIpsilatIndices = reshape(chIpsilatIndices,[4,2,3]);

		f(1) = figure(1);
		f(2) = figure(2);
		f(3) = figure(3);

		for conditionIdx = 1:numel(conditions)
			for subjectIdx = 1:numel(subjects)

						plotIdx = find(ismember(panelTitles',conditions{conditionIdx}));		
						theCondition = find(strcmp([{sInputs.Condition}],...
								conditions{conditionIdx}));
						theSubject = find(strcmp([{sInputs.SubjectName}],...
								subjects{subjectIdx}));

						fileIdx = intersect(theCondition,theSubject);

						n = numel(fileIdx);

						for file = fileIdx
								
								DataMat = in_bst_timefreq(sInputs(file).FileName);
								parentStruct = bst_process('GetInputStruct',DataMat.DataFile);

								F = DataMat.Freqs;
								betaRange = (F>= 13 & F<= 30);

								trial 	= regexp(parentStruct.Comment,'\d+','match');
								trial		= str2num(trial{:});

								medIdx  = cellfun(@strfind, repmat(conditions(conditionIdx)...
															,[2 1]),possMedLevel','uni',false);
								medIdx	= ~cellfun(@isempty,medIdx);

								chControlatIdx = chControlatIndices(trial,medIdx,subjectIdx);
								chIpsilatIdx = chIpsilatIndices(trial,medIdx,subjectIdx);

								data = DataMat.TF;

								[betaPeakVal betaPeakIdx]	= max(data([1:3],1,betaRange),[],3);
								betaBandIndices 					= find(betaRange);
								betaPeakFrequency					= F(betaBandIndices(betaPeakIdx));
								singleSubjBetaRange 			= repmat(betaPeakFrequency(:),[1,2])...
																						+ [-5 5; -5 5; -5 5];

								F = repmat(F,[3, 1]);
								betaRange = F >= repmat(singleSubjBetaRange(:,1),[1 40]) & ...
														F <= repmat(singleSubjBetaRange(:,2),[1 40]); 

							
								
								betaBandPower(1) = squeeze(sum(data(chControlatIdx,...
																				1,betaRange(chControlatIdx,:)),3));
								betaBandPower(2) = squeeze(sum(data(chIpsilatIdx,...
																				1,betaRange(chIpsilatIdx,:)),3));
									

								fprintf(fid,'%s, %s, %s, %d, %e, %e, %e, %e\n',...
																					subjects{subjectIdx},...
																					possMedLevel{medIdx},...
																				  conditions{conditionIdx},...
																					trial,...
																					betaPeakVal([1 3]),betaBandPower);

								figure(f(subjectIdx));
								subplot(4,2,plotIdx);
								hold on;

								plot(DataMat.Freqs,squeeze(data(1,:,:)),'r');
								
								plot([betaPeakFrequency(1), betaPeakFrequency(1)],...
										[0  betaPeakVal(1)],'k--');

								plot(DataMat.Freqs,squeeze(data(3,:,:)),'b');

								plot([betaPeakFrequency(3), betaPeakFrequency(3)],...
										[0  betaPeakVal(3)],'k--');

								ylim([0 5e-12]);

								title(conditions{conditionIdx});
								box off;


						end
				end

%				data2 = data2 ./n;
%				data1 = (data1 ./n);
%				dataStd = sqrt(data2-(data1).^2);
%				upBound	= data1 + dataStd./sqrt(n);
%				loBound = data1 - dataStd./sqrt(n);
%
%				opts1 = {'EdgeColor','none',...
%								'FaceColor',[1 0.5 .5]};
%				opts2 = {'EdgeColor','none',...
%								'FaceColor',[0.5 0.5 1]};
																	
%				subplot(4,2,plotIdx);
%				hold on;
%
%				plot(DataMat.Freqs,squeeze(data1(1,:,:)),'r','LineWidth',2);
%				fill_between(DataMat.Freqs,squeeze(upBound(1,:,:)),...
%						squeeze(loBound(1,:,:)),[],opts1{:});
%				
%				plot([betaPeakFrequency(1), betaPeakFrequency(1)],...
%						[0  betaPeakVal(1)],'k--');
%
%				fill_between(DataMat.Freqs,squeeze(upBound(2,:,:)),...
%						squeeze(loBound(2,:,:)),[],opts2{:});
%				plot(DataMat.Freqs,squeeze(data1(2,:,:)),'b','LineWidth',2);
%
%				plot([betaPeakFrequency(2), betaPeakFrequency(2)],...
%						[0  betaPeakVal(2)],'k--');
%
%
%				ylim([0 5e-12]);
%
%				title(conditions{conditionIdx});
%				box off;


				end

%				print(f1,'-dpng','-r600',fullfile('~/Desktop/',...
%								strcat(subjects{subjectIdx},'.png')));
		fclose(fid);
end
