function varargout = process_customavg( varargin )
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
    sProcess.Comment     = 'Split Events for Gait Ini';
    sProcess.FileTag     = '| imported';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Pre-process';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw','data'};
    sProcess.OutputTypes = {'raw','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

		sProcess.options.Time.Comment = ' Time Range ';
		sProcess.options.Time.Type 		= 'value';
		sProcess.options.Time.Value		= {300, 'ms', 2};
   
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
    
		% ===== LOAD THE DATA =====
		% Read the first file in the list, to initialize the loop
		DataMat = in_bst(sInputs(1).FileName, [], 0);

		Fs = 1/mean(diff(DataMat.Time));

		ChannelFlag = DataMat.ChannelFlag;
       
    % ===== PROCESS =====
		% this are costum events for gait ini only study
		% light flash and stance_foot_off for allt the trials
		nTrials				= numel(DataMat.Events(1).samples);
		nEvents				= numel(DataMat.Events);
		eventOnsets 	= [DataMat.Events.samples];
		eventOnsets 	= reshape(eventOnsets',[nTrials,nEvents])';

		% {'flash','reaction','stance_foot'}

		windowSamples	= round(sProcess.options.Time.Value{1} * Fs);

		for trial = 1:size(eventOnsets,2)

				% extract a corresponding amount of time-points 
				% during standing
				flashTimePoint  = eventOnsets(1,trial);
				standingSamples = DataMat.F(:,flashTimePoint-windowSamples:flashTimePoint);
				standingTime 		= DataMat.Time(flashTimePoint-windowSamples:...
																								flashTimePoint);

				% reaction time
				movementOnsets = eventOnsets(2,trial);
				reactionSamples= DataMat.F(:,movementOnsets-windowSamples:movementOnsets);
				reactionTime 	 = DataMat.Time(movementOnsets-windowSamples:movementOnsets);

				% during supposed gait ini
				gaitIniOnset	 = eventOnsets(3,trial);
				gaitIniSamples = DataMat.F(:,gaitIniOnset-windowSamples:gaitIniOnset);
				gaitIniTime 	 = DataMat.Time(gaitIniOnset-windowSamples:gaitIniOnset);

				% and during walking
				stanceFootOffTimePoint = eventOnsets(3,trial)+1000;
				walkingSamples = DataMat.F(:,stanceFootOffTimePoint-windowSamples:...
																			stanceFootOffTimePoint);
				walkingTime    = DataMat.Time(stanceFootOffTimePoint-windowSamples:...
																			stanceFootOffTimePoint);

				% ===== SAVE THE RESULTS =====
				% Get the output study (pick the one from the first file)
				% these are common to all output files
				iStudy 									= sInputs(1).iStudy;
				theStudy 								= bst_get('Study',iStudy);
				theProtocol							= bst_get('ProtocolInfo');
				[p f e] 								= fileparts(theStudy.FileName);
				parentPath							= fullfile(theProtocol.STUDIES,p);
				DataMatOut 							= db_template('datamat');
				DataMatOut.ChannelFlag 	= ChannelFlag;   
				DataMatOut.DataType    	= 'data';
				DataMatOut.nAvg        	= 1;         

				DataMatOut.F           	= standingSamples;
				DataMatOut.Comment     	= strcat('Standing (#',num2str(trial),')');
				DataMatOut.Time        	= standingTime;
		
				filenameOut 						= fullfile(parentPath,...
																		strjoin({'data_standing',...
																		strcat('trial',num2str(trial),'.mat')},'_'));

				OutputFiles = [OutputFiles {filenameOut}];

				% Save on disk
				save(filenameOut, '-struct', 'DataMatOut');
				% Register in database
				db_add_data(iStudy, filenameOut, DataMatOut);
				
				DataMatOut.F           = walkingSamples;
				DataMatOut.Time        = walkingTime;
				DataMatOut.Comment     = strcat('Walking (#',num2str(trial),')');

				filenameOut 					 = fullfile(parentPath,...
																		strjoin({'data_walking',...
																		strcat('trial',num2str(trial),'.mat')},'_'));
				OutputFiles = [OutputFiles {filenameOut}];

				% Save on disk
				save(filenameOut, '-struct', 'DataMatOut');
				% Register in database
				db_add_data(iStudy, filenameOut, DataMatOut);

				DataMatOut.F           = gaitIniSamples;
				DataMatOut.Time        = gaitIniTime;
				DataMatOut.Comment     = strcat('GaitIni (#',num2str(trial),')');

				filenameOut 						= fullfile(parentPath,...
																		strjoin({'data_gaitini',...
																		strcat('trial',num2str(trial),'.mat')},'_'));

				% Save on disk
				save(filenameOut, '-struct', 'DataMatOut');
				OutputFiles = [OutputFiles {filenameOut}];

				% Register in database
				db_add_data(iStudy, filenameOut, DataMatOut);

				DataMatOut.F           	= reactionSamples;
				DataMatOut.Comment     	= strcat('Reaction (#',num2str(trial),')');
				DataMatOut.Time        	= reactionTime;
		
				filenameOut 						= fullfile(parentPath,...
																		strjoin({'data_reaction',...
																		strcat('trial',num2str(trial),'.mat')},'_'));

				OutputFiles = [OutputFiles {filenameOut}];

				% Save on disk
				save(filenameOut, '-struct', 'DataMatOut');
				% Register in database
				db_add_data(iStudy, filenameOut, DataMatOut);

				standingSamples = [];
				standingTime = [];
				walkingSamples = [];
				walkingTime = [];
				gaitIniSamples = [];
				gaitIniTime = [];

		end


end
