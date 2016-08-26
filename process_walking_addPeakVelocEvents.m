function varargout = process_addPeakVelocEvents( varargin )
% PROCESS_ADDPEAKVELOCEVENTS: 
%
% DESCRIPTION: 


% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2015 University of Southern California & McGill University
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
% Authors: Gabriele Arnulfo
%
macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Add Peak Velocity Events';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Walking';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results', 'timefreq', 'matrix'};
    sProcess.OutputTypes = {'data', 'results', 'timefreq', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Default values for some option
    sProcess.isSourceAbsolute = 1;
    sProcess.processDim       = 1;    % Process channel by channel

		sProcess.options.markerString.Comment = 'Marker Name: ';
		sProcess.options.markerString.Type    = 'text';
		sProcess.options.markerString.Value   = 'ankle';


end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
		Comment = sProcess.Comment;
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>

    % Get inputs
		DATA_FOLDER = '/media/lgabri/My Passport/gait-v2';
		subjectIdx = sInput.SubjectName;
		conditionString = sInput.Condition;


		sMat = in_bst_channel(sInput.ChannelFile);

		trialIdx = cell2mat(regexp(conditionString,'trial\d+','match'));

		dataMat = in_bst(sInput.FileName);

		% define velocity data file
		velocityFileName = fullfile(DATA_FOLDER,subjectIdx,strjoin({subjectIdx,'walking','off',...
																trialIdx,sProcess.options.markerString.Value,'velocity.csv'},'_'));

		rawEventFileName = fullfile(DATA_FOLDER,subjectIdx,strjoin({subjectIdx,'walking','off',...
																									'raw','events.csv'},'_'));

		% read data 
		velocData = importdata(velocityFileName);

		% read events without TENS
		rawEventData = importdata(rawEventFileName);

		% get the joint order from variable/column names in velocData
		markerNames = velocData.colheaders;

		% guess the feet order from variable names
		leftMarkerIdx	  = ~cellfun(@isempty,regexp(markerNames, '_L'));
		rightMarkerIdx	= ~cellfun(@isempty,regexp(markerNames, '_R'));
			
		% extract the first event from all available which is 
		% the first heel_contact_left defined in the raw
		% time axis. This point will be used as reference point 
		% in the sense that we compute event delay from the first
		% detected heel contact
		trialIdxRowInMat = str2double(regexp(trialIdx,'\d+','match'));

		% get column names for raw events
		if ~isfield(rawEventData,'colheaders')
				error(' Something wrong with column names in raw events file');
		else
				rawEvMask = ~cellfun(@isempty,regexp(rawEventData.colheaders,'(heel|toe)'));
		end

		rawOffsetInSeconds = min(rawEventData.data(trialIdxRowInMat,rawEvMask));

		evGroupNames = {dataMat.Events.label};
		gaitEventGroups = ~cellfun(@isempty,regexp(evGroupNames,'(heel|toe)'));

		referenceOffsetInSeconds = min([dataMat.Events(gaitEventGroups).times]);

		thresholds = [1.5 1.5];
			

		% detect peak of the velocity
		[~,lhLoc] = findpeaks(velocData.data(:,leftMarkerIdx),'MinPeakHeight',...
															thresholds(1),'MinPeakDistance',10);
		[~,rhLoc] = findpeaks(velocData.data(:,rightMarkerIdx),'MinPeakHeight',...
															thresholds(2),'MinPeakDistance',10);

		times1 = velocData.data(lhLoc,1)'-rawOffsetInSeconds+referenceOffsetInSeconds;
		times2 = velocData.data(rhLoc,1)'-rawOffsetInSeconds+referenceOffsetInSeconds;

		newEvent = struct('label',{'peakVeloc_L','peakVeloc_R'},...
											'times',{times1,times2},...
											'samples',{round(times1.*400),round(times2.*400)},...
											'color',{[.75 0 .75],[1 1 0]},...
											'epochs',{ones(1,numel(lhLoc)),ones(1,numel(rhLoc))},...
											'reactTimes',[],...
											'select',1);

		dataMat.Events = [dataMat.Events, newEvent];

		tibMask = ~cellfun(@isempty,regexp({sMat.Channel.Name},'(sol|[L|R]S_emg)'));

		figure, plot(dataMat.Time,dataMat.F(tibMask,:)), 
				hold on, 

		tmpEvents = [dataMat.Events(gaitEventGroups).times];
		peakLeftEvents = [dataMat.Events(end-1).times];
		peakRightEvents = [dataMat.Events(end).times];

		plot(repmat(tmpEvents(:)',[2 1]),repmat([min(dataMat.F(:));.5],...
				[1 numel(tmpEvents)]),'k--');

		plot(repmat(peakLeftEvents(:)',[2 1]),repmat([min(dataMat.F(:));.5],...
				[1 numel(peakLeftEvents)]),'b--');

		plot(repmat(peakRightEvents(:)',[2 1]),repmat([min(dataMat.F(:));.5],...
				[1 numel(peakRightEvents)]),'g--');


		plot(velocData.data(:,1)-rawOffsetInSeconds+referenceOffsetInSeconds,...
					velocData.data(:,2:end).*1e-1,'LineWidth',3)


		STUDY_DIR = bst_get('ProtocolInfo');
		bst_save(fullfile(STUDY_DIR.STUDIES,sInput.FileName),dataMat);

end
