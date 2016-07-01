
function varargout = process_import_to_db( varargin )
% PROCESS_RESAMPLE: Resample matrix with a new sampling frequency.
%
% USAGE:      sProcess = process_resample('GetDescription')
%               sInput = process_resample('Run', sProcess, sInput, method)
%               sInput = process_resample('Run', sProcess, sInput)
%        [x, time_out] = process_resample('Compute', x, time_in, NewRate, method)
%        [x, time_out] = process_resample('Compute', x, time_in, NewRate)
%        [x,Pfac,Qfac] = process_resample('ResampleCascade', x, NewRate, OldRate, Method)

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
% Authors: Francois Tadel, 2010-2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Import to db';
    sProcess.FileTag     = '| imported';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'Pre-process';
    sProcess.Index       = 2001;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw','data'};
    sProcess.OutputTypes = {'raw','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Default values for some options
    sProcess.processDim  = 1;    % Process channel by channel
    sProcess.isSeparator = 1;
    
    % Definition of the options
    % === Resample frequency
    sProcess.options.splitlength.Comment = 'Split in time blocks of:  ';
    sProcess.options.splitlength.Type    = 'value';
    sProcess.options.splitlength.Value   = {2,'s',2};
    
    sProcess.options.split.Comment = 'Split in blocks';
    sProcess.options.split.Type    = 'checkbox';
    sProcess.options.split.Value   = 0;
    
    sProcess.options.time.Comment = 'Time Range';
    sProcess.options.time.Type    = 'timewindow';
    sProcess.options.time.Value   = {[],'s',2};

    sProcess.options.ssp.Comment = 'Apply ssp';
    sProcess.options.ssp.Type    = 'checkbox';
    sProcess.options.ssp.Value   = 1;

    sProcess.options.baseline.Comment = 'Remove DC';
    sProcess.options.baseline.Type    = 'checkbox';
    sProcess.options.baseline.Value   = 1;

    sProcess.options.DCtime.Comment = 'Time Range';
    sProcess.options.DCtime.Type    = 'timewindow';
    sProcess.options.DCtime.Value    = {[],'s',2};
    
    sProcess.options.createcond.Comment = 'Create new condition';
    sProcess.options.createcond.Type    = 'checkbox';
    sProcess.options.createcond.Value   = 1;
    

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>

ImportOptions = db_template('ImportOptions');
ImportOptions.ImportMode = 'Time';
ImportOptions.UseEvents = 0;
ImportOptions.TimeRange=      sProcess.options.time.Value{1};
ImportOptions.EventsTimeRange=[];
ImportOptions.GetAllEpochs=   0;
ImportOptions.iEpochs=        1;
ImportOptions.SplitRaw=       sProcess.options.split.Value;
ImportOptions.SplitLength=    sProcess.options.splitlength.Value{1};
ImportOptions.Resample=       0;
ImportOptions.ResampleFreq=   0;
ImportOptions.UseCtfComp=     0;
ImportOptions.UseSsp=         sProcess.options.ssp.Value;
if(sProcess.options.baseline.Value)
ImportOptions.RemoveBaseline= 'all';
else
ImportOptions.RemoveBaseline= 'no';
end
ImportOptions.BaselineRange=  sProcess.options.DCtime.Value{1};
ImportOptions.events=         [];
ImportOptions.CreateConditions= sProcess.options.createcond.Value;
ImportOptions.ChannelReplace= 1;
ImportOptions.ChannelAlign =   1;
ImportOptions.IgnoreShortEpochs = 1;
ImportOptions.EventsMode =     'ask';
ImportOptions.EventsTrackMode ='ask';
ImportOptions.DisplayMessages=0;

DataFile = sInput.FileName;

NewFiles = import_raw_to_db( DataFile,ImportOptions );

end


function NewFiles = import_raw_to_db( DataFile ,ImportOptions)
% IMPORT_RAW_TO_DB: Import in the database some blocks of recordings from a continuous file already linked to the database.
%
% USAGE:  NewFiles = import_raw_to_db( DataFile )

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
% Authors: Francois Tadel, 2011-2014


% ===== GET FILE INFO =====
% Get study description
[sStudy, iStudy, iData] = bst_get('DataFile', DataFile);
if isempty(sStudy)
    error('File is not registered in the database.');
end
% Is it a "link to raw file" or not
isRaw = strcmpi(sStudy.Data(iData).DataType, 'raw');
% Get subject index
[sSubject, iSubject] = bst_get('Subject', sStudy.BrainStormSubject);
% Progress bar
bst_progress('start', 'Import raw file', 'Processing file header...');
% Read file descriptor
DataMat = in_bst_data(DataFile);
% Read channel file
ChannelFile = bst_get('ChannelFileForStudy', DataFile);
ChannelMat = in_bst_channel(ChannelFile);
% Get sFile structure
if isRaw
    sFile = DataMat.F;
else
    sFile = in_fopen(DataFile, 'BST-DATA');
end

% Import file
    

NewFiles = import_data(sFile, ChannelMat, sFile.format, [], iSubject,ImportOptions);


end