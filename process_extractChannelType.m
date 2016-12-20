function varargout = process_extractChannelType( varargin )
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
    sProcess.Comment     = 'Extract Channels Type';
    sProcess.FileTag     = '| extracted';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Examples';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw','data'};
    sProcess.OutputTypes = {'raw','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
        % Sensor types
        sProcess.options.sensortypes.Comment = 'Sensor types or names: ';
        sProcess.options.sensortypes.Type    = 'text';
        sProcess.options.sensortypes.Value   = 'EEG';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
    % Get option values
    example3 = sProcess.options.sensortypes.Value;
    
    % ===== LOAD THE DATA =====
    % Read the first file in the list, to initialize the loop
    DataMat = in_bst(sInputs(1).FileName, [], 0);
    epochSize = size(DataMat.F);
     % Initialize the load matrix: [Nchannels x Ntime x Nepochs]
    AllMat = zeros(epochSize(1), epochSize(2)*length(sInputs));

	ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
	% Get channels to process
	iChannels = channel_find(ChannelMat.Channel, ...
		sProcess.options.sensortypes.Value);

    % Reading all the input files in a big matrix
    for i = 1:length(sInputs)
        % Read the file #i
        DataMat = in_bst(sInputs(i).FileName, [], 0);
        % Check the dimensions of the recordings matrix in this file
        if ~isequal(size(DataMat.F), epochSize)
            % Add an error message to the report
            bst_report('Error', sProcess, sInputs, 'One file has a different number of channels or a different number of time samples.');
            % Stop the process
            return;
        end

		Time([1:epochSize(2)]+(i-1)*epochSize(2)) = DataMat.Time;

        % Add the current file in the big load matrix
        AllMat(:,[1:epochSize(2)]+(i-1)*epochSize(2)) = DataMat.F;
    end
    
    % ===== PROCESS =====
    % Just doing a simple average of the trials, can be replaced with anything
    AllMat = AllMat(iChannels,:);
    
    % ===== SAVE THE RESULTS =====
    % Get the output study (pick the one from the first file)
    iStudy = sInputs(1).iStudy;
    % Create a new data file structure
    DataMat 						= db_template('datamat');
    DataMat.F           = AllMat;
    DataMat.Comment     = ['Channels extracted ' sProcess.options.sensortypes.Value];
    DataMat.ChannelFlag = ones(length(iChannels), 1);   % List of good/bad channels (1=good, -1=bad)
    DataMat.Time        = Time;
    DataMat.DataType    = 'data';
    DataMat.nAvg        = length(sInputs);         % Number of epochs that were averaged to get this file
    % Create a default output filename 
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'data_extracted_');
    % Save on disk
    save(OutputFiles{1}, '-struct', 'DataMat');
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, DataMat);
end





