function varargout = process_compare_PSDs2( varargin )
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
    sProcess.Comment     = 'Compare PSDs pippo';
    sProcess.FileTag     = '| compared';
    sProcess.Category    = 'Filter2';
    sProcess.SubGroup    = 'Examples';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq','data'};
    sProcess.OutputTypes = {'timefreq','data'};
    sProcess.nInputs     = 2;
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
function OutputFiles = Run(sProcess, sInputsA,sInputsB) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
    % Get option values
    example3 = sProcess.options.sensortypes.Value{:};
    
    % ===== LOAD THE DATA =====
    % Read the first file in the list, to initialize the loop
    DataMatA = in_bst(sInputsA(1).FileName, [], 0);
    epochSizeA = size(DataMat.F);
    DataMatB = in_bst(sInputsB(1).FileName, [], 0);
    epochSizeB = size(DataMatB.F);
     % Initialize the load matrix: [Nchannels x Ntime x Nepochs]
    AllMatA = DataMatA.F;
    AllMatB = DataMatB.F;

	ChannelMat = in_bst_channel(sInputsA(iFile).ChannelFile);
	% Get channels to process
	iChannels = channel_find(ChannelMat.Channel, ...
		sProcess.options.sensortypes.Value);

    % ===== PROCESS =====
    % Just doing a simple average of the trials, can be replaced with anything
    AllMatA = AllMatA(iChannels,:);
    AllMatB = AllMatB(iChannels,:);
    
	AllMatA = AllMatA./repmat(sum(AllMatA,2),[1 epochSizeA(2)]);
	AllMatB = AllMatB./repmat(sum(AllMatB,2),[1 epochSizeB(2)]);

	comp = (AllMatA-AllMatB)./AllMatB;

	plot(DataMatA.Time,comp)

    % ===== SAVE THE RESULTS =====
    % Get the output study (pick the one from the first file)
%    iStudy = sInputs(1).iStudy;
%    % Create a new data file structure
%    DataMat 			= db_template('datamat');
%    DataMat.F           = AllMat;
%    DataMat.Comment     = 'Channels extracted';
%    DataMat.ChannelFlag = ones(length(iChannels), 1);   % List of good/bad channels (1=good, -1=bad)
%    DataMat.Time        = Time;
%    DataMat.DataType    = 'data';
%    DataMat.nAvg        = length(sInputs);         % Number of epochs that were averaged to get this file
%    % Create a default output filename 
%    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'data_extracted_');
%    % Save on disk
%    save(OutputFiles{1}, '-struct', 'DataMat');
%    % Register in database
%    db_add_data(iStudy, OutputFiles{1}, DataMat);
end





