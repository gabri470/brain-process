function varargout = process_mtspectrum( varargin )
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
    sProcess.Comment     = 'MT Spectrum';
    sProcess.FileTag     = '| extracted';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Chronux';
    sProcess.Index       = 1001;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw','data'};
    sProcess.OutputTypes = {'timefreq','data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;


end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
    
    for i = 1:length(sInputs)

        % Read the file #i
        DataMatIn = in_bst(sInputs(i).FileName, [], 0);

		for ch = 1:DataMatIn.F.header.NumberOfChannels
			
			F = in_fread_brainamp(DataMatIn.F,fopen(DataMatIn.F.filename));
			tw = 4;
			[mtspect, f] = pmtm(F(ch,5*422:end),tw,...
				2^nextpow2(size(F,2)),422);
		
			A(ch,:) = mtspect;
	
		end; 

		chLabels = DataMatIn.F.header.chnames;
		chLabels = cellfun(@(x) x(1:end-5),chLabels,'UniformOutput',false);

		
		% ===== SAVE THE RESULTS =====
		% Get the output study (pick the one from the first file)
		iStudy 				= sInputs(i).iStudy;
		% Create a new data file structure
		DataMat 			= db_template('timefreq');
		DataMat.TF          = permute(A,[1 3 2]);
		DataMat.Comment     = 'PSD mtspect';
		DataMat.ChannelFlag = DataMatIn.ChannelFlag;
		DataMat.DataType 	= 'data';
		DataMat.Time        = DataMatIn.Time;
		DataMat.DataFile	= sInputs(i).FileName
		DataMat.Freqs		= f;
		DataMat.Method 		= 'psd';
		DataMat.Measure		= 'power';
		DataMat.RowNames 	= chLabels;

		% Create a default output filename 
		OutputFiles{1} = bst_process('GetNewFilename', ...
			fileparts(sInputs(i).FileName), 'timefreq_psd_mt');

		% Save on disk
		save(OutputFiles{1}, '-struct', 'DataMat');
		% Register in database
		db_add_data(iStudy, OutputFiles{1}, DataMat);

    end
end





