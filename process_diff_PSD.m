function varargout = process_diff_PSD( varargin )
% PROCESS_DIFF: Difference of each couple of samples (A-B).

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
% LIABILITY OR RESPONSIBILITY FOR THE USEt OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2010

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Compare PSD (B-A)/A';
    sProcess.FileTag     = [];
    sProcess.Category    = 'Filter2';
    sProcess.SubGroup    = 'Other';
    sProcess.Index       = 701;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results', 'timefreq', 'matrix'};
    sProcess.OutputTypes = {'data', 'results', 'timefreq', 'matrix'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;
    sProcess.isPaired    = 0;
    % Default values for some options
    sProcess.isSourceAbsolute = 1;




end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
    % Absolute values 
    if isfield(sProcess.options, 'source_abs') && sProcess.options.source_abs.Value
        Comment = [Comment, ', abs'];
    end
end


%% ===== RUN =====
function sOutput = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};

    
    % ===== LOAD THE DATA =====
    % Read the first file in the list, to initialize the loop
    DataMatA = in_bst(sInputsA(1).FileName, [], 0);
    epochSizeA = size(DataMatA.TF);
    DataMatB = in_bst(sInputsB(1).FileName, [], 0);
    epochSizeB = size(DataMatB.TF);
     % Initialize the load matrix: [Nchannels x Ntime x Nepochs]
    AllMatA = DataMatA.TF;
    AllMatB = DataMatB.TF;

	ChannelMat = in_bst_channel(sInputsA(1).ChannelFile);


    % ===== PROCESS =====
    % Just doing a simple average of the trials, can be replaced with anything
    AllMatA = squeeze(AllMatA);
    AllMatB = squeeze(AllMatB);
    
%	AllMatA = AllMatA./repmat(sum(AllMatA,2),[1 epochSizeA(3)]);
%	AllMatB = AllMatB./repmat(sum(AllMatB,2),[1 epochSizeB(3)]);

	comp = (AllMatB-AllMatA)./AllMatA;

	figure, plot(DataMatA.Freqs,AllMatA), legend(DataMatA.RowNames),title(['PSD OFF ' sInputsA.SubjectName]),xlim([0 40]);
	hold on,plot([DataMatA.Freqs(1) DataMatA.Freqs(end)],[0 0],'--b')
	sOutput = [];
	figure, plot(DataMatA.Freqs,AllMatB), legend(DataMatA.RowNames),title(['PSD ON ' sInputsA.SubjectName]),xlim([0 40]);
	hold on,plot([DataMatA.Freqs(1) DataMatA.Freqs(end)],[0 0],'--b')
end




