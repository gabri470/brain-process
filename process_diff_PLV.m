function varargout = process_diff_ab( varargin )
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

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compare PLV (B-A)/A';
    sProcess.FileTag     = [];
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Other';
    sProcess.Index       = 701;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw' ,'timefreq'};
    sProcess.OutputTypes = {'data', 'raw' ,'timefreq' };
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
    AllMatA = DataMatA.TF./abs(DataMatA.TF);
    AllMatB = DataMatB.TF./abs(DataMatB.TF);

	ChannelMat = in_bst_channel(sInputsA(1).ChannelFile);


    % ===== PROCESS =====
    % Just doing a simple average of the trials, can be replaced with anything
	% CH x T x F	- complex
	for f = 1:epochSizeA(3)
		cPLVA(:,:,f) = AllMatA(:,:,f)*AllMatA(:,:,f)'./epochSizeA(2);
		cPLVB(:,:,f) = AllMatB(:,:,f)*AllMatB(:,:,f)'./epochSizeB(2);
	end % frequency

	PLVA  = squeeze(abs(cPLVA(1,2,:)));
	iPLVA = squeeze(abs(imag(cPLVA(1,2,:))));

	PLVB  = squeeze(abs(cPLVB(1,2,:)));
	iPLVB = squeeze(abs(imag(cPLVB(1,2,:))));

	figure,

%	comp = (PLVB-PLVA)./PLVA;
	subplot(2,1,1)
	plot(DataMatA.Freqs,PLVA,DataMatB.Freqs,PLVB), ...
		legend([{strrep(sInputsA(1).Condition,'_',' ')},...
				{strrep(sInputsB(1).Condition,'_',' ')}]);...
		title(['diff PLV ' sInputsA.SubjectName]),xlim([0 40]);

	hold on,plot([DataMatA.Freqs(1) DataMatA.Freqs(end)],[0 0],'--b')
	sOutput = [];

%	comp = (iPLVB-iPLVA)./iPLVA;
	subplot(2,1,2)
	plot(DataMatA.Freqs,iPLVA,DataMatB.Freqs,iPLVB), ...
		legend([{strrep(sInputsA(1).Condition,'_',' ')},...
				{strrep(sInputsB(1).Condition,'_',' ')}]);...
		title(['diff iPLV ' sInputsA.SubjectName]),xlim([0 40]);

	hold on,plot([DataMatA.Freqs(1) DataMatA.Freqs(end)],[0 0],'--b')
	sOutput = [];

end




