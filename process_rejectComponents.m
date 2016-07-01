function varargout = process_rejectComponents( varargin )
		% PROCESS_REJECTCOMPONENTS
		
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
		% Authors:

		macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
		% Description the process
		sProcess.Comment     = 'Reject ICA Components';
		sProcess.FileTag     = '';
		sProcess.Category    = 'Custom';
		sProcess.SubGroup    = 'Pre-process';
		sProcess.Index       = 1000;

		% Definition of the input accepted by this process
		sProcess.InputTypes  = {'data'};
		sProcess.OutputTypes = {'data'};
		sProcess.nInputs     = 1;
		sProcess.nMinFiles   = 1;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
		Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

		% ===== LOAD THE DATA =====
		ICAdata = in_bst_data(sInputs(1).FileName);
      			
		% ===== PROCESS =====
		Compute(ICAdata,sInputs);

        OutputFiles = {};

end


function  Compute(ICAdata,sInputs)

		bst_removeICA_GUI(ICAdata,sInputs);

end
