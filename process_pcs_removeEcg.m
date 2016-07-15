function varargout = process_pcs_removeEcg( varargin )
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
    sProcess.Comment     = '[PCS] Remove Ecg';
    sProcess.FileTag     = '| ecg_pruned';
    sProcess.Category    = 'Filter';
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
function sInputs = Run(sProcess, sInputs) %#ok<DEFNU>


		for file = 1:numel(sInputs)

				% check whether we have already detected the cardiac artefact
				DataMat = in_bst(sInputs.FileName);
				ChannelMat = in_bst_channel(sInputs.ChannelFile);
				iChannels = channel_find(ChannelMat.Channel,'SEEG');

				evCardiacIdx = find(ismember({DataMat.Events.label},'cardiac'));

				fs = round(1/mean(diff(DataMat.Time)));

				offset = round(DataMat.Time(1)*fs);

				evCardiacSamples = DataMat.Events(evCardiacIdx).samples - offset;
				timeWnd= bsxfun(@plus,(round(-100*.4):round(500*.4)),evCardiacSamples')';
				kernel = tukeywin(size(timeWnd,1));


				signals = sInputs.A(iChannels,timeWnd(:));
				signals = reshape(signals',[241 numel(evCardiacSamples) 2]);
				artefact = squeeze(mean(signals,2));

				artefact= bsxfun(@times,artefact,kernel(:));

				sInputs.A(iChannels,timeWnd) = sInputs.A(iChannels,timeWnd) - ...
						repmat(artefact',1,numel(evCardiacSamples));

		end

end
