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
    sProcess.FileTag     = 'ecg_pruned';
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
				DataMat = in_bst(sInputs(file).FileName);
				ChannelMat = in_bst_channel(sInputs(file).ChannelFile);
				iChannels = channel_find(ChannelMat.Channel,'SEEG');

				fs = round(1/mean(diff(DataMat.Time)));

				% find indices of cardiac artefacts
				evCardiacIdx = ismember({DataMat.Events.label},'cardiac');

				offset = round(DataMat.Time(1)*fs);

				% compute cardiac samples 
				evCardiacSamples = DataMat.Events(evCardiacIdx).samples - offset;
				timeWnd= bsxfun(@plus,round(-.5*fs):(round(.5*fs)-1),evCardiacSamples')';

				% reject events that are outside the analyses window
				evCardiacMask = sum(timeWnd < 0,1)==0 & sum(timeWnd > max(size(sInputs(file).A)),1)==0; 
				timeWnd = timeWnd(:,evCardiacMask);
				evCardiacSamples = evCardiacSamples(evCardiacMask);


				signals = sInputs(file).A(iChannels,timeWnd(:));
				signals = reshape(signals',[size(timeWnd,1) numel(evCardiacSamples) 2]);
				artefact = squeeze(mean(signals,2));

                tk_win = tukeywin(size(timeWnd,1),0.2);
                
				timeWnd = timeWnd';

				for iCh = 1:2 
						% lfpseg is events x t
						lfpseg = squeeze(signals(:,:,iCh))';            
						meanlfpseg = artefact(:,iCh);
						clfpsegAdapt = zeros(size(timeWnd));

						for ecgIdx = 1:numel(evCardiacSamples)

								Amp = fminsearch(@(Amp)myminfun(Amp,lfpseg(ecgIdx,:)',...
										tk_win.*meanlfpseg),1);         

								clfpsegAdapt(ecgIdx,:) = lfpseg(ecgIdx,:)'- ...
										Amp*tk_win.*meanlfpseg;
						end

						sInputs(file).A(iChannels(iCh),timeWnd(:)) = clfpsegAdapt(:);
				end

		end

end


function F = myminfun(Amp,lfp,art)

	F = sqrt(sum((lfp - Amp*art).^2));
end
