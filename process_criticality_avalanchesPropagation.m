function varargout = process_criticality_avalanchesPropagation( varargin )

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
varargout =  {};
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Avalanches Propagation';
sProcess.FileTag     = '__';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Criticality';
sProcess.Index       = 901;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'timefreq'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(~, sInputs) %#ok<DEFNU>

% init output PLV
DataMat = db_template('timefreqmat');

for fileIdx = 1:numel(sInputs)
    % read one file at time
    data = in_bst_data(sInputs(fileIdx).FileName);
    channels = in_bst_channel(sInputs(fileIdx).ChannelFile);
    iChannels = channel_find(channels.Channel, 'SEEG');
    
    % [plv, nPLV, amp] = computePhaseMetric(data,nPermutation,method)
    avalanchesProp = avalanchesPropagation(20, 1.4, data.F(iChannels,:),1000);
    
    % ===== SAVE THE RESULTS =====
    % % Get the output study (pick the one from the first file)
    DataMat.TF			= avalanchesProp;
    
    iStudy 				= sInputs(fileIdx).iStudy;
    DataMat.Comment     = 'avalanchesPropagation';
    DataMat.ChannelFlag = data.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
    
    DataMat.Time		= [min(data.Time) max(data.Time)];
    DataMat.Freqs       = [1 2];
    DataMat.DataType    = 'data';
    DataMat.Method      = 'avalanches';
    DataMat.DataFile    = sInputs(fileIdx).FileName;
    DataMat.Measure     = 'other';
    DataMat.RefRowNames = [{channels.Channel.Name}];
    DataMat.RowNames    = [{channels.Channel.Name}];
    DataMat.nAvg        = length(sInputs); % Number of epochs that were averaged to get this file
    % Create a default output filename
    OutputFiles{1} = bst_process('GetNewFilename', ...
        fileparts(sInputs(fileIdx).FileName), 'timefreq_critic_avalProp');
    
    % Save on disk
    save(OutputFiles{1}, '-struct', 'DataMat');
    
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, DataMat);
end

end

%% Avalanches propagation
function [avalanchesPropMatrix, lengthAva] = avalanchesPropagation(timebin,tshold,Data,srate)
%
%

[row,col]=size(Data);
newcol = floor(1000*col/srate/timebin); %max number of columns divisible by time bins
npeaks = zeros(row,floor(timebin/1000*srate)*newcol);
active_channels=cell(1,newcol); %initialize active channels matrix

npeaks=fast_findpeaks(Data(:,1:size(npeaks,2)),tshold);
npeaks = reshape(npeaks,[row,floor(timebin/1000*srate),newcol]); %e.g. 101 x 40 x 15000
for tb=1:newcol  %save all the active channels (counted ones) for each timebin
    [active_channels{tb},~]=find(npeaks(:,:,tb));
end

%sum by columns and rows --> number of peaks within a single time bin
npeaks = sum((sum(npeaks,2) > 0),1); %e.g. 1 x 1 x 15000
avalanches=squeeze(npeaks)';
%zero padding to allow the first avalanche computing
avalanches=[0 avalanches];
index=find(avalanches==0); %find all the zeros in the avalanches vector
dimension_avalanches=diff(index); %dimensions of all avalanches
avalanches_propagation=zeros(row,row);

for i=1:length(dimension_avalanches)-1 %this is the number of avalanches
    if dimension_avalanches(i)>2 %if the avalanche involves at least two channels --> distance between zeros at least 3
        %save a ch*ch matrix with the columns which are the first channel
        %(SOURCE) and the rows which are the second (TARGET)
        avalanches_propagation(active_channels{index(i)+1},active_channels{index(i)})=avalanches_propagation(active_channels{index(i)+1},active_channels{index(i)})+1;
    end
end
%divide for the number of avalanches ---> probability
avalanchesPropMatrix=avalanches_propagation/length(dimension_avalanches);
lengthAva=length(dimension_avalanches);
end


