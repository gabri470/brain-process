function varargout = process_filter_fir( varargin )

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
sProcess.Comment     = 'filter w fir bank';
sProcess.FileTag     = @GetFileTag;
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Pre-process';
sProcess.Index       = 61;

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

sProcess.options.minFreq.Type	 					= 'text';
sProcess.options.minFreq.Comment 				= 'Min Freq';
sProcess.options.minFreq.Value	 				= '2';

sProcess.options.scalingFactor.Type	 		= 'text';
sProcess.options.scalingFactor.Comment 	= 'Scaling Factor';
sProcess.options.scalingFactor.Value	 	= num2str(sqrt(2));

sProcess.options.bandFlatTop.Type 			= 'text';
sProcess.options.bandFlatTop.Comment 		= 'Band Flat Top';
sProcess.options.bandFlatTop.Value 			= '2';

sProcess.options.stopBand.Type 					= 'text';
sProcess.options.stopBand.Comment 			= 'Stop Band';
sProcess.options.stopBand.Value 				= '2';

sProcess.options.bandStopWidth.Type 		= 'text';
sProcess.options.bandStopWidth.Comment 	= 'Band Stop Width';
sProcess.options.bandStopWidth.Value 		= '5';

sProcess.options.nFilters.Type 					= 'text';
sProcess.options.nFilters.Comment 			= 'Number of Filters';
sProcess.options.nFilters.Value 				= '13';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
	Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

		fs = 5000;
		[filterSettings,fn,filterStrings] = generateFilterSettings(sProcess.options,fs);
		nFilters = size(filterSettings,1);
		
		totElement = numel(sInputs)*nFilters;

		pBar = bst_progress('start','Filtering','Filtering',0,totElement);


		for fileIdx = 1:numel(sInputs)
				DataMat = in_bst_data(sInputs(fileIdx).FileName);

				for fIdx = 1:nFilters
						bst_progress('inc',(((fileIdx-1)*numel(sInputs))+fIdx));

						fs = 1/mean(diff(DataMat.Time));

						[filteredDataMat] = Compute(DataMat.F',filterSettings,fIdx,fs);

						% ===== SAVE THE RESULTS =====   
						filtDataMat 					= db_template('DataMat');
						filtDataMat						= DataMat;
						iStudy 								= sInputs(fileIdx).iStudy;    
						filtDataMat.Comment   = strcat('Filtered Band-',num2str(fn(fIdx),'%.2f'));    
						filtDataMat.F					= filteredDataMat';
						filtDataMat.DataType  = 'data';    
						filtDataMat.nAvg      = 1; 
						filtDataMat 					= bst_history('add',filtDataMat,'filtering',filterStrings{fIdx});

						% Create a default output filename     
						OutputFiles{1} = bst_process('GetNewFilename', ...
							fileparts(sInputs(fileIdx).FileName), 'data_band');    
						% Save on disk    
						save(OutputFiles{1}, '-struct', 'filtDataMat');    

						% Register in database    
						db_add_data(iStudy, OutputFiles{1}, filtDataMat);
						clear filtDataMat;

				end
		end

end


function data = Compute(data,filterSettings,fIdx,fs)

		% data T x CH Input => output
		data = bandpass_FIR_filter(data,fs,filterSettings{fIdx,1});
		data = bandstop_FIR_filter(data,fs,filterSettings{fIdx,2});


end


function [filterSettings,fn,filterStrings] = generateFilterSettings(options,fs) 
%GENERATEFILTERSETTINGS wrapper for designfilt 
%	FILTERSETTINGS = GENERATEFILTERSETTINGS(SPROCESS.OPTIONS,FS) 
%	This function should translate the sProcess.options into a cellarray of fields 
%	that can be easily passed to designfilt

	nFilters 		= str2num(options.nFilters.Value);
	fn 			 		= str2num(options.minFreq.Value);
	wb			 		= str2num(options.bandFlatTop.Value);
	ws			 		= str2num(options.stopBand.Value);
	m				 		= str2num(options.scalingFactor.Value);
	stopBandWb	= str2num(options.bandStopWidth.Value);

	f 			 		= zeros(nFilters,1);
	fnyq		 		= fs/2;
	lineNoise		= 50;
	harmonics		= linspace(50,fs/2,(fs/2)/50);
	nHarmonics	= numel(harmonics);

	filterSettings = cell(nFilters,2);

	fprintf('FILTER SETTINGS:\n');
	for fIdx = 1:nFilters

			if fn > 100
					m = sqrt(sqrt(2));
			end

			fprintf('\t[%.2f - %.2f] -> %.2f <- [%.2f - %.2f]\t',(fn / ws),(fn - (wb * fn)/2),fn,...
					(fn + (wb * fn)/2),(min([0.95*(fs/2) fn * ws])));

			filterStrings{fIdx} = sprintf('[%.2f - %.2f] -> %.2f <- [%.2f - %.2f]\t',(fn / ws),...
					(fn - (wb * fn)/2),fn,(fn + (wb * fn)/2),(min([0.95*(fs/2) fn * ws])));
	
			passBandLp = (fn + (wb * fn)/2);
			passBandHp = (fn - (wb * fn)/2);
			stopBandLp = (min([0.95*(fs/2) fn * ws]));
			stopBandHp = (fn / ws);

			transitionBand = [passBandHp - stopBandHp, stopBandLp - passBandLp];

			f(fIdx) 	 = fn;
			fn 				 = fn * m;

			% [x, FiltSpec, Messages] = bandpass_FIR_filter(x, Fs, HighPass,LowPass,transitionBand,attenuation,removeDC, isMirror, Function)
			[~,filtSpec] = bandpass_FIR_filter([],fs,passBandHp,passBandLp,transitionBand,[],[],1);
			filterSettings{fIdx,1} = filtSpec; 

			%	curHarmonics contains the harmonics of 50Hz in the band of interest
			curHarmonics = harmonics(harmonics > stopBandHp & harmonics < stopBandLp);

			% [x, FiltSpec, Messages] = bandstop_FIR_filter(x, Fs, notch_freq_list,transition_band,attenuation, isMirror, Function)
			[~,filtSpec] = bandstop_FIR_filter([],fs,curHarmonics,[stopBandWb,stopBandWb],80,1);
			filterSettings{fIdx,2} = filtSpec;
					
			fprintf('\n');

	end

end


%function [filterSettings,f,filterStrings] = generateFilterSettings(options,fs) 
%%GENERATEFILTERSETTINGS wrapper for designfilt 
%%	FILTERSETTINGS = GENERATEFILTERSETTINGS(SPROCESS.OPTIONS,FS) 
%%	This function should translate the sProcess.options into a cellarray of fields 
%%	that can be easily passed to designfilt
%
%
%	nFilters 	= str2num(options.nFilters.Value);
%	fn 			 	= str2num(options.minFreq.Value);
%	wb			 	= str2num(options.bandFlatTop.Value);
%	ws			 	= str2num(options.stopBand.Value);
%	m				 	= str2num(options.scalingFactor.Value);
%
%	f 			 	= zeros(nFilters,1);
%	fnyq		 	= fs/2;
%	lineNoise	= 50;
%	harmonics	= linspace(50,fs/2,(fs/2)/50)/(fs/2);
%	tmp 			= linspace(50,fs/2,(fs/2)/50);
%	nHarmonics= numel(harmonics);
%
%	filterSettings = cell(nFilters,nHarmonics+1);
%
%	fprintf('FILTER SETTINGS:\n');
%	for fIdx = 1:nFilters
%
%			if fn > 100
%					m = sqrt(sqrt(2));
%			end
%
%			fprintf('\t[%.2f - %.2f] -> %.2f <- [%.2f - %.2f]\t',(fn / ws),(fn - (wb * fn)/2),fn,...
%					(fn + (wb * fn)/2),(min([0.95*(fs/2) fn * ws])));
%
%			filterStrings{fIdx} = sprintf('[%.2f - %.2f] -> %.2f <- [%.2f - %.2f]\t',(fn / ws),...
%					(fn - (wb * fn)/2),fn,(fn + (wb * fn)/2),(min([0.95*(fs/2) fn * ws])));
%	
%			passBandLp = (fn + (wb * fn)/2)/(fs/2);
%			passBandHp = (fn - (wb * fn)/2)/(fs/2);
%			stopBandLp = (min([0.95*(fs/2) fn * ws]))/(fs/2);
%			stopBandHp = (fn / ws)/(fs/2);
%
%			f(fIdx) = fn;
%			fn = fn * m;
%
%			filterSettings{fIdx,end} = {'bandpassfir',...
%					'StopbandFrequency1',stopBandHp,...
%					'PassbandFrequency1',passBandHp,...
%					'PassbandFrequency2',passBandLp,...
%					'StopbandFrequency2',stopBandLp,...
%					'StopbandAttenuation1',80,...
%					'StopbandAttenuation2',80,...
%					'DesignMethod','kaiserwin'};
%
%			
%			curHarmonics = find((harmonics > stopBandHp & harmonics < stopBandLp));
%			tmpHarmonics = harmonics.*(fs/2);
%
%			for harmIdx = curHarmonics
%					fprintf(' %f ', tmp(harmIdx));
%					filterSettings{fIdx,harmIdx} = {'bandstopfir',...
%							'FilterOrder',1500, ...
%		          'CutoffFrequency1',tmpHarmonics(harmIdx)-(str2num(options.bandStopWidth.Value)),...
%							'CutoffFrequency2',tmpHarmonics(harmIdx)+(str2num(options.bandStopWidth.Value)), ...
%							'SampleRate',fs};
%					
%			end
%			fprintf('\n');
%
%	end
%
%end
