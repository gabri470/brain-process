function varargout = process_gaitIniPlotERSD( varargin )

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Gait plot';
sProcess.FileTag     = '';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Pre-process';
sProcess.Index       = 69;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'timefreq','data', 'results', 'matrix'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Default values for some options
sProcess.processDim  = 1;    % Process channel by channel
sProcess.isSeparator = 1;

% Definition of the options
SelectOptions = {...
    '', ...                               % Filename
    '', ...                               % FileFormat
    'open', ...                           % Dialog type: {open,save}
    'Select directory...', ...            % Window title
    'ImportData', ...                     % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
    'single', ...                         % Selection mode: {single,multiple}
    'dirs', ...                           % Selection mode: {files,dirs,files_and_dirs}
    bst_get('FileFilters', 'events'), ... % Get all the available file formats
    'EventsIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn
% Option: Event file
sProcess.options.evtfile.Comment = 'Saving directory:';
sProcess.options.evtfile.Type    = 'filename';
sProcess.options.evtfile.Value   = SelectOptions;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function sOutput = Run(sProcess, sInput) %#ok<DEFNU>
sOutput = {};
for iFile = 1:length(sInput)
    
    
%    channel = in_bst_channel(sInput(iFile).ChannelFile);
    
    data = in_bst_timefreq(sInput(iFile).FileName);
    
    freq = data.Freqs;
      
%    SavingDir  = sProcess.options.evtfile.Value{1};
%    
%    if(~exist(fullfile(SavingDir,sInput(iFile).SubjectName),'dir'))
%        mkdir(fullfile(SavingDir,sInput(iFile).SubjectName))
%    end
%    
%    SavingDir = fullfile(SavingDir,sInput(iFile).SubjectName);
%    Filename = strjoin({sInput(iFile).SubjectName, sInput(iFile).Condition},'-');


		betaBand = freq>= 13 & freq <= 30;
    
    hf=figure();
    set(hf,'PaperType','a4');
    set(hf,'PaperPosition',[0 0 get(hf,'PaperSize')]);
    set(hf,'PaperUnits','centimeters');
    set(hf,'Units','centimeters');
    set(hf,'Position',[0 0 0.75*get(hf,'PaperSize')]);
    set(hf,'Units','pixels');

		plot(freq,squeeze(10.*log10(data.TF)));
		title(Filename);
		legend(data.RowNames);
		xlabel('Frequency (Hz) ');
		ylabel('Log(PSD)');
		xlim([2 80]);
    
    filename = fullfile(SavingDir,strjoin({Filename,'1.pdf'},'-'));
    print(hf,'-dpdf',filename);
    close(hf)
end
end

function [curfig,com] = timetopoplot( data,chanlocs, time_latencies, topotitle, rowcols, curfig)

com = '';
if nargin < 1
    help pop_topoplot;
    return;
end;
if isempty(chanlocs)
    disp('Error: cannot plot topography without channel location file'); return;
end;


% additional options
% ------------------
options    = { 'masksurf' 'on' };
outoptions = { options{:} }; % for command

nbgraph = size(time_latencies(:),1);
if ~exist('topotitle')
    topotitle = '';
end;
if ~exist('rowcols') | isempty(rowcols) | rowcols == 0
    paperSize = get(curfig,'PaperSize');
    sx = floor(sqrt(nbgraph*paperSize(1)/paperSize(2)));
    sy = ceil(nbgraph/sx);
    rowcols(2) = sx;
    rowcols(1) = sy;
end;

SIZEBOX = 150;

fprintf('Plotting...\n');
if isempty(chanlocs)
    fprintf('Error: set has no channel location file\n');
    return;
end;

% plot the graphs
% ---------------
counter = 1;
countobj = 1;
allobj = zeros(1,1000);
if(isempty(curfig))
    curfig = figure;
end

clim = ([min(data(:)) max(data(:))]);

for index = 1:length(time_latencies)
    if nbgraph > 1
        if mod(index, rowcols(1)*rowcols(2)) == 1
%             if index> 1, figure(curfig); a = textsc(0.5, 0.05, topotitle); set(a, 'fontweight', 'bold'); end;
            
%             pos = get(curfig,'Position');
%             posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
%             posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
%             set(curfig,'Position', [posx posy  SIZEBOXx*rowcols(2)  SIZEBOXy*rowcols(1)]);
            try, icadefs; set(curfig, 'color', BACKCOLOR); catch, end;
        end;
        curax = subplot( rowcols(1), rowcols(2), mod(index-1, rowcols(1)*rowcols(2))+1);
        set(curax, 'visible', 'on','xtick',[],'ytick',[])
        
    end;
    
    % plot scalp map
    % --------------
    %fprintf('Printing to figure %d.\n',curfig);
    if ~isnan(time_latencies(index))
        
        if nbgraph > 1, axes(curax); end;
        tmpobj =topoplot(data(:,index), chanlocs,'style','map');
        set(gca,'clim',clim);
        
        %         tmpobj = myssp_eeg_view_topography_topoplot( data(:,index), chanlocs, shading, numcontour,4,clim);
        if nbgraph == 1,
            figure(curfig); if nbgraph > 1, axes(curax); end;
            title( [ 'Latency ' num2str(time_latencies(index)) ' ms from ' topotitle] );
        else
            figure(curfig); if nbgraph > 1, axes(curax); end;
            title([num2str(time_latencies(index)) ' ms']);
        end;
        
        allobj(countobj:countobj+length(tmpobj)-1) = tmpobj;
        countobj = countobj+length(tmpobj);
        drawnow;
        axis equal;
        if index == length(time_latencies)
            if nbgraph == 1
                clim = get(gca, 'clim');
                pos = get(gca,'position');
                q = [pos(1) pos(2) 0 0];
                s = [pos(3) pos(4) pos(3) pos(4)];
                col = colormap;
                ax = axes('position', [0.95 0 .05 1].*s+q);
                cbar(ax,[1:64],clim);
            else
                cbar('vert');
            end;
        end;
    else
        axis off
    end;
end;
if nbgraph> 1,
    figure(curfig); a = textsc(0.5, 0.05, topotitle);
    set(a, 'fontweight', 'bold');
end;
if nbgraph== 1,
    com = 'figure;';
end;
set(allobj(1:countobj-1), 'visible', 'on');

figure(curfig);
axcopy(curfig, 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position''''); set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); clear postmp;');

com = [com sprintf('pop_topoplot(%s,%d, %s);', ...
    inputname(1), vararg2str({time_latencies topotitle rowcols outoptions{:} }))];
return;
end
