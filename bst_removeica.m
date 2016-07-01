function [varargout] = bst_removeica(varargin)

modality=varargin{1};
handles=varargin{2};
varargout=[];

switch(lower(modality))
    
    case 'upcomp'
        v=get(handles.Display.CurAnalysis.CompDispEdit,'Value');
        str=get(handles.Display.CurAnalysis.CompDispEdit,'String');
        if v<size(str,1)
            set(handles.Display.CurAnalysis.CompDispEdit,'Value',v+1);
            set(handles.Display.CurAnalysis.NextComp,'enable','off');
            bst_removeica('resetaxes',handles);
            bst_removeica('displaycomp',handles);
            set(handles.Display.CurAnalysis.NextComp,'enable','on');
            set(handles.Display.CurAnalysis.AllComponents,'Value',v+1);
        end
        
    case 'downcomp'
        v=get(handles.Display.CurAnalysis.CompDispEdit,'Value');
        if v>1
            set(handles.Display.CurAnalysis.CompDispEdit,'Value',v-1);
            set(handles.Display.CurAnalysis.PreviousComp,'enable','off');
            bst_removeica('resetaxes',handles);
            bst_removeica('displaycomp',handles);
            set(handles.Display.CurAnalysis.PreviousComp,'enable','on');
            set(handles.Display.CurAnalysis.AllComponents,'Value',v-1);
        end
        
    case 'zoomin'
        if(~isfield(handles,'zoom'))
            h = zoom;
            handles.zoom = h;
        end
        
        h = handles.zoom;
        
        if(strcmp(get(h,'enable'),'on'))
            set(h,'enable','off')
        else
            set(h,'enable','on')
            set(h,'direction','in')
        end
        
    case 'zoomout'
        if(~isfield(handles,'zoom'))
            h = zoom;
            handles.zoom = h;
        end
        
        h = handles.zoom;
        
        if(strcmp(get(h,'enable'),'on'))
            set(h,'enable','off')
        else
            set(h,'enable','on')
            set(h,'direction','out')
        end
        
    case 'thrslider'
        
        v=get(handles.Display.CurAnalysis.CompDispEdit,'Value');
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        
        thr = (get(handles.Display.CurAnalysis.ThrSlider,'Value'));
        
        set(handles.Display.CurAnalysis.ThrEdit,'String',num2str(thr));
        set(handles.Display.CurAnalysis.ThrEdit,'Value',thr);
        
        EEG.wave_decomp.THR(v) = thr;
        
        set(handles.Display.CurAnalysis.CompDispEdit,'UserData',EEG);
        bst_removeica('plotcomps2',handles);
        
    case 'editthr'
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        thr = str2double(get(handles.Display.CurAnalysis.ThrEdit,'String'));
        set(handles.Display.CurAnalysis.ThrEdit,'Value',(thr));
        set(handles.Display.CurAnalysis.ThrSlider,'Value',(thr));
        
        EEG.wave_decomp.THR(selected_comp) = thr;
        
        set(handles.Display.CurAnalysis.CompDispEdit,'UserData',EEG);
        bst_removeica('plotcomps2',handles);
        
    case 'start'
        set(handles.Display.CurAnalysis.ButtRmv,'Enable','off','UserData',[]);
        set(handles.Display.CurAnalysis.ButtRmvandSave,'Enable','off');
        set(handles.Display.CurAnalysis.ButtCancel,'Enable','off','UserData',[]);
        set(handles.Display.CurAnalysis.AllComponents,'Enable','off','String',{' '},'Value',1,'UserData',[]);
        set(handles.Display.CurAnalysis.SelComponents,'Enable','off','String',{' '},'Value',1,'UserData',[]);
        handles=bst_removeica('resetfig',handles);
        EEG = get(handles.Display.CurAnalysis.CompDispEdit,'Userdata');
        
        components=1:size(EEG.icaweights,1);
        
        if(~isfield(EEG,'componentsremoved'))
            EEG.componentsremoved = [];
        end
        if(~isfield(EEG,'wave_decomp'))
            EEG.wave_decomp.THR = [];
        end

        set(handles.Display.CurAnalysis.AllComponents,'Enable','on','String',{' '},'Value',1,'UserData',components);
        set(handles.Display.CurAnalysis.ButtRmv,'Enable','on','UserData',[]);
        set(handles.Display.CurAnalysis.ButtRmvandSave,'Enable','on')
        set(handles.Display.CurAnalysis.ButtCancel,'Enable','on','UserData',[]);
        
        set(handles.Display.CurAnalysis.SelComponents,'Enable','on','String',{' '},'Value',1,'UserData',[]);        
        set(handles.Display.CurAnalysis.CompDispEdit,'Enable','on','String',{' ' },'Value',1)
        
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data;
        
        [fft_sig, f] = mypwelch(EEG.icaact,512,512,EEG.srate);
        
        EEG.spectrum.f = f;
        EEG.spectrum.data = fft_sig;
        
        [SWA,SWD, ~, NSAMPLES] = wave_decomp_ICA(EEG.icaact,8,'db4');
        EEG.wave_decomp.SWD = SWD;
        EEG.wave_decomp.SWA = SWA;
        if(isempty(EEG.wave_decomp.THR))
            EEG.wave_decomp.THR = squeeze(max(max(abs(SWD)))+1);
        end
        EEG.wave_decomp.NSAMPLES = NSAMPLES;
        
        set(handles.Display.CurAnalysis.CompDispEdit,'UserData',EEG);        
        set(handles.Display.CurAnalysis.SelComponents,'UserData',EEG.componentsremoved);

        bst_removeica('resetaxes',handles);
        bst_removeica('updatelistboxes',handles);
        bst_removeica('displaycomp',handles);
        set(handles.Display.CurAnalysis.NextComp,'Enable','on','UserData',[]);
        set(handles.Display.CurAnalysis.EEGLAB,'Enable','on','UserData',[]);
        set(handles.Display.CurAnalysis.EEGLAB2,'Enable','on','UserData',[]);
        set(handles.Display.CurAnalysis.EEGLAB3,'Enable','on','UserData',[]);
        set(handles.Display.CurAnalysis.PreviousComp,'Enable','on','UserData',[]);
        set(handles.Display.CurAnalysis.ZoomIn,'Enable','on','Visible','on','UserData',[]);
        set(handles.Display.CurAnalysis.ZoomOut,'Enable','on','Visible','on','UserData',[]);
        set(handles.Display.CurAnalysis.ThrSlider,'Enable','on','Visible','on','UserData',[]);
        set(handles.Display.CurAnalysis.ThrEdit,'Enable','on','Visible','on','UserData',[]);
        
    case 'resetfig'
        cla(handles.Display.CurAnalysis.AxTopo,'reset')
        cla(handles.Display.CurAnalysis.AxSpec,'reset')
        cla(handles.Display.CurAnalysis.AxChs3,'reset')
        cla(handles.Display.CurAnalysis.AxSpec2,'reset')
        cla(handles.Display.CurAnalysis.AxChs2,'reset')
        set(handles.Display.CurAnalysis.AxChs2,'Visible','off')
        set(handles.Display.CurAnalysis.AxTopo,'Visible','off')
        set(handles.Display.CurAnalysis.AxSpec,'Visible','off')
        set(handles.Display.CurAnalysis.AxChs3,'Visible','off')
        set(handles.Display.CurAnalysis.AxSpec2,'Visible','off')
        set(handles.Display.CurAnalysis.CompDispEdit,'Enable','off','String',{' '},'Value',1)
        set(handles.Display.CurAnalysis.NextComp,'Enable','off','UserData',[])
        set(handles.Display.CurAnalysis.PreviousComp,'Enable','off','UserData',[])
        set(handles.Display.CurAnalysis.ZoomIn,'Enable','off','Visible','off','UserData',[])
        set(handles.Display.CurAnalysis.ZoomOut,'Enable','off','Visible','off','UserData',[])
        set(handles.Display.CurAnalysis.ThrSlider,'Enable','off','Visible','off','UserData',[])
        set(handles.Display.CurAnalysis.ThrEdit,'Enable','off','Visible','off','UserData',[])
        
        set(handles.Display.CurAnalysis.EEGLAB,'Enable','off','UserData',[])
        set(handles.Display.CurAnalysis.EEGLAB2,'Enable','off','UserData',[])
        set(handles.Display.CurAnalysis.EEGLAB3,'Enable','off','UserData',[])
        varargout{1}=handles;
        set(handles.Fmenu,'UserData',handles); %added 13/10/2008
        
    case 'resetaxes'
        cla(handles.Display.CurAnalysis.AxTopo,'reset')
        cla(handles.Display.CurAnalysis.AxSpec,'reset')
        cla(handles.Display.CurAnalysis.AxChs3,'reset')
        cla(handles.Display.CurAnalysis.AxSpec2,'reset')
        cla(handles.Display.CurAnalysis.AxChs2,'reset')
        
        set(handles.Display.CurAnalysis.AxChs2,'Visible','off')
        set(handles.Display.CurAnalysis.AxTopo,'Visible','off')
        set(handles.Display.CurAnalysis.AxSpec,'Visible','off')
        set(handles.Display.CurAnalysis.AxChs3,'Visible','off')
        set(handles.Display.CurAnalysis.AxSpec2,'Visible','off')
        
    case 'eeglab'
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        eegplot( EEG.icaact, 'srate', EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
            'limits', [EEG.xmin EEG.xmax]*1000 , 'command', []);
        
    case 'eeglab2'
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        mypop_topoplot( EEG, 0, 1:1:size(EEG.icawinv,2), '2D Scalp Maps (Spectra)' ,[], 'flat', 3);
        
    case 'eeglab3'
        set(handles.Fmenu,'Pointer','watch')
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        try
            hf=bst_removeica_plotspectralbands(EEG);
            set(handles.Fmenu,'Pointer','arrow')
            set(hf,'Visible','on')
        catch
            set(handles.Fmenu,'Pointer','arrow')
            close(hf)
        end
        
    case 'rmv'
        v=get(handles.Display.CurAnalysis.AllComponents,'Value');
        allcomps=get(handles.Display.CurAnalysis.AllComponents,'UserData');
        selchns=get(handles.Display.CurAnalysis.SelComponents,'UserData');
        selchns=sort(union(selchns,allcomps(v)));
        set(handles.Display.CurAnalysis.SelComponents,'UserData',selchns);
        bst_removeica('updatelistboxes',handles);
        
    case 'keep'
        v=get(handles.Display.CurAnalysis.SelComponents,'Value');
        selchns=get(handles.Display.CurAnalysis.SelComponents,'UserData');
        if ~isempty(selchns)
            selchns(v)=[];
            set(handles.Display.CurAnalysis.SelComponents,'UserData',selchns);
            
            bst_removeica('updatelistboxes',handles);
        end
        
    case 'updatelistboxes'
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        
        components=get(handles.Display.CurAnalysis.SelComponents,'UserData');
        
        [ varegg ] = myssp_eeglab_computevariance( EEG.data,...
            { EEG.icasphere EEG.icaweights }, EEG.icawinv,...
            setdiff(1:size(EEG.icaweights,1), components));
        
        resetmsg(handles);
        updatemsg(num2str(varegg),handles);
        
        selchns=get(handles.Display.CurAnalysis.SelComponents,'UserData');
        allcomps=get(handles.Display.CurAnalysis.AllComponents,'UserData');
        strcomp={' '};
        for i=1:1:length(allcomps)
            strcomp{i}=['N ',num2str(allcomps(i))];
        end
        strcomp2={' '};
        for i=1:1:length(selchns)
            strcomp2{i}=['N ',num2str(selchns(i))];
        end
        set(handles.Display.CurAnalysis.SelComponents,'String',strcomp2,'Value',1)
        set(handles.Display.CurAnalysis.AllComponents,'String',strcomp)
        set(handles.Display.CurAnalysis.CompDispEdit,'String',strcomp)
        
    case 'plotcomps2'
        v=get(handles.Display.CurAnalysis.CompDispEdit,'Value');
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        swd=EEG.wave_decomp.SWD(:,:,v);
        swa=EEG.wave_decomp.SWA(:,:,v);
        thr=EEG.wave_decomp.THR(v);
        nsamples=EEG.wave_decomp.NSAMPLES;
        
        artefact = compute_artefact(swd,swa,thr,nsamples,'db4');
        new_icaact = EEG.icaact(v,:) - artefact;
        [new_spectrum, freqs] = mypwelch([new_icaact; artefact],512,512,EEG.srate);
        
        if(isempty(get(handles.Display.CurAnalysis.AxChs2,'userdata')))
            axes(handles.Display.CurAnalysis.AxChs2)
            hold on
            h_artefact_plot=plot(handles.Display.CurAnalysis.AxChs2,EEG.Time,artefact,'r');
            set(gca,'userdata',h_artefact_plot);
            hold off
        else
            set(get(handles.Display.CurAnalysis.AxChs2,'userdata'),'YData',artefact);
        end
        
        if(isempty(get(handles.Display.CurAnalysis.AxSpec,'userdata')))
            axes(handles.Display.CurAnalysis.AxSpec)
            hold on
            h_spectrum_plot(1)=plot(freqs,10*log10(new_spectrum(2,:)),'-r');
            h_spectrum_plot(2)=plot(freqs,10*log10(new_spectrum(1,:)),'-g');
            set(h_spectrum_plot,'LineWidth',2)
            
            set(gca,'userdata',h_spectrum_plot);
            hold off
        else
            h_spectrum_plot = get(handles.Display.CurAnalysis.AxSpec,'userdata');
            set(h_spectrum_plot(1),'YData',10*log10(new_spectrum(2,:)));
            set(h_spectrum_plot(2),'YData',10*log10(new_spectrum(1,:)));
            
        end
        
        if(isempty(get(handles.Display.CurAnalysis.AxChs3,'userdata')))
            axes(handles.Display.CurAnalysis.AxChs3)
            
            h_new_icaact_plot=plot(EEG.Time,new_icaact,'g');
            set(gca, 'xlim', [EEG.Time(1) EEG.Time(end)],'userdata',h_new_icaact_plot);
            set(gca,'ytick',[])
            xlabel('Time(s)')
            title(['Original - Artefactual Component N',num2str(v)])
        else
            
            set(get(handles.Display.CurAnalysis.AxChs3,'userdata'),'YData',new_icaact);
        end
        
        set(handles.Display.CurAnalysis.AxChs3,'Visible','on')
        
        
    case 'displaycomp'
        set(handles.Fmenu,'Pointer','watch')
        try
            v=get(handles.Display.CurAnalysis.CompDispEdit,'Value');
            set(handles.Display.CurAnalysis.AllComponents,'Value',v);
            
            EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
            
            data=EEG.wave_decomp.SWD(:,:,v);
            
            set(handles.Display.CurAnalysis.ThrSlider,'Max',max(max(abs(data)))+1);
            
            thr = EEG.wave_decomp.THR(v);
            
            set(handles.Display.CurAnalysis.ThrSlider,'Value',thr);
            set(handles.Display.CurAnalysis.ThrEdit,'Value',thr);
            set(handles.Display.CurAnalysis.ThrEdit,'String',num2str(thr));
            
            spectra=EEG.spectrum.data(v,:);
            freqs = EEG.spectrum.f;
            
            axes(handles.Display.CurAnalysis.AxSpec)
            he=plot(freqs,10*log10(spectra),'-b');
            set(he,'LineWidth',2)
            set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)');
            xlabel('Frequency (Hz)')
            title('Power Spectrum')
            xlim([0 40])
            
            data=EEG.icaact(v,:);
            axes(handles.Display.CurAnalysis.AxChs2)
            plot(EEG.Time,data);
            set(gca, 'xlim', [EEG.Time(1) EEG.Time(end)]);
            ylabel('\mu V')
            
            axes(handles.Display.CurAnalysis.AxSpec2)
            pnts = length(EEG.times);
            trials = length(EEG.Time)/pnts;
            imagesc(EEG.times,1:trials,reshape(data,[pnts trials])');
            set(gca,'clim',[-max(abs(min((data))), (abs(max((data))))) ...
                max(abs(min((data))), (abs(max((data)))))]);
            
            axis xy
            xlabel('Time (s)');
            ylabel('trial')
            
            axes(handles.Display.CurAnalysis.AxTopo);
            
            topoplot(EEG.icawinv(:,v), EEG.chanlocs,'style','map');
            
%             bst_topography_topoplot(EEG.icawinv(:,v), EEG.chanlocs,'flat', 3,8); %axis square;
            
            set(handles.Fmenu,'Pointer','arrow')
            set(handles.Display.CurAnalysis.AxTopo,'Visible','off')
            set(handles.Display.CurAnalysis.AxSpec,'Visible','on')
            set(handles.Display.CurAnalysis.AxSpec2,'Visible','on')
            set(handles.Display.CurAnalysis.AxChs2,'Visible','on')
            
        catch
            bst_removeica('resetfig',handles);
            set(handles.Fmenu,'Pointer','arrow')
        end
        
    case 'removeandsave'
        
        updatemsg(' recomposing signals',handles);
        EEG				= get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        components= get(handles.Display.CurAnalysis.SelComponents,'UserData');
        sInputs = get(handles.Display.CurAnalysis.ButtRmvandSave,'UserData');
        
        SWD = EEG.wave_decomp.SWD;
        SWA = EEG.wave_decomp.SWA;
        NSAMPLES = EEG.wave_decomp.NSAMPLES;
        THR = EEG.wave_decomp.THR;
        
        [ data, ~] = myssp_eeglab_compvar( EEG.icaact,...
            { EEG.icasphere EEG.icaweights }, EEG.icawinv,...
            SWD,SWA,NSAMPLES,THR, setdiff(1:size(EEG.icaweights,1), components));
                
        EEG.componentsremoved=components;
        EEG.icaact = [];
        EEG.wave_decomp = rmfield(EEG.wave_decomp,{'SWD','SWA','NSAMPLES'});
        EEG.Comment     = ['ICA: ' datestr(clock)];
        
        ProtocolInfo = bst_get('ProtocolInfo');
        [p,f] = bst_fileparts(sInputs.FileName);
        fPath = strrep(p, ProtocolInfo.STUDIES, '');
		% Create a default output filename 
		OutputFiles = fullfile(ProtocolInfo.STUDIES,fPath,f);		% Save on disk
		bst_save(OutputFiles, EEG,'v6');
		% Register in database
		db_add_data(sInputs.iStudy, OutputFiles, EEG);
        
        DataMatOut          = db_template('datamat');
		DataMatOut.DataType = 'data';
        DataMatOut.Comment  = ['ICA_pruned: ' datestr(clock)];
        DataMatOut.Time     = EEG.Time;
        DataMatOut.F        = zeros(length(EEG.ChannelFlag),length(EEG.Time));
        DataMatOut.F(EEG.iChannels,:) = data;
        DataMatOut.ChannelFlag = EEG.ChannelFlag;
		% Create a default output filename 
		OutputFiles = fullfile(ProtocolInfo.STUDIES,fPath,'ICA_pruned_data.mat');		% Save on disk
        
		% Save on disk
		bst_save(OutputFiles, DataMatOut,'v6');
		% Register in database
		db_add_data(sInputs.iStudy, OutputFiles, DataMatOut);

        set(handles.Fmenu,'UserData',handles);
        set(handles.Display.CurAnalysis.ButtRmvandSave,'UserData',true);
        varargout{1}=handles;

        updatemsg(' done ',handles);
        db_reload_database();
        
    case 'save'
        updatemsg(' saving ICA',handles);
        EEG=get(handles.Display.CurAnalysis.CompDispEdit,'UserData');
        components=get(handles.Display.CurAnalysis.SelComponents,'UserData');
        sInputs = get(handles.Display.CurAnalysis.ButtRmvandSave,'UserData');
                
        EEG.componentsremoved=components;
        EEG.Comment     = ['ICA: ' datestr(clock)];
        EEG.icaact = [];

        ProtocolInfo = bst_get('ProtocolInfo');
        [p,f] = bst_fileparts(sInputs.FileName);
        fPath = strrep(p, ProtocolInfo.STUDIES, '');
		% Create a default output filename 
		OutputFiles = fullfile(ProtocolInfo.STUDIES,fPath,f);
		% Save on disk
		bst_save(OutputFiles, EEG,'v6');
		% Register in database
		db_add_data(sInputs.iStudy, OutputFiles, EEG);
        updatemsg(' done',handles);
                
        db_reload_database();

        set(handles.Fmenu,'UserData',handles);
        varargout{1}=handles;
        
end

function []=updatemsg(curmsg,handles)
msg=get(handles.Display.listboxWAR,'String');
msg{end+1}=['   ... ',curmsg,' ...'];
set(handles.Display.listboxWAR,'String',msg);
set(handles.Display.listboxWAR,'Value',length(msg));
set(handles.Display.listboxWAR,'UserData',curmsg);

function []=resetmsg(handles)
msg{1}='';
set(handles.Display.listboxWAR,'String',msg);
set(handles.Display.listboxWAR,'Value',length(msg));
