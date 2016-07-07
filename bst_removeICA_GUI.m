function handles = bst_removeICA_GUI(EEG,iStudy)

WS   					= myssp('WinScale');				%-Window scaling factors
FS   					= myssp('FontSizes');			%-Scaled font sizes
PF   					= myssp_platform('fonts');			%-Font names (for this platform)
Rect 					= myssp('WinSize','Main','raw').*WS;		%-Main window
[SSPver,SSPc] = myssp('Ver','',1,1);

handles.Fmenu = figure('IntegerHandle','off',...
    'Name',sprintf('%s%s',' SSP ',myssp('ver'),' modality: EEG ',myssp('GetUser',' (%s)')),...
    'NumberTitle','off',...
    'Tag','Menu',...
    'Position',Rect,...
    'Resize','on',...
    'Color',[1 1 1]*.8,...
    'UserData',struct('SPMver',SSPver,'SPMc',SSPc),...
    'MenuBar','none',...
    'DefaultTextFontName',PF.helvetica,...
    'DefaultTextFontSize',FS(10),...
    'DefaultUicontrolFontName',PF.helvetica,...
    'DefaultUicontrolFontSize',FS(12),...
    'DefaultUicontrolInterruptible','on',...
    'Renderer','painters',...
    'Visible','on');

set(handles.Fmenu,'Units','normalized');

handles.Display.Panel=uipanel('Parent',handles.Fmenu,'Tag','PanelDataStructure','BackgroundColor',[1 1 1]*.8,...
    'Position',[0 0 0.35 1],'Units','normalized',...
    'Visible','on');

uicontrol(handles.Display.Panel,'Style','Text','String',['Display Panel'],...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1]*.8,...
    'FontName',PF.times,'FontSize',14,...
    'FontWeight','Bold',...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.01 0.96 0.98 0.03]);

handles.Display.AnalyFrame=uipanel('Parent',handles.Display.Panel,'BackgroundColor',[1 1 1]*.8,...
    'Position',[0.01 0.3 0.98 0.42],'Units','normalized',...
    'Visible','on');

set(handles.Display.AnalyFrame,'Title','Remove Component(s)','FOntSize',12,...
    'FontName',PF.times,'Fontweight','bold');

handles.Display.listboxWAR=uicontrol(handles.Display.Panel,'Style','listbox',...
    'String',[{' '}],...
    'Tag','listboxWAR','ForegroundColor','red','BackgroundColor',[1 1 1]*0.8,...
    'FontName',PF.times,'FontSize',10,'Enable','Inactive',...
    'HorizontalAlignment','left','Units','normalized',...
    'Callback','','Position',[0.01 0.05 0.98 0.2]);

handles.Display.DisplayPanel=uipanel('Parent',handles.Fmenu,...
    'Tag','DisplayDataStructure','BackgroundColor',[1 1 1]*.8,...
    'Position',[0.35 0 0.65 1],'Units','normalized',...
    'Title','','Visible','on','FOntWeight','bold');

handles.Display.CurAnalysis.ButtRmv=uicontrol(handles.Display.AnalyFrame,...
    'String',['>>'],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tooltipstring','Remove Component(s)',...
    'Tag','CurAnalysis',...
    'Callback','bst_removeica(''rmv'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.37 0.47 0.1 0.08]);

handles.Display.CurAnalysis.AllComponentstext=uicontrol(handles.Display.AnalyFrame,'Style','Text',...
    'String',['All Components'],...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1]*.8,...
    'Tag','CurAnalysis',...
    'FontName',PF.times,'FontSize',10,'Fontweight','bold',...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.15 0.65 0.2 0.1]);

handles.Display.CurAnalysis.AllComponents=uicontrol(handles.Display.AnalyFrame,'Style','listbox',...
    'String',[' '],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'Tag','CurAnalysis',...
    'FontName',PF.times,'FontSize',10,'Min',1,'Max',1000,...
    'Callback','bst_removeica(''pickcomponent'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.15 0.22 0.2 0.42]);

handles.Display.CurAnalysis.ButtCancel=uicontrol(handles.Display.AnalyFrame,...
    'String',['<<'],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'Tooltipstring','Keep Component(s)',...
    'Tag','CurAnalysis',...
    'FontName',PF.times,'FontSize',10,...
    'Callback','bst_removeica(''keep'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.37 0.37 0.1 0.08]);

handles.Display.CurAnalysis.SelComponentstext=uicontrol(handles.Display.AnalyFrame,'Style','Text',...
    'String',['Removed Components'],...
    'Tag','CurAnalysis',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1]*.8,...
    'FontName',PF.times,'FontSize',10,'Fontweight','bold',...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.5 0.65 0.2 0.1]);

handles.Display.CurAnalysis.SelComponents=uicontrol(handles.Display.AnalyFrame,'Style','listbox',...
    'String',[' '],'UserData',[],'Enable','off',...
    'Tag','CurAnalysis',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,'Min',1,'Max',1000,...
    'Callback','bst_removeica(''pickcomponent'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.5 0.22 0.2 0.42]);

handles.Display.CurAnalysis.ButtRmvandSave=uicontrol(handles.Display.AnalyFrame,...
    'String',['Remove and Save'],'UserData',iStudy,'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'Tooltipstring','Remove Selected Component(s) and Save',...
    'Tag','CurAnalysis',...
    'FontName',PF.times,'FontSize',10,...
    'Callback','bst_removeica(''removeandsave'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.15 0.05 0.3 0.1]);

handles.Display.ButtSave=uicontrol(handles.Display.AnalyFrame,...
    'String',['Save!'],'UserData',[],'Enable','on',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'Tooltipstring','Save changes',...
    'FontName',PF.times,'FontSize',10,...
    'Callback','bst_removeica(''save'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.55 0.05 0.3 0.1]);

handles.Display.CurAnalysis.PanelAllComponents=uipanel('Parent',handles.Display.DisplayPanel,'Tag','Panelcomponents','BackgroundColor',[1 1 1]*.8,...
    'Position',[0.01 0.88 0.98 0.075],'Units','normalized',...
    'Visible','on');

handles.Display.CurAnalysis.EEGLAB=uicontrol(handles.Display.DisplayPanel,...
    'String',['Components Activity (scroll)'],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tooltipstring','Display Raw Components (EEGLAB)',...
    'Tag','CurAnalysis',...
    'Callback','bst_removeica(''eeglab'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.25 0.91 0.2 0.04]);

handles.Display.CurAnalysis.EEGLAB3=uicontrol(handles.Display.DisplayPanel,...
    'String',['Spectral Bands'],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tooltipstring','Display Components Spectral Bands',...
    'Tag','CurAnalysis',...
    'Callback','bst_removeica(''eeglab3'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.5 0.91 0.2 0.04]);

handles.Display.CurAnalysis.EEGLAB2=uicontrol(handles.Display.DisplayPanel,...
    'String',['2D Scalp Maps (Spectra)'],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tooltipstring','Display Spectral 2D-Maps',...
    'Tag','CurAnalysis',...
    'Callback','bst_removeica(''eeglab2'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.75 0.91 0.2 0.04]);

handles.Display.CurAnalysis.AllCompTitle=uicontrol(handles.Display.DisplayPanel,'Style','Text',...
    'String',['All components:'],...
    'Tag','CurAnalysis',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1]*.8,...
    'FontName',PF.times,'FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.02 0.905 0.2 0.04]);

handles.Display.CurAnalysis.PanelSingleComponents=uipanel('Parent',handles.Display.DisplayPanel,'Tag','Panelsingcomponents','BackgroundColor',[1 1 1]*.8,...
    'Position',[0.01 0.01 0.98 0.86],'Units','normalized',...
    'Visible','on');

handles.Display.CurAnalysis.AxTopo=axes('Parent',handles.Display.DisplayPanel,'pos',[0.55 0.025 0.38 0.25],...
    'ydir','normal','box','on','Visible','off','YTick',[],'XTick',[],'Color',[1 1 1]*.8,'Tag','CurAnalysis');

handles.Display.CurAnalysis.AxChs2=axes('Parent',handles.Display.DisplayPanel,'pos',[0.1 0.60 0.4 0.2],...
    'ydir','normal','box','on','Visible','off','YTick',[],'XTick',[],'Color',[1 1 1]*.8,'Tag','CurAnalysis');

handles.Display.CurAnalysis.AxSpec=axes('Parent',handles.Display.DisplayPanel,'pos',[0.65 0.35 0.3 0.17],...
    'ydir','normal','box','on','Visible','off','YTick',[],'XTick',[],'Tag','CurAnalysis');

handles.Display.CurAnalysis.AxChs3=axes('Parent',handles.Display.DisplayPanel,'pos',[0.55 0.60 0.4 0.2],...
    'ydir','normal','box','on','Visible','off','YTick',[],'XTick',[],'Color',[1 1 1]*.8,'Tag','CurAnalysis');

handles.Display.CurAnalysis.AxSpec2=axes('Parent',handles.Display.DisplayPanel,'pos',[0.1 0.07 0.4 0.45],...
    'ydir','normal','box','on','Visible','off','Tag','CurAnalysis');

linkedaxes1(1) = handles.Display.CurAnalysis.AxChs2;
linkedaxes1(2) = handles.Display.CurAnalysis.AxChs3;

linkaxes(linkedaxes1,'xy')

handles.Display.CurAnalysis.ThrSlider=uicontrol(handles.Display.DisplayPanel,'Style','slider',...
    'String',[{'1'}],'Value',1,'Enable','off','Visible','off','Min',0,...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'Tooltipstring','Component in the display',...
    'FontName',PF.times,'FontSize',10,...
    'Tag','CurAnalysis',...
    'Callback','bst_removeica(''thrslider'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.05 0.60 0.015 0.2]);

handles.Display.CurAnalysis.ThrEdit=uicontrol(handles.Display.DisplayPanel,'Style','text',...
    'String',[{' '}],'Value',[],'Enable','off','Visible','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'Tooltipstring','Component in the display',...
    'FontName',PF.times,'FontSize',10,...
    'Tag','CurAnalysis',...
    'Callback','bst_removeica(''editthr'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.16 0.8 0.1 0.03]);

handles.Display.CurAnalysis.CompTitle=uicontrol(handles.Display.DisplayPanel,'Style','Text',...
    'String',['Components Properties:'],...
    'Tag','CurAnalysis',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1]*.8,...
    'FontName',PF.times,'FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.02 0.835 0.24 0.04]);

handles.Display.CurAnalysis.CompDisp=uicontrol(handles.Display.DisplayPanel,'Style','Text',...
    'String',['Component: '],...
    'Tag','CurAnalysis',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1]*.8,...
    'FontName',PF.times,'FontSize',10,...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.39 0.825 0.1 0.03]);

handles.Display.CurAnalysis.CompDispEdit=uicontrol(handles.Display.DisplayPanel,'Style','popupmenu',...
    'String',[{' '}],'Value',1,'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'Tooltipstring','Component in the display',...
    'FontName',PF.times,'FontSize',10,...
    'Tag','CurAnalysis',...
    'Callback','bst_removeica(''displaycomp'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','center','Units','normalized',...
    'Position',[0.49 0.83 0.1 0.03],...
    'Userdata',EEG);

handles.Display.CurAnalysis.NextComp=uicontrol(handles.Display.DisplayPanel,...
    'String',['>'],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tooltipstring','Display Next Component',...
    'Callback','bst_removeica(''upcomp'',get(gcf,''UserData''),gco);',...
    'Tag','CurAnalysis',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.595 0.845 0.03 0.015]);

handles.Display.CurAnalysis.PreviousComp=uicontrol(handles.Display.DisplayPanel,...
    'String',['<'],'UserData',[],'Enable','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tag','CurAnalysis',...
    'Tooltipstring','Display Previous Component',...
    'Callback','bst_removeica(''downcomp'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.595 0.83 0.03 0.015]);

handles.Display.CurAnalysis.ZoomIn=uicontrol(handles.Display.DisplayPanel,...
    'String',['+'],'UserData',[],'Enable','off','Visible','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tooltipstring','Zoom In',...
    'Callback','bst_removeica(''zoomin'',get(gcf,''UserData''),gco);',...
    'Tag','CurAnalysis',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.1 0.8 0.03 0.03]);

handles.Display.CurAnalysis.ZoomOut=uicontrol(handles.Display.DisplayPanel,...
    'String',['-'],'UserData',[],'Enable','off','Visible','off',...
    'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],...
    'FontName',PF.times,'FontSize',10,...
    'Tag','CurAnalysis',...
    'Tooltipstring','Zoom Out',...
    'Callback','bst_removeica(''zoomout'',get(gcf,''UserData''),gco);',...
    'HorizontalAlignment','left','Units','normalized',...
    'Position',[0.13 0.8 0.03 0.03]);

bst_removeica('start',handles)