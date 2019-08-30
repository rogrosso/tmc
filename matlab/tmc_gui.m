function tmc_gui
    % SIMPLE_GUI2 Select a data set from the pop-up menu, then
    % click one of the plot-type push buttons. Clicking the button
    % plots the selected data in the axes.

    close all;

    % keep track of current state values
    cstate_ = containers.Map;
    cstate_('StandardMC') = 0;
    cstate_('Isovalue') = 0.0293;
    cstate_('Levelset') = 1;
    cstate_('CellTri')  = 0;
    cstate_('BlueLines') = 0;
    cstate_('InnerHexagon') = 0;
    cstate_('Contours')  = 0;
    cstate_('VertexConfig') = 0;
    cstate_('MCcase')   = 105; % this is MC case 13
    cstate_('MCcaseName') = 'case13';
    cstate_('PredefCase') = 14; % no predefined configuration
    cstate_('LevelsetHandle') = 0;
    cstate_('LinesHandles') = [];
    cstate_('InnerHexaHandles') = [];
    cstate_('ContourHandles') = [];
    cstate_('VertexConfigHandles') = [];
    cstate_('TrisHandles') = [];
    cstate_('MCTrisHandles') = [];
    cstate_('SymmTraHandle') = 0;
    cstate_('PositiveVertHandle') = 0;
    cstate_('SaveFig') = 0;
    cstate_('Recompute') = 0;
    cstate_('ColorId') = 1;
    cstate_('FunctionValues') = [];
    cstate_('AsymptoticDecider') = 0;
    cstate_('Asymptotes') = [];
    cstate_('Video') = 0;
    cstate_('Animation') = 0;
    cstate_('HyperbolicArcsU') = [];
    cstate_('HyperbolicArcsV') = [];
    cstate_('HyperbolicArcsW') = [];
    cstate_('PlotHyperbolicArcsW') = 0;

    %  Create and then hide the GUI as it is being constructed.
    w_uicont = 250;
    h_uicont = 700;
    f = figure('Visible','off','Position',[10,10,w_uicont,h_uicont]);
    cstate_('FigureHandle') = f;
    hold on;
    %  Construct the components.
    p_uicont = 50;
    p_vpos    = 675;
    p_vcheckb = 30;
    p_vpopup  = 32;
    p_vtext1  = 35;
    p_vtext2  = 42;
    p_vpushb  = 32;
    hStandardMC = uicontrol('Parent',f,'Style','checkbox','String','Standard MC',...
        'Position',[p_uicont,p_vpos,100,25],...
        'Callback',@StandardMCcheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hLevelset = uicontrol('Parent',f,'Style','checkbox','String','Level set',...
        'Position',[p_uicont,p_vpos,70,25],...
        'Callback',@levelsetcheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hCellTri = uicontrol('Parent',f,'Style','checkbox','String','Cell triangulation',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@celltricheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hBlueLines = uicontrol('Parent',f,'Style','checkbox','String','Blue lines',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@bluelinescheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hInnerHexagon = uicontrol('Parent',f,'Style','checkbox','String','Inner hexagon',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@innerhexagoncheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hContours = uicontrol('Parent',f,'Style','checkbox','String','Contours',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@contourscheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hVertexConfig = uicontrol('Parent',f,'Style','checkbox','String','Vertex +/-',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@vertexconfigcheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hAsymptoticDecider = uicontrol('Parent',f,'Style','checkbox','String','Asymptotes',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@asymptoticdecidercheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hHyperbolicArcs = uicontrol('Parent',f, 'Style','checkbox','String', 'Hyperbolic Arcs',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@hyperbolicarcscheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hVideo = uicontrol('Parent',f,'Style','checkbox','String','Video',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@videocheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hAnimation = uicontrol('Parent',f,'Style','checkbox','String','Animation',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@animationcheckbox_Callback);

    p_vpos = p_vpos - p_vtext1;
    hMCcase_text = uicontrol('Parent',f,'Style','text','String','Select MC case',...
        'Position',[p_uicont,p_vpos,130,25],'HorizontalAlignment','Left');

    p_vpos = p_vpos - p_vpopup;
    hMCcase_popup = uicontrol('Parent',f,'Style','popupmenu',...
        'String',{'case 0','case 1','case 2','case 3','case 4','case 5','case 6',...
        'case 7','case 8','case 9','case 10','case 11','case 12','case 13','case 14'},...
        'Position',[p_uicont,p_vpos,100,25],...
        'Callback',@mccasepopupmenu_Callback);

    p_vpos = p_vpos - p_vtext2;
    hSymm_text = uicontrol('Parent',f,'Style','text','String','Apply a randomly selected symmetry transform?', ...
        'Position',[p_uicont,p_vpos,130,30],'HorizontalAlignment','Left');

    p_vpos = p_vpos - p_vcheckb;
    hSymmTra_checkbox = uicontrol('Parent',f,'Style','checkbox','String','Symmetry transformation',...
        'Position',[p_uicont,p_vpos,200,30],...
        'Callback',@symmtracheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hPositiveVert_checkbox = uicontrol('Parent',f,'Style','checkbox','String','Positive vertices',...
        'Position',[p_uicont,p_vpos,200,30],...
        'Callback',@positivevertcheckbox_Callback);

    p_vpos = p_vpos - p_vtext2;
    hPredefCase_text = uicontrol('Parent',f,'Style','text','String','Choose a predifined configuration', ...
        'Position',[p_uicont,p_vpos,130,25],'HorizontalAlignment','Left');

    p_vpos = p_vpos -p_vtext1;
    hPredefCase_popup = uicontrol('Parent',f,'Style','popupmenu','String',...
        {'none','contour with 6 vts','contour with 7 vts, 1st conf.','contour with 7 vts, 2nd conf.', ...
        'contour with 8 vts, type A','contour with 8 vts, type B','contour with 9 vts','contour with 9 vts',...
        'MC 4 with tunnel','MC 7 with tunnel','MC 10 with tunnel','MC 12 with tunnel',...
        'MC 13 with tunnel','MC 13, contour 12 vts','MC 13, contour 12 vts','MC 13 simple','Pathological case'},...
        'Position',[p_uicont,p_vpos,150,30],...
        'Callback',@predefpopup_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hSave_checkbox = uicontrol('Parent',f,'Style','checkbox','String','Save fig file?',...
        'Position',[p_uicont,p_vpos,200,30],...
        'Callback',@savecheckbox_Callback);


    % recompute with same config
    p_vpos = p_vpos - p_vpushb;
    hRecompute_pushbutton = uicontrol('Parent',f,'Style','pushbutton','String','Recompute',...
        'Position',[p_uicont,p_vpos,80,20],'Callback',@recomppushbutton_Callback);

    % help window
    p_vpos = p_vpos - p_vpushb;
    hHelpDialog_pushbutton = uicontrol('Parent',f,'Style','pushbutton','String','Help dialog',...
        'Position',[p_uicont,p_vpos,100,20],'Callback',@helppushbutton_Callback);

    ltUIc = [hStandardMC,...
             hLevelset,...
             hCellTri,...
             hBlueLines,...
             hInnerHexagon,...
             hContours,...
             hVertexConfig,...
             hAsymptoticDecider,...
             hHyperbolicArcs,...
             hVideo,...
             hMCcase_text,...
             hMCcase_popup,...
             hSymm_text,...
             hSymmTra_checkbox,...
             hPredefCase_text,...
             hPredefCase_popup,...
             hRecompute_pushbutton,...
             hHelpDialog_pushbutton];
    align(ltUIc,'Left','none');

    % Initialize the GUI.
    % Change units to normalized so components resize
    % automatically.
    f.Units = 'normalized';
    %ha.Units = 'normalized';
    hStandardMC.Units = 'normalized';
    hLevelset.Units = 'normalized';
    hCellTri.Units = 'normalized';
    hBlueLines.Units = 'normalized';
    hContours.Units  = 'normalized';
    hVertexConfig.Units = 'normalized';
    hAsymptoticDecider.Units = 'normalized';
    hHyperbolicArcs.Units = 'normalized';
    hVideo.Units = 'normalized';
    hMCcase_text.Units = 'normalized';
    hMCcase_popup.Units = 'normalized';
    hSymm_text.Units = 'normalized';
    hSymmTra_checkbox.Units = 'normalized';
    hPositiveVert_checkbox.Units = 'normalized';
    hPredefCase_text.Units = 'normalized';
    hPredefCase_popup.Units = 'normalized';
    hRecompute_pushbutton.Units = 'normalized';
    hHelpDialog_pushbutton.Units = 'normalized';

    % set defualt values
    hStandardMC.set('Value',0);
    hLevelset.set('Value',1);
    hCellTri.set('Value',0);
    hMCcase_popup.set('Value',14);

    % Finish initialization
    axis off;
    set(f,'menubar','none');
    % Assign the GUI a name to appear in the window title.
    f.Name = 'tcm GUI';
    % Move the GUI to the center of the screen.
    movegui(f,'northeast')
    % Make the GUI visible.
    f.Visible = 'on';

    % draw unit cell
    f3d = figure('Position',[10,10,700,700],'Color',[1,1,1]);
    f3d.Name = 'tcm';
    set(0,'currentfigure',f3d);
    cstate_('CurrentFig') = f3d;
    movegui(f3d,'center');
    axes( 'Position',[.22 .21 .58 .58],'CameraViewAngleMode','manual');

    hold on;
    vm = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    h_bbox = patch('Vertices',vm,'Faces',fm,'FaceColor','none');
    daspect([1 1 1]);
    view(25,50);
    axis tight;
    axis square;
    %axis off;
    camlight;
    lighting gouraud;
    xlabel('u');
    ylabel('v');
    zlabel('w');
    %set(gca,'xtick',[]);
    %set(gca,'ytick',[]);
    %set(gca,'ztick',[]);
    %set(gca,'Position',[0.05 0.05 0.9 0.9]);


    % draw isosurface for default case 13
    f = zeros(1,8);
    t_mc(cstate_);
    hPredefCase_popup.set('Value',1);
    cstate_('PredefCase') = 1;


    %  Callbacks for simple_gui. These callbacks automatically
    %  have access to component handles and initialized data
    %  because they are nested at a lower level.

    %  Pop-up menu callback. Read the pop-up menu Value property
    %  to determine which item is currently displayed and make it
    %  the current data.
    function mccasepopupmenu_Callback(source,eventdata)
        % Determine the selected data set.
        % Set current data to the selected data set.
        hPredefCase_popup.set('Value',1);
        cstate_('PredefCase') = 1;
        cstate_('Isovalue') = 0;
        switch source.String{source.Value}
            case 'case 0'
                cstate_('MCcase') =  0;
                cstate_('MCcaseName') = 'case00';
            case 'case 1'
                cstate_('MCcase') =  1;
                cstate_('MCcaseName') = 'case01';
            case 'case 2'
                cstate_('MCcase') =  3;
                cstate_('MCcaseName') = 'case02';
            case 'case 3'
                cstate_('MCcase') =  6;
                cstate_('MCcaseName') = 'case03';
            case 'case 4'
                cstate_('MCcase') =  24;
                cstate_('MCcaseName') = 'case04';
            case 'case 5'
                cstate_('MCcase') =  7;
                cstate_('MCcaseName') = 'case05';
            case 'case 6'
                cstate_('MCcase') =  25;
                cstate_('MCcaseName') = 'case06';
            case 'case 7'
                cstate_('MCcase') =  22;
                cstate_('MCcaseName') = 'case07';
            case 'case 8'
                cstate_('MCcase') =  15;
                cstate_('MCcaseName') = 'case08';
            case 'case 9'
                cstate_('MCcase') =  23;
                cstate_('MCcaseName') = 'case09';
            case 'case 10'
                cstate_('MCcase') =  60;
                cstate_('MCcaseName') = 'case10';
            case 'case 11'
                cstate_('MCcase') =  27;
                cstate_('MCcaseName') = 'case11';
            case 'case 12'
                cstate_('MCcase') =  30;
                cstate_('MCcaseName') = 'case12';
            case 'case 13'
                cstate_('MCcase') =  105;
                cstate_('MCcaseName') = 'case13';
            case 'case 14'
                cstate_('MCcase') =  29;
                cstate_('MCcaseName') = 'case14';
            otherwise
                cstate_('MCcase') =  105;
                cstate_('MCcaseName') = 'case13';
        end
        cstate_('Recompute') = 1;
        t_mc(cstate_);
    end

    % Push button callbacks. Each callback plots current_data in
    % the specified plot type.

    function StandardMCcheckbox_Callback(source,eventdata)
        val = hStandardMC.get('Value');
        if val == 1
            cstate_('StandardMC') = 1;
        else
            cstate_('StandardMC') = 0;
        end
        t_mc(cstate_);
    end

    function levelsetcheckbox_Callback(source,eventdata)
        % Display surf plot of the currently selected data.
        val = hLevelset.get('Value');
        if val == 1
            %hLevelset.set('Value',1);
            cstate_('Levelset') = 1;
        else
            %hLevelset.set('Value',0);
            cstate_('Levelset') = 0;
        end
            t_mc(cstate_);
    end

    function celltricheckbox_Callback(source,eventdata)
        % Display mesh plot of the currently selected data.
        val = hCellTri.get('Value');
        if val == 1
            cstate_('CellTri')  = 1;
        else
            cstate_('CellTri')  = 0;
        end
        t_mc(cstate_);
    end

    function bluelinescheckbox_Callback(source,eventdata)
        val = hBlueLines.get('Value');
        if val == 1
            cstate_('BlueLines') = 1;
        else
            cstate_('BlueLines') = 0;
        end
        t_mc(cstate_);
    end

    function innerhexagoncheckbox_Callback(source,eventdata)
        val = hInnerHexagon.get('Value');
        if val == 1
            cstate_('InnerHexagon') = 1;
        else
            cstate_('InnerHexagon') = 0;
        end
        t_mc(cstate_);
    end

    function contourscheckbox_Callback(source,eventdata)
        val = hContours.get('Value');
        if val == 1
            cstate_('Contours') = 1;
        else
            cstate_('Contours') = 0;
        end
        t_mc(cstate_);
    end

    function vertexconfigcheckbox_Callback(source,eventdata)
        val = hVertexConfig.get('Value');
        if val == 1
            cstate_('VertexConfig') = 1;
        else
            cstate_('VertexConfig') = 0;
        end
        t_mc(cstate_);
    end

    function asymptoticdecidercheckbox_Callback(source,eventdata)
        val = hAsymptoticDecider.get('Value');
        if val == 1
            cstate_('AsymptoticDecider') = 1;
        else
            cstate_('AsymptoticDecider') = 0;
        end
        t_mc(cstate_);
    end

    function hyperbolicarcscheckbox_Callback(source,eventdata)
        val = hHyperbolicArcs.get('Value');
        if val==1
            cstate_('PlotHyperbolicArcsW') = 1;
        else
            cstate_('PlotHyperbolicArcsW') = 0;
        end
        t_mc(cstate_);
    end

    function videocheckbox_Callback(source,eventdata)
        val = hVideo.get('Value');
        if val == 1
            cstate_('Video') = 1;
        else
            cstate_('Video') = 0;
        end
        t_mc(cstate_);
    end

    function animationcheckbox_Callback(source,eventdata)
        val = hAnimation.get('Value');
        if val == 1
            cstate_('Animation') = 1;
        else
            cstate_('Animation') = 0;
        end
        t_mc(cstate_);
    end

    function predefpopup_Callback(source,eventdata)
        % select one of the predefined cases
        val = hPredefCase_popup.get('Value');
        switch val
            case 1 % no config. selected
                cstate_('PredefCase') = 1;
                cstate_('MCcase') = 105;
                cstate_('Isovalue') = 0;
            case 2 % contour with 6 vertices, MC 3, bit pattern 6
                cstate_('PredefCase') = 2;
                cstate_('MCcase') = 6;
                cstate_('Isovalue') = 0.0387;
                hMCcase_popup.set('Value',4);
            case 3 % contour with 7 vertices, Type A with three straight lines, MC 6, bit pattern 25
                cstate_('PredefCase') = 3;
                cstate_('MCcase') = 25;
                cstate_('Isovalue') = 9.8588;
                hMCcase_popup.set('Value',7);
            case 4 % contour with 7 vertices, Type B with two straight lines, MC 6, bit pattern 25
                cstate_('PredefCase') = 4;
                cstate_('MCcase') = 25;
                cstate_('Isovalue') = 9.9994608191478135;
                hMCcase_popup.set('Value',7);
            case 5 % contour with 8 vertices, case a, MC 12, bit pattern 30
                cstate_('PredefCase') = 5;
                cstate_('MCcase') = 30;
                cstate_('Isovalue') = 0.0097;
                hMCcase_popup.set('Value',13);
            case 6 % contour with 8 vertices, case b, MC 10, bit pattern 60
                cstate_('PredefCase') = 6;
                cstate_('MCcase') = 60;
                cstate_('Isovalue') = 9.99946;
                hMCcase_popup.set('Value',11);
            case 7 % contour with 9 vertices, MC 7, bit pattern 22
                cstate_('PredefCase') = 7;
                cstate_('MCcase') = 22;
                cstate_('Isovalue') = -1.8061;
                hMCcase_popup.set('Value',8);
            case 8 % contour with 9 vertices, MC 7, bit pattern 22
                cstate_('PredefCase') = 8
                cstate_('MCcase') = 22;
                cstate_('Isovalue') = 1233.5999999999999;
                hMCcase_popup.set('Value',8);
            case 9 % MC 4 with tunnel, bit pattern 24
                cstate_('PredefCase') = 9;
                cstate_('MCcase') = 24;
                cstate_('Isovalue') = -1.7660;
                hMCcase_popup.set('Value',5);
            case 10 % MC 7 with tunnel, bit pattern 22
                cstate_('PredefCase') = 10;
                cstate_('MCcase') = 22;
                cstate_('Isovalue') = -0.6;
                hMCcase_popup.set('Value',8);
            case 11 % MC 10 with tunnel, bit pattern 60
                cstate_('PredefCase') = 11;
                cstate_('MCcase') = 60;
                cstate_('Isovalue') = 1.10918;
                hMCcase_popup.set('Value',11);
            case 12 % MC 12 with tunnel, bit pattern 30
                cstate_('PredefCase') = 12;
                cstate_('MCcase') = 30;
                cstate_('Isovalue') = 0.0708;
                hMCcase_popup.set('Value',13);
            case 13 % MC 13 with tunnel, MC 3, bit pattern 105
                cstate_('PredefCase') = 13;
                cstate_('MCcase') = 105;
                cstate_('Isovalue') = -1.3064;
                hMCcase_popup.set('Value',14);
            case 14 % MC 13 with contour with 12 vertices, bit pattern 105
                cstate_('PredefCase') = 14;
                cstate_('MCcase') = 105;
                cstate_('Isovalue') = 0.0293;
                hMCcase_popup.set('Value',14);
            case 15 % MC 13 with contour with 12 vertices, bit pattern 105
                cstate_('PredefCase') = 15;
                cstate_('MCcase') = 105;
                cstate_('Isovalue') = 1007.4;
                hMCcase_popup.set('Value',14);
            case 16 % MC 13 without tunnel and without contour with 12 vertices
                cstate_('PredefCase') = 16;
                cstate_('MCcase') = 105;
                cstate_('Isovalue') = 0.0293;
                hMCcase_popup.set('Value',14);
            case 17 % MC  pathologic case
                cstate_('PredefCase') = 17;
                cstate_('MCcase') = 199;
                cstate_('Isovalue') = 1000;
                hMCcase_popup.set('Value',7);
        end
        t_mc(cstate_);
    end


    function savecheckbox_Callback(source,eventdata)
        val = hSave_checkbox.get('Value');
        if val == 1
            cstate_('SaveFig') = 1;
        else
            cstate_('SaveFig') = 0;
        end
        t_mc(cstate_);
    end

    function recomppushbutton_Callback(source,eventdata)
        cstate_('Recompute') = 1;
        t_mc(cstate_);
    end

    function symmtracheckbox_Callback(source,eventdata)
        val = hSymmTra_checkbox.get('Value');
        if val == 1
            cstate_('SymmTraHandle') = 1;
            t_mc(cstate_);
        elseif val == 0
            cstate_('SymmTraHandle') = 0;
        end
        %t_mc(cstate_);
    end

    function positivevertcheckbox_Callback(source,eventdata)
        val = hPositiveVert_checkbox.get('Value');
        if val == 1
            cstate_('PositiveVertHandle') = 1;
        elseif val == 0
            cstate_('PositiveVertHandle') = 0;
        end
        t_mc(cstate_);
    end

    function helppushbutton_Callback(source,eventdata)
        m_text = 'how to use this stuff';
        message = {m_text,m_text,m_text,m_text,m_text};
        Title = 'Help window';
        msgbox(message,'title');
    end

end
