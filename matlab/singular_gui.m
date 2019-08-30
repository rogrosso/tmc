function singular_gui
    % SIMPLE_GUI2 Select a data set from the pop-up menu, then
    % click one of the plot-type push buttons. Clicking the button
    % plots the selected data in the axes.

    close all;

    % keep track of current state values
    cstate_ = containers.Map;
    cstate_('StandardMC') = 0;
    cstate_('Isovalue') = 1233.6;
    cstate_('Levelset') = 1;
    cstate_('CellTri')  = 0;
    cstate_('BlueLines') = 0;
    cstate_('InnerHexagon') = 0;
    cstate_('Contours')  = 0;
    cstate_('VertexConfig') = 0;
    cstate_('MCcase')   = 105; % this is MC case 13
    cstate_('MCcaseName') = 'case13';
    cstate_('PredefCase') = 14; % no predefined configuration
    cstate_('LevelsetHandle') = [];
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
    cstate_('HyperbolicArcsU') = [];
    cstate_('HyperbolicArcsV') = [];
    cstate_('HyperbolicArcsW') = [];
    cstate_('BBox') = [];
    cstate_('View') = [25,50];        

    %  Create and then hide the GUI as it is being constructed.
    w_uicont = 250;
    h_uicont = 500;
    fi_ = figure('Visible','off','Position',[10,10,w_uicont,h_uicont]);
    cstate_('FigureHandle') = fi_;
    hold on;
    %  Construct the components.
    p_uicont = 50;
    p_vpos    = 470;
    p_vcheckb = 30;
    p_vtext1  = 35;
    p_vtext2  = 42;
    p_vpushb  = 32;

    p_vpos = p_vpos - p_vcheckb;
    hLevelset = uicontrol('Parent',fi_,'Style','checkbox','String','Level set',...
        'Position',[p_uicont,p_vpos,70,25],...
        'Callback',@levelsetcheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hCellTri = uicontrol('Parent',fi_,'Style','checkbox','String','Cell triangulation',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@celltricheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hBlueLines = uicontrol('Parent',fi_,'Style','checkbox','String','Blue lines',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@bluelinescheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hContours = uicontrol('Parent',fi_,'Style','checkbox','String','Contours',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@contourscheckbox_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hAsymptoticDecider = uicontrol('Parent',fi_,'Style','checkbox','String','Asymptotes',...
        'Position',[p_uicont,p_vpos,120,25],...
        'Callback',@asymptoticdecidercheckbox_Callback);
    
    p_vpos = p_vpos - p_vtext2;
    hPredefCase_text = uicontrol('Parent',fi_,'Style','text','String','Choose a predifined configuration', ...
        'Position',[p_uicont,p_vpos,130,25],'HorizontalAlignment','Left');

    p_vpos = p_vpos -p_vtext1;
    hPredefCase_popup = uicontrol('Parent',fi_,'Style','popupmenu','String',...
        {'Pair one','Pair two','Pair three', 'Pair four','cruel case','via dolorosa'},...
        'Position',[p_uicont,p_vpos,150,30],...
        'Callback',@predefpopup_Callback);

    p_vpos = p_vpos - p_vcheckb;
    hSave_checkbox = uicontrol('Parent',fi_,'Style','checkbox','String','Save fig file?',...
        'Position',[p_uicont,p_vpos,200,30],...
        'Callback',@savecheckbox_Callback);

    % help window
    p_vpos = p_vpos - p_vpushb;
    hHelpDialog_pushbutton = uicontrol('Parent',fi_,'Style','pushbutton','String','Help dialog',...
        'Position',[p_uicont,p_vpos,100,20],'Callback',@helppushbutton_Callback);

    ltUIc = [hLevelset,...
             hCellTri,...
             hBlueLines,...
             hContours,...
             hAsymptoticDecider,...
             hPredefCase_popup,...
             hHelpDialog_pushbutton];
    align(ltUIc,'Left','none');

    % Initialize the GUI.
    % Change units to normalized so components resize
    % automatically.
    fi_.Units = 'normalized';
    hLevelset.Units = 'normalized';
    hCellTri.Units = 'normalized';
    hBlueLines.Units = 'normalized';
    hContours.Units  = 'normalized';
    %hVertexConfig.Units = 'normalized';
    hAsymptoticDecider.Units = 'normalized';
    hPredefCase_popup.Units = 'normalized';
    hHelpDialog_pushbutton.Units = 'normalized';

    % set defualt values
    hLevelset.set('Value',1);
    hCellTri.set('Value',0);

    % Finish initialization
    axis off;
    set(fi_,'menubar','none');
    % Assign the GUI a name to appear in the window title.
    fi_.Name = 'tcm GUI';
    % Move the GUI to the center of the screen.
    movegui(fi_,'northeast')
    % Make the GUI visible.
    fi_.Visible = 'on';


    % draw isosurface for default case 13
    % draw unit cell
    f3d = figure('Position',[10,10,900,900],'Color',[1,1,1]);
    f3d.Name = 'tcm';
    set(0,'currentfigure',f3d);
    cstate_('CurrentFig') = f3d;
    movegui(f3d,'center');
    axes( 'Position',[.22 .21 .58 .58],'CameraViewAngleMode','manual');

    hold on;
    daspect([1 1 1]);
    view_ = cstate_('View');
    view(view_(1),view_(2));
    axis equal tight;
    %axis square; %normal or square
    %axis off;
    camlight;
    lighting gouraud;
    %xlabel('u');
    %ylabel('v');
    %zlabel('w');
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);
    %set(gca,'fontsize',14);
    %set(gca,'Position',[0.05 0.05 0.9 0.9]);
    
    % init plot
    hPredefCase_popup.set('Value',1);
    cstate_('PredefCase') = 1;
    t_singular(cstate_);


    %  Callbacks for simple_gui. These callbacks automatically
    %  have access to component handles and initialized data
    %  because they are nested at a lower level.

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
            t_singular(cstate_);
    end

    function celltricheckbox_Callback(source,eventdata)
        % Display mesh plot of the currently selected data.
        val = hCellTri.get('Value');
        if val == 1
            cstate_('CellTri')  = 1;
        else
            cstate_('CellTri')  = 0;
        end
        t_singular(cstate_);
    end

    function bluelinescheckbox_Callback(source,eventdata)
        val = hBlueLines.get('Value');
        if val == 1
            cstate_('BlueLines') = 1;
        else
            cstate_('BlueLines') = 0;
        end
        t_singular(cstate_);
    end

    function contourscheckbox_Callback(source,eventdata)
        val = hContours.get('Value');
        if val == 1
            cstate_('Contours') = 1;
        else
            cstate_('Contours') = 0;
        end
        t_singular(cstate_);
    end

    function asymptoticdecidercheckbox_Callback(source,eventdata)
        val = hAsymptoticDecider.get('Value');
        if val == 1
            cstate_('AsymptoticDecider') = 1;
        else
            cstate_('AsymptoticDecider') = 0;
        end
        t_singular(cstate_);
    end

    function predefpopup_Callback(source,eventdata)
        % select one of the predefined cases
        val = hPredefCase_popup.get('Value');
        cstate_('Isovalue') = 1233.6;
        switch val
            case 1 % no config. selected
                cstate_('PredefCase') = 1;
            case 2 % contour with 6 vertices, MC 3, bit pattern 6
                cstate_('PredefCase') = 2;
            case 3 % contour with 7 vertices, Type A with three straight lines, MC 6, bit pattern 25
                cstate_('PredefCase') = 3;
            case 4 % contour with 7 vertices, Type B with two straight lines, MC 6, bit pattern 25
                cstate_('PredefCase') = 4;
            case 5 
                cstate_('PredefCase') = 5; % cruel case
            case 6
                cstate_('PredefCase') = 6; % via dolorosa, 
        end
        t_singular(cstate_);
    end


    function savecheckbox_Callback(source,eventdata)
        val = hSave_checkbox.get('Value');
        if val == 1
            cstate_('SaveFig') = 1;
        else
            cstate_('SaveFig') = 0;
        end
        t_singular(cstate_);
    end

    function helppushbutton_Callback(source,eventdata)
        m_text = 'how to use this stuff';
        message = {m_text,m_text,m_text,m_text,m_text};
        Title = 'Help window';
        msgbox(message,'title');
    end

end
