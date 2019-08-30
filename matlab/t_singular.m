function t_singular(cstate_)
    % use for triangle mesh [128,128,128];
    % read current state values
    i0 = cstate_('Isovalue');
    levelset_ = cstate_('Levelset');
    celltris_ = cstate_('CellTri');
    bluelines_ = cstate_('BlueLines');
    innerhexa_ = cstate_('InnerHexagon');
    contours_  = cstate_('Contours');
    vertexconfig_ = cstate_('VertexConfig');
    standardMC_ = cstate_('StandardMC');
    mc = cstate_('MCcase');
    symm_ = cstate_('SymmTraHandle');
    positiveVert_ = cstate_('PositiveVertHandle');
    c_pos = cstate_('ColorId');
    asympt_ = cstate_('AsymptoticDecider');

    % check if reset values of f
    recomputeF = cstate_('Recompute');
    cstate_('Recompute') = 0;
    
    
    
    f3d =  cstate_('CurrentFig');
    set(0,'currentfigure',f3d);
    f = cstate_('FunctionValues');
    f1 = [];
    fsht = 0;
    %
    if recomputeF
        c_pos = randi(14);
        cstate_('ColorId') = c_pos;
    end
    face_shift = [0, 0, 0];
    switch cstate_('PredefCase')
        case 1
            fsht = 1;
            face_shift = [1, 0, 0];
            f = [899.00000000000000,1339.0000000000000,1210.0000000000000,1220.0000000000000,754.00000000000000,998.00000000000000, 928.00000000000000, 1264.0000000000000];
            f1 = [1339.0000000000000,2543.0000000000000,1220.0000000000000,2040.0000000000000,998.00000000000000,2443.0000000000000,1264.0000000000000,2421.0000000000000];
        case 2
            fsht = 2;
            face_shift = [0, 1, 0];
            f = [1003.0000000000000,1523.0000000000000,1230.0000000000000,1246.0000000000000,1062.0000000000000,1448.0000000000000,1293.0000000000000,1029.0000000000000];
            f1 = [1230.0000000000000,1246.0000000000000, 1133.0000000000000, 978.00000000000000, 1293.0000000000000, 1029.0000000000000, 1228.0000000000000,673.00000000000000];
        case 3
            fsht = 2;
            face_shift = [0,1,0];
            f = [1165.0000000000000, 1073.0000000000000, 1311.0000000000000, 1225.0000000000000, 1046.0000000000000, 1018.0000000000000, 1203.0000000000000, 1237.0000000000000];
            f1 = [1311.0000000000000,1225.0000000000000,1255.0000000000000,1256.0000000000000,1203.0000000000000,1237.0000000000000,1118.0000000000000,1288.0000000000000];
        case 4
            fsht = 2;
            face_shift = [0, 1, 0];
            f = [1763.0000000000000,1052.0000000000000,1815.0000000000000,1050.0000000000000,960.00000000000000,1325.0000000000000,1150.0000000000000,1260.0000000000000];
            f1 = [1815.000000000000,1050.0000000000000,2134.0000000000000,1106.0000000000000,1150.0000000000000,1260.0000000000000,1248.0000000000000,1144.0000000000000];
        case 5
            fsht = 2;
            face_shift = [0,1,0];
            F = [2.74742516087490,-3.39187542578189,-12.5297639669456,0.431517989649243,-6.92460546400188,2.52228314017858,14.6950568276448,-10.0732624062474];
            i0 = -1.3064;
            [f,f1] = newCellVertices(F,i0);
        case 6 % via dolorosa
            fsht = 2;
            face_shift = [0,1,0];
            F = [-15.6504952739285,2.90290077342601,24.5454566157887,-24.5274127623786,21.6741877710053,-4.49696327433901,-19.7891575872492,-15.5588482753161];
            i0 = -1.8061;
            [f, f1] = newCellVertices2(F,i0);
    end

    % keep current cell configuration
    cstate_('FunctionValues') = f;

    % apply symmetry transform?
    if symm_ == 1
        t_index = randi(24);
        fres = transform(t_index,f,positiveVert_);
        f = fres;
        if positiveVert_ == 1
            i0 = -i0;
        end
    end

    
    % set colors
    color = getColor(c_pos);
    c_pos = 2;
    color1 = getColor(c_pos);
    color2 = getColor(c_pos + 1);
   
    % plot bounding boxes
    clearBBox(cstate_);
    vm = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    h_bbox = patch('Vertices',vm,'Faces',fm,'FaceColor','none');
    bbox_ = [];
    bbox_ = [bbox_, h_bbox];
    for i = 1:8
        vm(i,:) = vm(i,:) + face_shift;
    end
    h_bbox = patch('Vertices',vm,'Faces',fm,'FaceColor','none');
    bbox_ = [bbox_, h_bbox];
    cstate_('BBox') = bbox_;
    
    % material 
    ka = 0.6;
    kd = 0.8;
    ks = 0.7;
    sh = 40;
    material([ka kd ks sh]);
        
    % generate isosurface
    if levelset_ == 1
        cc =0:0.01:1;
        sz = length(cc);
        [x,y,z] = meshgrid(1:sz,1:sz,1:sz);
        fd =zeros(sz,sz,sz);
        for ii=1:sz
            for jj=1:sz
                for kk=1:sz
                    u = cc(ii);
                    v = cc(jj);
                    w = cc(kk);
                    x(ii,jj,kk) = cc(ii);
                    y(ii,jj,kk) = cc(jj);
                    z(ii,jj,kk) = cc(kk);
                    fd(ii,jj,kk) = trilinear(f,u,v,w);
                end
            end
        end

        clearLevelSetPlot(cstate_);
        lt_phandle = [];
        p_handle = patch(isosurface(x,y,z,fd,i0),'FaceColor',color1,'EdgeColor','none');
        % material 
        ka = 0.6;
        kd = 0.8;
        ks = 0.7;
        sh = 40;
        material([ka kd ks sh]);
        lt_phandle = [lt_phandle, p_handle];
        if fsht == 1
                x = x + 1;
        elseif fsht == 2
                y = y + 1;
        elseif fsht == 3
            z = z + 1;
        end
        for ii=1:sz
            for jj=1:sz
                for kk=1:sz
                    u = cc(ii);
                    v = cc(jj);
                    w = cc(kk);
                    fd(ii,jj,kk) = trilinear(f1,u,v,w);
                end
            end
        end
        p_handle = patch(isosurface(x,y,z,fd,i0),'FaceColor',color2,'EdgeColor','none');
        lt_phandle = [lt_phandle, p_handle];
        cstate_('LevelsetHandle') = lt_phandle;
        ka = 0.6;
        kd = 0.8;
        ks = 0.7;
        sh = 40;
        material([ka kd ks sh]);
    else
        clearLevelSetPlot(cstate_);
    end
    
    % check if blue lines have to be plotted
    if bluelines_
        clearBlueLines(cstate_);
        plotStraightLines(cstate_,f,f1,i0,fsht);
    else
        clearBlueLines(cstate_);
    end

    % check for cell triangulation
    if celltris_
        clearAllTriangles(cstate_);
        if standardMC_ == 1
            s_mc(cstate_,fsht,f,f1,mc,i0);
        else
            triangulateCell(cstate_,f,f1,fsht,mc,i0)
        end
    else
        clearAllTriangles(cstate_);
    end

    % plot contours
    if contours_
        clearContours(cstate_);
        computeCellContours(cstate_,fsht,f,f1,mc,i0);
    else
        clearContours(cstate_);
    end

    % plot vertex sign used to clasify MC
    if vertexconfig_
        clearVertexConfig(cstate_);
        plotVertexConfig(cstate_,fsht,f,f1,mc,i0);
    else
        clearVertexConfig(cstate_);
    end

    % plot asymptotes
    if asympt_
        clearAsymptotes(cstate_);
        plotAsymptotes(cstate_,fsht,f,f1,mc,i0);
    else
        clearAsymptotes(cstate_);
    end

    if innerhexa_
        clearInnerHexagon(cstate_);
        plotInnerHexagon(cstate_,f,mc,i0);
    else
        clearInnerHexagon(cstate_);
    end
end


function setLinesVisibility(cstate_,vis_flag)
    l_handles = cstate_('LinesHandles');
    for l = 1:length(l_handles)
        if vis_flag
            set(l_handles(l),'Visible','on');
        else
            set(l_handles(l),'Visible','off');
        end
    end
end

function setInnerHexVisibility(cstate_,vis_flag)
    l_handles = cstate_('InnerHexaHandles');
    for l = 1:length(l_handles)
        if vis_flag
            set(l_handles(l),'Visible','on');
        else
            set(l_handles(l),'Visible','off');
        end

        set(l_handles(l),'Color','y');
        %set(l_handles(l),'EdgeLighting','none');
        set(l_handles(l),'LineWidth',4);
    end
end

function setTrisVisibility(cstate_,vis_flag)
    l_handles = cstate_('TrisHandles');
    for l = 1:length(l_handles)
        if vis_flag
            set(l_handles(l),'Visible','on');
        else
            set(l_handles(l),'Visible','off');
        end
    end
end

% Trilinear interpolation
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = trilinear(f,u,v,w)
    val = (1-w)*((1-v)*(f(1)*(1-u)+f(2)*u)+v*(f(3)*(1-u)+f(4)*u))+ w*((1-v)*(f(5)*(1-u)+f(6)*u)+v*(f(7)*(1-u)+f(8)*u));
end

% compute new cell vertices to produce extreme case
function [f1,f2] = newCellVertices(f,i0)
    [ic,cvt,fs,ui,vi,wi] = quadraticEquations(f,i0);
    
    % obtain 8 vertices
    small = 0.00000000000001;
    f1 = zeros(1,8);
    f2 = zeros(1,8);
    v1 = cvt(1,2);
    v2 = cvt(4,2);
    u1 = 0;
    u2 = 1;
    w1 = 0;
    w2 = 1;
    f1(1) = trilinear(f,u1,v1,w1);
    f1(2) = trilinear(f,u2,v1,w1);
    f1(3) = trilinear(f,u1,v2,w1);
    f1(4) = trilinear(f,u2,v2,w1);
    f1(5) = trilinear(f,u1,v1,w2);
    f1(6) = trilinear(f,u2,v1,w2);
    f1(7) = trilinear(f,u1,v2,w2);
    f1(8) = trilinear(f,u2,v2,w2);
    
    [ic,cvt,fs,ui,vi,wi] = quadraticEquations(f1,i0);
    
    f2(1) = f1(3);
    f2(2) = f1(4);
    f2(5) = f1(7);
    f2(6) = f1(8);
    f2(3) = f2(1) + 1;
    f2(4) = f2(1) + 1.5;
    f2(7) = f2(1) + 0.1;
    f2(8) = f2(1) + 1.5;
end

% compute new cell vertices to produce extreme case
function [f2,f1] = newCellVertices2(f,i0)
    [ic,cvt,fs,ui,vi,wi] = quadraticEquations(f,i0);
    
    % obtain 8 vertices
    small = 0.0000000000000005;
    f1 = zeros(1,8);
    f2 = zeros(1,8);
    v1 = cvt(4,2) - small;
    v2 = 1;
    u1 = cvt(6,1) - small;
    u2 = 1;
    w1 = 0;
    w2 = 1;
    f1(1) = trilinear(f,u1,v1,w1);
    f1(2) = trilinear(f,u2,v1,w1);
    f1(3) = trilinear(f,u1,v2,w1);
    f1(4) = trilinear(f,u2,v2,w1);
    f1(5) = trilinear(f,u1,v1,w2);
    f1(6) = trilinear(f,u2,v1,w2);
    f1(7) = trilinear(f,u1,v2,w2);
    f1(8) = trilinear(f,u2,v2,w2);
    
    [ic,cvt,fs,ui,vi,wi] = quadraticEquations(f1,i0);
    
    f2(3) = f1(1);
    f2(4) = f1(2);
    f2(7) = f1(5);
    f2(8) = f1(6);
    f2(1) = f2(3) + 1;
    f2(2) = f2(3) + 1.5;
    f2(5) = f2(3) + 0.1;
    f2(6) = f2(3) + 1.5;
end


%
% Compute the bilinear coefficients for a pair of opposite faces
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f1, f2, f3, f4, h1, h2, h3, h4] = setBilinearform(F,face)
    switch(face)
        case 1 % face w = 0
            f1 = F(1); f2 = F(2); f3 = F(3); f4 = F(4);
            h1 = F(5); h2 = F(6); h3 = F(7); h4 = F(8);
        case 2 % face w = 1
            f1 = F(5); f2 = F(6); f3 = F(7); f4 = F(8);
            h1 = F(1); h2 = F(2); h3 = F(3); h4 = F(4);
        case 3 % face v = 0
            f1 = F(1); f2 = F(2); f3 = F(5); f4 = F(6);
            h1 = F(3); h2 = F(4); h3 = F(7); h4 = F(8);
        case 4 % face v = 1
            f1 = F(3); f2 = F(4); f3 = F(7); f4 = F(8);
            h1 = F(1); h2 = F(2); h3 = F(5); h4 = F(6);
        case 5 % face u = 0
            f1 = F(1); f2 = F(3); f3 = F(5); f4 = F(7);
            h1 = F(2); h2 = F(4); h3 = F(6); h4 = F(8);
        case 6 % face u = 1
            f1 = F(2); f2 = F(4); f3 = F(6); f4 = F(8);
            h1 = F(1); h2 = F(3); h3 = F(5); h4 = F(7);
    end % end switch

end % set bilinearform


%
% Plot vertex configuration, i.e. positive and negative
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotVertexConfig(cstate_,coord,f,f1,mc,i0)
    vt = [[0,0,0];[1,0,0];[0,1,0];[1,1,0];[0,0,1];[1,0,1];[0,1,1];[1,1,1]];
    v_handles = [];
    for v = 1:8
        if f(v) > i0
            v_ = plot3(vt(v,1),vt(v,2),vt(v,3),'m.','MarkerSize',40);
            v_handles = [v_handles,v_];
        end
    end
    xsh = 0;
    ysh = 0;
    zsh = 0;
    if coord == 1
        xsh = 1;
    elseif coord == 2
        ysh = 1;
    elseif coord == 3
        zsh = 1;
    end
    for v = 1:8
        if f1(v) > i0
            v_ = plot3(vt(v,1)+xsh,vt(v,2)+ysh,vt(v,3)+zsh,'r.','MarkerSize',40);
            v_handles = [v_handles,v_];
        end
    end
    cstate_('VertexConfigHandles') = v_handles;
end
% Plot straight lines lying on the level set joining
% opposite faces
function plotStraightLines(cstate_,F,F1,i0,coord)
    [ic,cvt,fs,ui,vi,wi] = quadraticEquations(F,i0);
    % check if there is a tunnel
    ict = ic(1) + ic(2) + ic(3);

    % plot straight lines joining opposite faces
    lt_handles = [];
    for i = 1:2
        if (0 < ui(i) && ui(i) < 1) && (0 < vi(i) && vi(i) < 1) 
            cu1 = [ui(i),ui(i)];
            cv1 = [vi(i),vi(i)];
            cw1 = [0, 1];
            l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
            lt_handles = [lt_handles,l_];
        end
    end
    for i = 1:2
        if (0 < ui(i) && ui(i) < 1) && (0 < wi(i) && wi(i) < 1)
            cu1 = [ui(i),ui(i)];
            cw1 = [wi(i),wi(i)];
            cv1 = [0, 1];
            l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
            lt_handles = [lt_handles,l_];
        end
    end
    for i = 1:2
        wk = zeros(1,2);
        wk(1) = wi(2);
        wk(2) = wi(1);
        if (0 < vi(i) && vi(i) < 1) && (0 < wk(i) && wk(i) < 1)
            cv1 = [vi(i),vi(i)];
            cw1 = [wk(i),wk(i)];
            cu1 = [0, 1];
            l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
            lt_handles = [lt_handles,l_];
        end
    end
    
    [ic,cvt,fs,ui,vi,wi] = quadraticEquations(F1,i0);
    % check if there is a tunnel
    ict = ic(1) + ic(2) + ic(3);

    % plot straight lines joining opposite faces
    xsh = [0,0];
    ysh = [0,0];
    zsh = [0,0];
    if coord == 1
        xsh = [1,1];
    elseif coord == 2
        ysh = [1,1];
    elseif coord == 3
        zsh = [1,1];
    end
        
    for i = 1:2
        if (0 <= ui(i) && ui(i) <= 1) && (0 <= vi(i) && vi(i) <= 1) 
            cu1 = [ui(i),ui(i)] + xsh;
            cv1 = [vi(i),vi(i)] + ysh;
            cw1 = [0, 1] + zsh;
            l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
            lt_handles = [lt_handles,l_];
        end
    end
    for i = 1:2
        if (0 <= ui(i) && ui(i) <= 1) && (0 <= wi(i) && wi(i) <= 1)
            cu1 = [ui(i),ui(i)] + xsh;
            cw1 = [wi(i),wi(i)] + zsh;
            cv1 = [0, 1] + ysh;
            l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
            lt_handles = [lt_handles,l_];
        end
    end
    for i = 1:2
        wk = zeros(1,2);
        wk(1) = wi(2);
        wk(2) = wi(1);
        if (0 <= vi(i) && vi(i) <= 1) && (0 <= wk(i) && wk(i) <= 1)
            cv1 = [vi(i),vi(i)] + ysh;
            cw1 = [wk(i),wk(i)] + zsh;
            cu1 = [0, 1] + xsh;
            l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
            lt_handles = [lt_handles,l_];
        end
    end
    
    cstate_('LinesHandles') = lt_handles;
           
end % plotStraightLines(cstate_,F,F1,i0)


% Delete graphic objects from figure
%
function clearLevelSetPlot(cstate_)
    lt_handles = cstate_('LevelsetHandle');
	for l =1:length(lt_handles)
        delete(lt_handles(l));
    end
    lt_handles = [];
    cstate_('LevelsetHandle');
end

function clearBlueLines(cstate_)
    lt_handles = cstate_('LinesHandles');
    for l =1:length(lt_handles)
        delete(lt_handles(l));
    end
    lt_handles = [];
    cstate_('LinesHandles') = lt_handles;
end

function clearInnerHexagon(cstate_)
    lt_handles = cstate_('InnerHexaHandles');
    for l =1:length(lt_handles)
        delete(lt_handles(l));
    end
    lt_handles = [];
    cstate_('InnerHexaHandles') = lt_handles;
end

function clearContours(cstate_)
    lt_handles = cstate_('ContourHandles');
    for l =1:length(lt_handles)
        delete(lt_handles(l));
    end
    lt_handles = [];
    cstate_('ContourHandles') = lt_handles;
end

function clearBBox(cstate_)
    lt_handles = cstate_('BBox');
    for l =1:length(lt_handles)
        delete(lt_handles(l));
    end
    lt_handles = [];
    cstate_('BBox') = lt_handles;
end

function clearTriangles(cstate_)
    t_handles = cstate_('TrisHandles');
    for t = 1:length(t_handles)
        delete(t_handles(t));
    end
    cstate_('TrisHandles') = [];
end

function clearMCTriangles(cstate_)
    t_handles = cstate_('MCTrisHandles');
    for t = 1:length(t_handles)
        delete(t_handles(t));
    end
    cstate_('MCTrisHandles') = [];
end

function clearAllTriangles(cstate_)
    clearTriangles(cstate_);
    clearMCTriangles(cstate_);
end

function clearVertexConfig(cstate_)
    t_handles = cstate_('VertexConfigHandles');
    for t = 1:length(t_handles)
        delete(t_handles(t));
    end
    cstate_('VertexConfigHandles') = [];
end

function clearAsymptotes(cstate_)
    t_handles = cstate_('Asymptotes');
    for t = 1:length(t_handles)
        delete(t_handles(t));
    end
    cstate_('Asymptotes') = [];
end

%
%
function c_ = getColor(c_index)
    c_table = [[20,137,240];...
               [43,206,72];...
               [255,204,153];...
               [148,255,181];...
               [143,124,0];...
               [194,0,136];...
               [255,164,5];...
               [255,168,187];...
                [255,0,16];...
                [94,241,242]; ...
                [0,153,143]; ...
                [224,255,102]; ...
                [255,255,128];...
                [255,80,5]];

    if c_index > 14
        c_index = 1;
    end
    c_ = c_table(c_index,:)/255;
end

%
%
% Compute contours of isosurface at faces
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triangulateCell(cstate_,f,f1,coord,mc,isov)
    [nrvO,vtO,nrcO,cnsO,contours0,face,fvt,s_flag,a_cnt] = computeOrientedContours(f,mc,isov);
    %triangulateContours6InnerVertices(cstate_,0,f,mc,isov,nrvO,vtO,nrcO,cnsO,contours0,face,fvt,s_flag,a_cnt);
    triangulateSingularFace(cstate_,0,f,mc,isov,nrvO,vtO,nrcO,cnsO,contours0,face,fvt,s_flag,a_cnt);
    [nrvO,vtO,nrcO,cnsO,contours0,face,fvt,s_flag,a_cnt] = computeOrientedContours(f1,mc,isov);
    %triangulateContours6InnerVertices(cstate_,coord,f1,mc,isov,nrvO,vtO,nrcO,cnsO,contours0,face,fvt,s_flag,a_cnt);
    triangulateSingularFace(cstate_,coord,f1,mc,isov,nrvO,vtO,nrcO,cnsO,contours0,face,fvt,s_flag,a_cnt);
end

%
% Compute and plot contours
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function computeCellContours(cstate_,coord,f,f1,mc,isov)
    [nrv,vt,nrc,cns,contours] = computeContorus(f,mc,isov);

    c_handles = [];
    for c = 1:nrc
        uc = zeros(1,cns(c)+1);
        vc = zeros(1,cns(c)+1);
        wc = zeros(1,cns(c)+1);
        for i = 1:cns(c)+1
            v_index = mod(i-1,cns(c))+1;
            uc(i) = vt(contours(c,v_index),1);
            vc(i) = vt(contours(c,v_index),2);
            wc(i) = vt(contours(c,v_index),3);
        end
        l_ = plot3(uc,vc,wc,'black','LineWidth',2);
        c_handles = [c_handles,l_];
    end
    [nrv,vt,nrc,cns,contours] = computeContorus(f1,mc,isov);
    xsh = 0;
    ysh = 0;
    zsh = 0;
    if coord == 1
        xsh = 1;
    elseif coord == 2
        ysh = 1;
    elseif coord == 3
        zsh = 1;
    end
    for c = 1:nrc
        uc = zeros(1,cns(c)+1);
        vc = zeros(1,cns(c)+1);
        wc = zeros(1,cns(c)+1);
        for i = 1:cns(c)+1
            v_index = mod(i-1,cns(c))+1;
            uc(i) = vt(contours(c,v_index),1) + xsh;
            vc(i) = vt(contours(c,v_index),2) + ysh;
            wc(i) = vt(contours(c,v_index),3) + zsh;
        end
        l_ = plot3(uc,vc,wc,'magenta','LineWidth',2);
        c_handles = [c_handles,l_];
    end
    cstate_('ContourHandles') = c_handles;


end


% Compute contours
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [nrv,vt,nrc,cns,contours] = computeContorus(f,mc,isov)
    % Use the fact that the isosurface does not start or end in the cell.
    % It must cut the cell at the edges building closed hyperbolic arcs.
    % Triangulatiom
    %   i) compute intersection of level set with edges
    %  ii) compute closed polygons, we call these polygon contours
    % iii) compute the solution to the three quadratic equations
    %  iv) classify and triangulate contours

    % cell topology
    % vertices accroding to the coordinate axes, x-axes fastes
    % v1 = 000, v2 = 001, v3 = 010, v4 = 011
    % v5 = 100, v6 = 101, v7 = 110, v8 = 111
    % edges:
    % e1 = (1,2), e2 = (2,4), e3 = (3,4), e4 = (1,3)
    % e5 = (5,6), e6 = (6,8), e7 = (7,8), e8 = (5,7)
    % e9 = (1,5), e10 = (2,6), e11 = (4,8), e12 = (3,7)
    % faces:
    % f1 = (1,2,3,4), f2 = (5,6,7,8),
    % f3 = (1,2,5,6), f4 = (3,4,7,8),
    % f5 = (1,3,5,7), f6 = (2,4,6,8)

    % compute offsets along axis,
    % these are local coordinates used for interpolating vertices later on
    edges = [[1 2 -1  0  0]; ...
        [2 4  1 -1  0]; ...
        [3 4 -1  1  0]; ...
        [1 3  0 -1  0]; ...
        [5 6 -1  0  1]; ...
        [6 8  1 -1  1]; ...
        [7 8 -1  1  1]; ...
        [5 7  0 -1  1]; ...
        [1 5  0  0 -1]; ...
        [2 6  1  0 -1]; ...
        [4 8  1  1 -1]; ...
        [3 7  0  1 -1]];

    vt = zeros(12,4); % 4. coordinate means if vertex is set
    nrv = 0;
    for i = 1:12
        v1 = edges(i,1);
        v2 = edges(i,2);
        vals = [f(v1),f(v2)];
        vals = sort(vals);
        if vals(1) < isov && isov <= vals(2)
            % there is an intersection
            nrv = nrv + 1;
            iv = (isov - f(v1)) / (f(v2) - f(v1));
            if (edges(i,3) < 0)
                vt(i,1) = iv;
            else
                vt(i,1) = edges(i,3);
            end
            if (edges(i,4) < 0)
                vt(i,2) = iv;
            else
                vt(i,2) = edges(i,4);
            end
            if (edges(i,5) < 0)
                vt(i,3) = iv;
            else
                vt(i,3) = edges(i,5);
            end
            vt(i,4) = 1;
        else
            vt(i,4) = -1;
        end
    end
    % build contours, porcess facewise
    % faces
    fedge = [[1 2 3 4];[5 6 7 8];[1 10 5 9];[3 11 7 12];[4 12 8 9];[2 11 6 10]];
    segments = zeros(12,2);
    nrs = 1;
    ambiguous = zeros(1,6);
    ambcase   = zeros(1,6);
    for i = 1:6
        pos = 0;
        s = zeros(1,4);
        for j = 1:4
            if vt(fedge(i,j),4) > -1
                pos = pos + 1;
                s(pos) = fedge(i,j);
            end
        end
        % get segments from ambiguous case
        if pos > 2
            ambiguous(i) = 1; % face has ambigous case
            [u0, v0, idx] = getAsymptotes(i,f);
            % get vt1 of face
            vt1 = vt(fedge(i,1),1:3);
            if vt1(idx) < u0
                ambcase(i) = 1;
                segments(nrs,1) = s(1);
                segments(nrs,2) = s(4);
                segments(nrs+1,1) = s(2);
                segments(nrs+1,2) = s(3);
            else
                segments(nrs,1) = s(1);
                segments(nrs,2) = s(2);
                segments(nrs+1,1) = s(3);
                segments(nrs+1,2) = s(4);
            end
            nrs = nrs + 2;
        else
            segments(nrs,1) = s(1);
            segments(nrs,2) = s(2);
            nrs = nrs + 1;
        end
    end

    % compute closed contours
    contours = connect(nrv,vt,segments);
    % count nr. of countours and size
    cns = zeros(1,4);
    nrc = 0;
    for i = 1:4
        pos = 0;
        for j = 1:12
            if contours(i,j) > 0
                pos = j;
            end
        end
        cns(i) = pos;
        if pos > 0
            nrc = nrc+1;
        end
    end
end


% connect vertices into a closed contour
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function contours = connect(nrv,vt,segment)
    %
    %
    % there might be at most 4 contours with at most 12 vertices
    tcont = zeros(4,13);
    % each vertex appears at most once in only one contour
    % a contour has at least three vertices, i.e. it is a triangle or a larger
    % polygon
    ltconnect = zeros(12,2);
    vflag = ones(1,12);
    for i = 1:12
        s1 = segment(i,1);
        s2 = segment(i,2);
        if s1 > 0
            vflag(s1) = 0;
            vflag(s2) = 0;
            if ltconnect(s1,1) > 0
                ltconnect(s1,2) = s2;
            else
                ltconnect(s1,1) = s2;
            end
            if ltconnect(s2,1) > 0
                ltconnect(s2,2) = s1;
            else
                ltconnect(s2,1) = s1;
            end
        end
    end
    % create contours
    s1 = 0;
    s2 = 0;
    uvt = 0;
    % use mor data structures to be faster
    ccnt = 1;
    for v = 1:12
        if vflag(v) < 1 % vertex still not used
            s1 = v;
            tcont(ccnt,1) = v;
            pos = 2;
            orig = s1;
            next = ltconnect(orig,1);
            vflag(v) = 1; % set vertex used!
            flag = 1;
            while (flag)
                v1 = ltconnect(next,1);
                v2 = ltconnect(next,2);
                if v1 == orig
                    orig = next;
                    next = v2;
                else
                    orig = next;
                    next = v1;
                end
                tcont(ccnt,pos) = orig;
                vflag(orig) = 1; % vertex used
                pos = pos+1;
                if next == s1
                    flag = 0;
                    ccnt = ccnt + 1;
                end
            end % end of while
        end % end of if vertex not used
    end % end for over connectivity list
    contours = tcont;
end


function n_vertices = findVertexAtFace(F,cvt,face)
    f_set = false(1,6);
    vtSaddle = zeros(6,3);

    % at this point there is only one contour
    for i = 1:12
        if face(1,i) > 0
            f_set(face(1,i)) = true;
        end
    end
    % process faces
    ep = 0.000000001;
    sh = 0.03;
    %check if inner hexagon vertex is on this face 
    % modify coords of those which do not coincide with saddle point
    nr_v = 1;
    for i = 1:6
        % face 1, at this face is w=0
        [u0,v0,idx] = getAsymptotes(1,F);
        if abs(cvt(i,3)) < ep
            % check if this point is at the center of the asymptotes
            if abs(u0 - cvt(i,1)) < ep && abs(v0 - cvt(i,2)) < ep
                vtSaddle(nr_v,:) = cvt(i,1:3);
                nr_v = nr_v + 1;
            else 
                vtSaddle(nr_v,:) = cvt(i,1:3) + [0,0,sh];
                nr_v = nr_v + 1;
            end
        elseif abs(cvt(i,3)-1) < ep
            % face 2, at this face w = 1
            [u0,v0,idx] = getAsymptotes(2,F);
            if abs(u0 - cvt(i,1)) < ep && abs(v0 - cvt(i,2)) < ep
                vtSaddle(nr_v,:) = cvt(i,1:3);
                nr_v = nr_v + 1;
            else 
                vtSaddle(nr_v,:) = cvt(i,1:3) + [0,0,-sh];
                nr_v = nr_v + 1;
            end
        elseif abs(cvt(i,2)) < ep    
            % face 3
            [u0,w0,idx] = getAsymptotes(3,F);
            % check if this point is at the center of the asymptotes
            if abs(u0 - cvt(i,1)) < ep && abs(w0 - cvt(i,3)) < ep
                vtSaddle(nr_v,:) = cvt(i,1:3);
                nr_v = nr_v + 1;
            else 
                vtSaddle(nr_v,:) = cvt(i,1:3) + [0,sh,0];
                nr_v = nr_v + 1;
            end
        elseif abs(cvt(i,2) - 1) < ep
            % face 4, v = 1
            [u0,w0,idx] = getAsymptotes(4,F);
            % check if this point is at the center of the asymptotes
            if abs(u0 - cvt(i,1)) < ep && abs(w0 - cvt(i,3)) < ep
                vtSaddle(nr_v,:) = cvt(i,1:3);
                nr_v = nr_v + 1;
            else
                vtSaddle(nr_v,:) = cvt(i,1:3) + [0,-sh,0];
                nr_v = nr_v + 1;
            end
        elseif abs(cvt(i,1)) < ep
            % face 5, u = 0
            [v0,w0,idx] = getAsymptotes(5,F);
            % check if this point is at the center of the asymptotes
            if abs(v0 - cvt(i,1)) < ep && abs(w0 - cvt(i,3)) < ep
                vtSaddle(nr_v,:) = cvt(i,1:3);
                nr_v = nr_v + 1;
            else
                vtSaddle(nr_v,:) = cvt(i,1:3) + [sh,0,0];
                nr_v = nr_v + 1;
            end
        elseif abs(cvt(i,1) - 1) < ep
            % face 6, u = 1
            [v0,w0,idx] = getAsymptotes(6,F);
            % check if this point is at the center of the asymptotes
            if abs(v0 - cvt(i,1)) < ep && abs(w0 - cvt(i,3)) < ep
                vtSaddle(nr_v,:) = cvt(i,1:3);
                nr_v = nr_v + 1;
            else
                vtSaddle(nr_v,:) = cvt(i,1:3) + [-sh,0,0];
                nr_v = nr_v + 1;
            end
        else
            vtSaddle(nr_v,:) = cvt(i,1:3);
            nr_v = nr_v + 1;
        end
    end
    nr_v = nr_v - 1;
    n_vertices = zeros(nr_v,3);
    for i = 1:nr_v
        n_vertices(i,:) = vtSaddle(i,:);
    end
end

function n_vertices = findAllVertexAtFace(F,cvt,face)
    vtSaddle = zeros(6,3);
    % process faces
    ep = 0.00000001;
    %check if inner hexagon vertex is on this face 
    for i = 1:6
        if abs(cvt(i,3)) < ep
            % face 1, at this face is w=0
            [u0,v0,idx] = getAsymptotes(1,F);
            vtSaddle(i,:) = [u0,v0,0];
        elseif abs(cvt(i,3)-1) < ep
            % face 2, at this face w = 1
            [u0,v0,idx] = getAsymptotes(2,F);
            vtSaddle(i,:) = [u0,v0,1];
        elseif abs(cvt(i,2)) < ep    
            % face 3
            [u0,w0,idx] = getAsymptotes(3,F);
            vtSaddle(i,:) = [u0,0,w0];
        elseif abs(cvt(i,2) - 1) < ep
            % face 4
            [u0,w0,idx] = getAsymptotes(4,F);
            vtSaddle(i,:) = [u0,1,w0];
        elseif abs(cvt(i,1)) < ep
            % face 5
            [v0,w0,idx] = getAsymptotes(5,F);
            vtSaddle(i,:) = [0,v0,w0];
        elseif abs(cvt(i,1) - 1) < ep
            % face 6
            [v0,w0,idx] = getAsymptotes(6,F);
            vtSaddle(i,:) = [1,v0,w0];
        else
            vtSaddle(i,:) = cvt(i,1:3);
        end
    end
    n_vertices = zeros(6,3);
    for i = 1:6
        n_vertices(i,:) = vtSaddle(i,:);
    end
end

% Triangulate contours. For Tunnel and controu with 12 vertices
% use all 6 vertices of inner hexagon to produce better triangulations
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triangulateContours6InnerVertices(cstate_,coord, F,mc,i0,nrv,vt,nrc,cns,contours,face,fvt,s_flag,a_cnt)
    % find out which contours build up a single piece of surface
    nr = 0;
    ltc = 0;

    coordsh = [0,0,0];
    if coord == 1
        coordsh = [1,0,0];
    elseif coord == 2
        coordsh = [0,1,0];
    elseif coord == 3
        coordsh = [0,0,1];
    end
    % handles to geometry
    t_handles = cstate_('TrisHandles');

    % solve quadratic equations
    [ic,cvt,fs,ui,vi,wi] = quadraticEquations(F,i0);
    % check if there is a tunnel
    ict = ic(1) + ic(2) + ic(3);
    vtSaddle = cvt;
    if ict < 6
        % there is no tunnel, needs to recompute 
        % this can be done better
        [ic,cvt,fs,ui,vi,wi] = quadraticEquationsSingular(F,i0);
        ict = ic(1) + ic(2) + ic(3);
    else
        % find vertex of inner hexagon which coincides with
        % the center of the asymptotes
        vtSaddle = findAllVertexAtFace(F,cvt,face);
    end


    % triangulate tunnel
    if ict == 6
        % if there are a single contour
        % if there are three contours, one does not belong to the tunnel
        c_flag = [false,false,false];
        min_u = min(ui);
        max_u = max(ui);
        % conside only contours with 3 vertices
        for c = 1:3
            if cns(c) == 3
                umin = 2;
                umax = -2;
                for i= 1:3
                    % collect min and max u
                    if vt(contours(c,i),1) < umin
                        umin = vt(contours(c,i),1);
                    end
                    if vt(contours(c,i),1) > umax
                        umax = vt(contours(c,i),1);
                    end
                end
                % check at which side of the solutions of the quadratic
                % equations you are
                if min_u > umax || max_u < umin
                    c_flag(c) = false;
                else
                    c_flag(c) = true;
                end
            else
                c_flag(c) = true;
            end
        end
        
        % collect vertex pairs
        for cq = 1:nrc
            if c_flag(cq) == false
                trOout = zeros(1,3);
                trOout(1,:) = [1,2,3];
                vtOut = zeros(3,3);
                vtOut(1,:) = vt(contours(cq,1),1:3) + xhs;
                vtOut(2,:) = vt(contours(cq,2),1:3) + ysh;
                vtOut(3,:) = vt(contours(cq,3),1:3) + zsh;
                t_hdl = patch('Faces',trOout,'Vertices',vtOut,'FaceColor',[0.8,0.8,0.8]);
                t_handles = [t_handles,t_hdl];
                
            else
                %cns_cq = size(vtSaddle,1);
                gconne = zeros(1,12);
                for i = 1:cns(cq)
                    index = -1;
                    d = 3;
                    for j = 1:6
                        val = norm(vt(contours(cq,i),1:3)-cvt(j,:));
                        if d > val
                            index = j;
                            d = val;
                        end
                    end
                    gconne(i) = index;
                end
                % merged triangles
                % use modified nr. of vertices
                
                for i = 1:cns(cq)
                    is = i;
                    if s_flag(cq,is) == false
                        ie = mod(i,cns(cq)) + 1;
                        d1 = gconne(is);
                        d2 = gconne(ie);
                        [nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),vtSaddle);
                        for tvi = 1:size(verts,1)
                            verts(tvi,:) = verts(tvi,:) + coordsh;
                        end
                        t_hdl = patch('Faces',tris,'Vertices',verts,'FaceColor',[0.8,0.8,0.8]); %,'EdgeColor','b');
                        t_handles = [t_handles,t_hdl];
                    end
                end
                % triangulate contour with 12 vertices
                if cns(1) == 12
                    tri12 = zeros(4,3);
                    tri12(1,:) = [1,2,3];
                    tri12(2,:) = [3,4,5];
                    tri12(3,:) = [5,6,1];
                    tri12(4,:) = [1,3,5];
                    for tvi = 1:size(cvt,1)
                        cvt(tvi,:) = cvt(tvi,:) + coordsh;
                    end
                    t_hdl = patch('Faces',tri12,'Vertices',cvt,'FaceColor',[0.8,0.8,0.8]); %,'EdgeColor','b');
                    t_handles = [t_handles,t_hdl];
                end
            end
        end
    else % there is no tunnel
        % modify vertex position in case of cell 2
        for tvi = 1:size(vt,1)
            vt(tvi,1:3) = vt(tvi,1:3)+ coordsh;
        end
        for fc = 1:6
            fvt(fc,:) = fvt(fc,:) + coordsh;
        end
        if ((ict == 2 && ic(1) == 2) || (ict == 2 && ic(2) == 2)  ||  (ict == 2 && ic(3) == 2) || ict < 2)
            % there is no saddle point, this is a simple polygon
            for s = 1:nrc
                %check if contour is singular
                if a_cnt(s)
                    % star like triangulation
                    nr_c = cns(s);
                    nr_aedges = 0;
                    for aedges = 1:nr_c
                        if face(s,aedges) > 0
                            nr_aedges = nr_aedges + 1;
                        end
                    end
                    s_tris = zeros(nr_c-nr_aedges,3);
                    vta = zeros(nr_c+1,3);
                    count = 1;
                    
                    for st = 1:nr_c
                        vta(st,:) = vt(contours(s,st),1:3);
                        i0 = st;
                        i1 = mod(st,nr_c) + 1;
                        if face(s,st) < 0
                            s_tris(count,1) = i0;
                            s_tris(count,2) = i1;
                            s_tris(count,3) = nr_c+1;    
                            count = count + 1;
                        else
                            vta(nr_c+1,:) = fvt(face(s,st),:);
                        end
                    end
                    t_hdl = patch('Faces',s_tris,'Vertices',vta,'FaceColor',[0.99,0.8,0.8]);
                    t_handles = [t_handles,t_hdl];
                else
                    if cns(s) == 3
                        tri3 = zeros(1,3);
                        tri3(1,1) = 1;
                        tri3(1,2) = 2;
                        tri3(1,3) = 3;
                        vt3  = zeros(3,3);
                        vt3(1,:) = vt(contours(s,1),1:3);
                        vt3(2,:) = vt(contours(s,2),1:3);
                        vt3(3,:) = vt(contours(s,3),1:3);
                        t_hdl = patch('Faces',tri3,'Vertices',vt3,'FaceColor',[0.8,0.8,0.8]);
                        t_handles = [t_handles,t_hdl];
                    elseif cns(s) == 4
                        tri4 = zeros(2,3);
                        tri4(1,1) = 1; tri4(1,2) = 2; tri4(1,3) = 3;
                        tri4(2,1) = 1; tri4(2,2) = 3; tri4(2,3) = 4;
                        vt4  = zeros(4,3);
                        vt4(1,:) = vt(contours(s,1),1:3);
                        vt4(2,:) = vt(contours(s,2),1:3);
                        vt4(3,:) = vt(contours(s,3),1:3);
                        vt4(4,:) = vt(contours(s,4),1:3);
                        t_hdl = patch('Faces',tri4,'Vertices',vt4,'FaceColor',[0.8,0.8,0.8]);
                        t_handles = [t_handles,t_hdl];
                    elseif cns(s) == 5
                        tri5 = zeros(3,3);
                        tri5(1,1) = 1; tri5(1,2) = 2; tri5(1,3) = 3;
                        tri5(2,1) = 1; tri5(2,2) = 3; tri5(2,3) = 4;
                        tri5(3,1) = 1; tri5(3,2) = 4; tri5(3,3) = 5;
                        vt5  = zeros(5,3);
                        vt5(1,:) = vt(contours(s,1),1:3);
                        vt5(2,:) = vt(contours(s,2),1:3);
                        vt5(3,:) = vt(contours(s,3),1:3);
                        vt5(4,:) = vt(contours(s,4),1:3);
                        vt5(5,:) = vt(contours(s,5),1:3);
                        t_hdl = patch('Faces',tri5,'Vertices',vt5,'FaceColor',[0.8,0.8,0.8]);
                        t_handles = [t_handles,t_hdl];
                    elseif cns(s) == 6
                        tri6 = zeros(4,3);
                        tri6(1,1) = 1; tri6(1,2) = 2; tri6(1,3) = 3;
                        tri6(2,1) = 3; tri6(2,2) = 4; tri6(2,3) = 5;
                        tri6(3,1) = 5; tri6(3,2) = 6; tri6(3,3) = 1;
                        tri6(4,1) = 1; tri6(4,2) = 3; tri6(4,3) = 5;

                        vt6  = zeros(6,3);
                        vt6(1,:) = vt(contours(s,1),1:3);
                        vt6(2,:) = vt(contours(s,2),1:3);
                        vt6(3,:) = vt(contours(s,3),1:3);
                        vt6(4,:) = vt(contours(s,4),1:3);
                        vt6(5,:) = vt(contours(s,5),1:3);
                        vt6(6,:) = vt(contours(s,6),1:3);
                        t_hdl = patch('Faces',tri6,'Vertices',vt6,'FaceColor',[0.8,0.8,0.8]);
                        t_handles = [t_handles,t_hdl];
                    end
                end % contour is not singular
            end
        else % there are some where a saddle point
            % plot bary center
            % find out case
            fc1 = fs(1,1)*fs(2,1) + fs(1,2)*fs(2,2);
            fc2 = fs(1,1)*fs(3,1) + fs(1,2)*fs(3,2);
            fc3 = fs(2,1)*fs(3,2) + fs(2,2)*fs(3,1);
            ucoord = 0;
            vcoord = 0;
            wcoord = 0;
            c_faces = fc1 + fc2 + fc3;
            if c_faces == 2
                if fc1 == 0
                    ucoord = fs(1,1) * ui(1) + fs(1,2) * ui(2);
                    vcoord = fs(2,1) * vi(1) + fs(2,2) * vi(2);
                    wcoord = fs(2,1) * wi(2) + fs(2,2) * wi(1);
                elseif fc2 == 0
                    ucoord = fs(1,1) * ui(1) + fs(1,2) * ui(2);
                    vcoord = fs(1,1) * vi(1) + fs(1,2) * vi(2);
                    wcoord = fs(1,1) * wi(2) + fs(1,2) * wi(1);
                elseif fc3 == 0
                    ucoord = fs(2,1) * ui(1) + fs(2,2) * ui(2);
                    vcoord = fs(2,1) * vi(1) + fs(2,2) * vi(2);
                    wcoord = fs(2,1) * wi(1) + fs(2,2) * wi(2);
                end
            elseif c_faces == 3
                ucoord = (fs(1,1) * ui(1) + fs(1,2) * ui(2)) / (fs(1,1) + fs(1,2));
                vcoord = (fs(2,1) * vi(1) + fs(2,2) * vi(2)) / (fs(2,1) + fs(2,2));
                wcoord = (fs(3,1) * wi(1) + fs(3,2) * wi(2)) / (fs(3,1) + fs(3,2));
            elseif c_faces == 4
                % check which coordinate has only one solution
                nr_u = fs(1,1) + fs(1,2);
                nr_v = fs(2,1) + fs(2,2);
                nr_w = fs(3,1) + fs(3,2);
                if nr_w == 1
                    ucoord = fs(3,1) * ui(1) + fs(3,2) * ui(2);
                    vcoord = fs(3,2) * vi(1) + fs(3,1) * vi(2);
                    wcoord = fs(3,1) * wi(1) + fs(3,2) * wi(2);
                elseif nr_v == 1
                    ucoord = fs(2,1) * ui(1) + fs(2,2) * ui(2);
                    vcoord = fs(2,1) * vi(1) + fs(2,2) * vi(2);
                    wcoord = fs(2,2) * wi(1) + fs(2,1) * wi(2);
                elseif nr_u == 1
                    ucoord = fs(1,1) * ui(1) + fs(1,2) * ui(2);
                    vcoord = fs(1,1) * vi(1) + fs(1,2) * vi(2);
                    wcoord = fs(1,1) * wi(1) + fs(1,2) * wi(2);
                end
            else
                disp('ERROR: wrong nr. of solutions');
            end
            
            %
            % At this point the code is not so clear
            % we are triangulating also non ambiguous cases
            % which are handled different in the c-code
            % this is just a demonstration software
            % star triangulation
            for s = 1:nrc
                if a_cnt(s)
                    % star like triangulation
                    nr_c = cns(s);
                    nr_aedges = 0;
                    for aedges = 1:nr_c
                        if face(s,aedges) > 0
                            nr_aedges = nr_aedges + 1;
                        end
                    end
                    s_tris = zeros(nr_c-nr_aedges,3);
                    vta = zeros(nr_c+1,3);
                    count = 1;
                    
                    for st = 1:nr_c
                        vta(st,:) = vt(contours(s,st),1:3);
                        i0 = st;
                        i1 = mod(st,nr_c) + 1;
                        if face(s,st) < 0
                            s_tris(count,1) = i0;
                            s_tris(count,2) = i1;
                            s_tris(count,3) = nr_c+1;    
                            count = count + 1;
                        else
                            vta(nr_c+1,:) = fvt(face(s,st),:);
                        end
                    end
                    t_hdl = patch('Faces',s_tris,'Vertices',vta,'FaceColor',[0.99,0.8,0.8]);
                    t_handles = [t_handles,t_hdl];
                else
                    if cns(s) == 3
                        tri3 = zeros(1,3);
                        tri3(1,1) = 1;
                        tri3(1,2) = 2;
                        tri3(1,3) = 3;
                        vt3  = zeros(3,3);
                        vt3(1,:) = vt(contours(s,1),1:3);
                        vt3(2,:) = vt(contours(s,2),1:3);
                        vt3(3,:) = vt(contours(s,3),1:3);
                        t_hdl = patch('Faces',tri3,'Vertices',vt3,'FaceColor',[0.8,0.8,0.8]);
                        t_handles = [t_handles,t_hdl];

                    else
    %                     if (cns(s) == 5)
    %                         disp('triangulating case with 5 vertices');
    %                     end
                        triStar = zeros(cns(s),3);
                        vtStar  = zeros(cns(s)+1,3);
                        for xq = 1:cns(s)
                            vtStar(xq,:) = vt(contours(s,xq),1:3);
                        end
                        vtStar(cns(s)+1,:) =  [ucoord,vcoord,wcoord] + coordsh;

                        for xq=1:cns(s)
                            i0 = xq;
                            i1 = mod(xq,cns(s)) + 1;
                            triStar(xq,1) = i0;
                            triStar(xq,2) = i1;
                            triStar(xq,3) = cns(s)+1;
                        end
                        t_hdl = patch('Faces',triStar,'Vertices',vtStar,'FaceColor',[0.8,0.8,0.8]);
                        t_handles = [t_handles,t_hdl];
                    end
                end
            end

        end
        % store handles to manage geometry objects

    end

    % keep track of graphic obejcts
    cstate_('TrisHandles') = t_handles;
    % set lighting off
    for t = 1:length(t_handles)
        set(t_handles(t),'FaceLighting','none');
    end
end % triangulateContours6InnerVertices


% Triangulate contours. For Tunnel and controu with 12 vertices
% use all 6 vertices of inner hexagon to produce better triangulations
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triangulateSingularFace(cstate_,coord, F,mc,i0,nrv,vt,nrc,cns,contours,face,fvt,s_flag,a_cnt)
    % find out which contours build up a single piece of surface
    nr = 0;
    ltc = 0;

    c_pos = 2;
    if coord > 0
       c_pos = c_pos + 1;
    end
    color = getColor(c_pos);
    
    coordsh = [0,0,0];
    if coord == 1
        coordsh = [1,0,0];
    elseif coord == 2
        coordsh = [0,1,0];
    elseif coord == 3
        coordsh = [0,0,1];
    end
    % handles to geometry
    t_handles = cstate_('TrisHandles');

    % solve quadratic equations
    [ic,cvt,fs,ui,vi,wi] = quadraticEquationsAdvanced(F,i0);
    % check if there is a tunnel
    ict = ic(1) + ic(2) + ic(3);
%     vtSaddle = cvt;
%     if ict < 6
%         % there is no tunnel, needs to recompute 
%         % this can be done better
%         [ic,cvt,fs,ui,vi,wi] = quadraticEquationsSingular(F,i0);
%         ict = ic(1) + ic(2) + ic(3);
%     else
%         % find vertex of inner hexagon which coincides with
%         % the center of the asymptotes
%         vtSaddle = findVertexAtFace(F,cvt,face);
%     end

    vtSaddle = findVertexAtFace(F,cvt,face);

    % triangulate tunnel
    if ict == 6
        % if there are a single contour
        % if there are three contours, one does not belong to the tunnel
        c_flag = [false,false,false];
        min_u = min(ui);
        max_u = max(ui);
        % conside only contours with 3 vertices
        for c = 1:3
            if cns(c) == 3
                umin = 2;
                umax = -2;
                for i= 1:3
                    % collect min and max u
                    if vt(contours(c,i),1) < umin
                        umin = vt(contours(c,i),1);
                    end
                    if vt(contours(c,i),1) > umax
                        umax = vt(contours(c,i),1);
                    end
                end
                % check at which side of the solutions of the quadratic
                % equations you are
                if min_u > umax || max_u < umin
                    c_flag(c) = false;
                else
                    c_flag(c) = true;
                end
            else
                c_flag(c) = true;
            end
        end
        
        % collect vertex pairs
        for cq = 1:nrc
            if c_flag(cq) == false
                trOout = zeros(1,3);
                trOout(1,:) = [1,2,3];
                vtOut = zeros(3,3);
                vtOut(1,:) = vt(contours(cq,1),1:3) + xhs;
                vtOut(2,:) = vt(contours(cq,2),1:3) + ysh;
                vtOut(3,:) = vt(contours(cq,3),1:3) + zsh;
                t_hdl = patch('Faces',trOout,'Vertices',vtOut,'FaceColor',color);
                t_handles = [t_handles,t_hdl];
                
            else
                %cns_cq = size(vtSaddle,1);
                gconne = zeros(1,12);
                for i = 1:cns(cq)
                    index = -1;
                    d = 30;
                    for j = 1:6
                        val = norm(vt(contours(cq,i),1:3)-cvt(j,:));
                        %val = norm(vt(contours(cq,i),1:3)-vtSaddle(j,:));
                        if d > val
                            index = j;
                            d = val;
                        end
                    end
                    gconne(i) = index;
                end
                % merged triangles
                % use modified nr. of vertices
                
                for i = 1:cns(cq)
                    is = i;
                    ie = mod(i,cns(cq)) + 1;
                    d1 = gconne(is);
                    d2 = gconne(ie);
                    %[nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),cvt);
                    %[nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),cvt,false);
                    [nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),vtSaddle,s_flag(cq,is));
                    %[nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),vtSaddle,false);
                    if nrT > 0
                        %[nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),vtSaddle);
                        for tvi = 1:size(verts,1)
                            verts(tvi,:) = verts(tvi,:) + coordsh;
                        end
                        t_hdl = patch('Faces',tris,'Vertices',verts,'FaceColor',color); %,'EdgeColor','b');
                        t_handles = [t_handles,t_hdl];
                    end
                end
                % triangulate contour with 12 vertices
                if cns(1) == 12
                    tri12 = zeros(4,3);
                    tri12(1,:) = [1,2,3];
                    tri12(2,:) = [3,4,5];
                    tri12(3,:) = [5,6,1];
                    tri12(4,:) = [1,3,5];
                    for tvi = 1:size(cvt,1)
                        cvt(tvi,:) = cvt(tvi,:) + coordsh;
                    end
                    t_hdl = patch('Faces',tri12,'Vertices',cvt,'FaceColor',color); %,'EdgeColor','b');
                    t_handles = [t_handles,t_hdl];
                end
            end
        end
    else % there is no tunnel
        % modify vertex position in case of cell 2
        %for tvi = 1:size(vt,1)
        %    vt(tvi,1:3) = vt(tvi,1:3)+ coordsh;
        %end
        %for fc = 1:6
        %    fvt(fc,:) = fvt(fc,:) + coordsh;
        %end
        if ((ict == 2 && ic(1) == 2) || (ict == 2 && ic(2) == 2)  ||  (ict == 2 && ic(3) == 2) || ict < 2)
            % there is no saddle point, this is a simple polygon
            for s = 1:nrc
                %check if contour is singular
                if a_cnt(s)
                    % star like triangulation
                    nr_c = cns(s);
                    nr_aedges = 0;
                    for aedges = 1:nr_c
                        if face(s,aedges) > 0
                            nr_aedges = nr_aedges + 1;
                        end
                    end
                    s_tris = zeros(nr_c-nr_aedges,3);
                    vta = zeros(nr_c+1,3);
                    count = 1;
                    
                    for st = 1:nr_c
                        vta(st,:) = vt(contours(s,st),1:3);
                        i0 = st;
                        i1 = mod(st,nr_c) + 1;
                        if face(s,st) < 0
                            s_tris(count,1) = i0;
                            s_tris(count,2) = i1;
                            s_tris(count,3) = nr_c+1;    
                            count = count + 1;
                        else
                            vta(nr_c+1,:) = fvt(face(s,st),:);
                        end
                    end
                    t_hdl = patch('Faces',s_tris,'Vertices',vta,'FaceColor',color);
                    t_handles = [t_handles,t_hdl];
                else
                    if cns(s) == 3
                        tri3 = zeros(1,3);
                        tri3(1,1) = 1;
                        tri3(1,2) = 2;
                        tri3(1,3) = 3;
                        vt3  = zeros(3,3);
                        vt3(1,:) = vt(contours(s,1),1:3);
                        vt3(2,:) = vt(contours(s,2),1:3);
                        vt3(3,:) = vt(contours(s,3),1:3);
                        t_hdl = patch('Faces',tri3,'Vertices',vt3,'FaceColor',color);
                        t_handles = [t_handles,t_hdl];
                    elseif cns(s) == 4
                        tri4 = zeros(2,3);
                        tri4(1,1) = 1; tri4(1,2) = 2; tri4(1,3) = 3;
                        tri4(2,1) = 1; tri4(2,2) = 3; tri4(2,3) = 4;
                        vt4  = zeros(4,3);
                        vt4(1,:) = vt(contours(s,1),1:3);
                        vt4(2,:) = vt(contours(s,2),1:3);
                        vt4(3,:) = vt(contours(s,3),1:3);
                        vt4(4,:) = vt(contours(s,4),1:3);
                        t_hdl = patch('Faces',tri4,'Vertices',vt4,'FaceColor',color);
                        t_handles = [t_handles,t_hdl];
                    elseif cns(s) == 5
                        tri5 = zeros(3,3);
                        tri5(1,1) = 1; tri5(1,2) = 2; tri5(1,3) = 3;
                        tri5(2,1) = 1; tri5(2,2) = 3; tri5(2,3) = 4;
                        tri5(3,1) = 1; tri5(3,2) = 4; tri5(3,3) = 5;
                        vt5  = zeros(5,3);
                        vt5(1,:) = vt(contours(s,1),1:3);
                        vt5(2,:) = vt(contours(s,2),1:3);
                        vt5(3,:) = vt(contours(s,3),1:3);
                        vt5(4,:) = vt(contours(s,4),1:3);
                        vt5(5,:) = vt(contours(s,5),1:3);
                        t_hdl = patch('Faces',tri5,'Vertices',vt5,'FaceColor',color);
                        t_handles = [t_handles,t_hdl];
                    elseif cns(s) == 6
                        tri6 = zeros(4,3);
                        tri6(1,1) = 1; tri6(1,2) = 2; tri6(1,3) = 3;
                        tri6(2,1) = 3; tri6(2,2) = 4; tri6(2,3) = 5;
                        tri6(3,1) = 5; tri6(3,2) = 6; tri6(3,3) = 1;
                        tri6(4,1) = 1; tri6(4,2) = 3; tri6(4,3) = 5;

                        vt6  = zeros(6,3);
                        vt6(1,:) = vt(contours(s,1),1:3);
                        vt6(2,:) = vt(contours(s,2),1:3);
                        vt6(3,:) = vt(contours(s,3),1:3);
                        vt6(4,:) = vt(contours(s,4),1:3);
                        vt6(5,:) = vt(contours(s,5),1:3);
                        vt6(6,:) = vt(contours(s,6),1:3);
                        t_hdl = patch('Faces',tri6,'Vertices',vt6,'FaceColor',color);
                        t_handles = [t_handles,t_hdl];
                    end
                end % contour is not singular
            end
        else % there are some where a saddle point
            % find out how many singular saddle points there are
            nr_singular = 0;
            for i = 1:6
                if fvt(i,1) > -1
                    nr_singular = nr_singular + 1;
                end
            end
            singular_s = zeros(nr_singular,3);
            sig_pos = 1;
            for i = 1:6
                if fvt(i,1) > -1
                    singular_s(sig_pos,:) = fvt(i,:);
                    sig_pos = sig_pos + 1;
                end
            end
            
            %loop over all contours, this is wrong
            for cq = 1:nrc
                gconne = zeros(1,12);
                for i = 1:cns(cq)
                    index = -1;
                    d = 3;
                    for j = 1:nr_singular
                        val = norm(vt(contours(cq,i),1:3)-singular_s(j,:));
                        if d > val
                            index = j;
                            d = val;
                        end
                    end
                    gconne(i) = index;
                end
                % merged triangles
                % use modified nr. of vertices

                for i = 1:cns(cq)
                    is = i;
                    ie = mod(i,cns(cq)) + 1;
                    d1 = gconne(is);
                    d2 = gconne(ie);
                    %[nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),cvt);
                    %[nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),cvt,s_flag(cq,is));
                    [nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),singular_s,s_flag(cq,is));
                    if nrT > 0
                        %[nrT,tris,verts] = triangulateVts(is,ie,d1,d2,vt,contours(cq,:),vtSaddle);
                        for tvi = 1:size(verts,1)
                            verts(tvi,:) = verts(tvi,:) + coordsh;
                        end
                        t_hdl = patch('Faces',tris,'Vertices',verts,'FaceColor',color); %,'EdgeColor','b');
                        t_handles = [t_handles,t_hdl];
                    end
                end
            end
        end
        % store handles to manage geometry objects

    end

    % keep track of graphic obejcts
    cstate_('TrisHandles') = t_handles;
    % set lighting off
    for t = 1:length(t_handles)
        set(t_handles(t),'FaceLighting','none');
    end
end % triangulateSingularFace



% compute the asymptotic decider
function alp = aDecider(f1,f2,f3,f4)
    eta = f1+f4-f2-f3;
    alp = (f1*f4 - f2*f3) / eta;
end

%
% Compute asymptotes at faces
% The functions accept arrays as input
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uw,vw,ew,aw] = wSurface(f,w)
    f0 = f(1);
    f1 = f(2);
    f2 = f(3);
    f3 = f(4);
    f4 = f(5);
    f5 = f(6);
    f6 = f(7);
    f7 = f(8);
    if w == 0
        ew = (f0+f3-f1-f2);
        uw = (f0-f2) ./ ew;
        vw = (f0-f1) ./ ew;
        aw = ((f0 .* f3) - (f1 .* f2)) ./ ew;
    elseif w == 1
        ew = f4+f7-f5-f6;
        uw = (f4-f6) ./ ew;
        vw = (f4-f5) ./ ew;
        aw = ((f4 .* f7) - (f5 .* f6)) ./ ew;
    else 
        ew = (f0+f3-f1-f2)*(1-w) + (f4+f7-f5-f6)*w;
        uw = ((f0-f2)*(1-w)+(f4-f6)*w) ./ ew;
        vw = ((f0-f1)*(1-w)+(f4-f5)*w) ./ ew;
        aw = ((f0*(1-w)+f4*w).*(f3*(1-w)+f7*w) - (f1*(1-w)+f5*w).*(f2*(1-w)+f6*w)) ./ ew;
    end
end
function [vu,wu,eu,au] = uSurface(f,u)
    f0 = f(1);
    f1 = f(2);
    f2 = f(3);
    f3 = f(4);
    f4 = f(5);
    f5 = f(6);
    f6 = f(7);
    f7 = f(8);
    if u == 0
        eu = (f0+f6-f2-f4);
        vu = (f0-f4) ./ eu;
        wu = (f0-f2) ./ eu;
        au = ((f0.*f6) - (f2.*f4)) ./ eu;
    elseif u == 1
        eu = (f1+f7-f3-f5);
        vu = (f1-f5) ./ eu;
        wu = (f1-f3) ./ eu;
        au = ((f1 .* f7) - (f3 .* f5)) ./ eu;
    else
        eu = (f0+f6-f2-f4)*(1-u) + (f1+f7-f3-f5)*u;
        vu = ((f0-f4)*(1-u)+(f1-f5)*u) ./ eu;
        wu = ((f0-f2)*(1-u)+(f1-f3)*u) ./ eu;
        au = ((f0*(1-u)+f1*u).*(f6*(1-u)+f7*u) - (f2*(1-u)+f3*u).*(f4*(1-u)+f5*u)) ./ eu;
    end
end
function [uv,wv,ev,av] = vSurface(f,v)
    f0 = f(1);
    f1 = f(2);
    f2 = f(3);
    f3 = f(4);
    f4 = f(5);
    f5 = f(6);
    f6 = f(7);
    f7 = f(8);
    if v == 0
        ev = (f0+f5-f1-f4);
        uv = (f0-f4) ./ ev;
        wv = (f0-f1) ./ ev;
        av = ((f0 .* f5) - (f1 .* f4)) ./ ev;
    elseif v == 1
        ev = (f2+f7-f3-f6);
        uv = (f2-f6) ./ ev;
        wv = (f2-f3) ./ ev;
        av = ((f2 .* f7) - (f3 .*f6)) ./ ev;
    else
        ev = (f0+f5-f1-f4)*(1-v) + (f2+f7-f3-f6)*v;
        uv = ((f0-f4)*(1-v)+(f2-f6)*v) ./ ev;
        wv = ((f0-f1)*(1-v)+(f2-f3)*v) ./ ev;
        av = ((f0*(1-v)+f2*v).*(f5*(1-v)+f7*v) - (f1*(1-v)+f3*v).*(f4*(1-v)+f6*v)) ./ ev;
    end
end

% computes u0 and v0 where the asymptotes intersect
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u0,v0,idx] = getAsymptotes(face,f)
    switch face
        case 1
            [u0, v0, eta, alpha] = wSurface(f,0);
            idx = 1;
        case 2
            [u0, v0, eta, alpha] = wSurface(f,1);
            idx = 1;
        case 3
            [u0, v0, eta, alpha] = vSurface(f,0);
            idx = 1;
        case 4
            [u0, v0, eta, alpha] = vSurface(f,1);
            idx = 1;
        case 5
            [u0, v0, eta, alpha] = uSurface(f,0);
            idx = 2;
        case 6
            [u0, v0, eta, alpha] = uSurface(f,1);
            idx = 2;
        otherwise disp('wrong face id')
    end
end

%
% select a rotation symmetry
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fres = transform(index,f,posSym_)
    % set of 24 rotations which are symmetries in MC
    rot = zeros(24,8);
    rot(1,:) = [1,2,3,4,5,6,7,8];
    rot(2,:) = [3,4,7,8,1,2,5,6];
    rot(3,:) = [7,8,5,6,3,4,1,2];
    rot(4,:) = [5,6,1,2,7,8,3,4];
    rot(5,:) = [5,1,7,3,6,2,8,4];
    rot(6,:) = [6,5,8,7,2,1,4,3];
    rot(7,:) = [2,6,4,8,1,5,3,7];
    rot(8,:) = [2,4,1,3,6,8,5,7];
    rot(9,:) = [4,3,2,1,8,7,6,5];
    rot(10,:) = [3,1,4,2,7,5,8,6];
    rot(11,:) = [2,1,6,5,4,3,8,7];
    rot(12,:) = [8,7,4,3,6,5,2,1];
    rot(13,:) = [3,7,1,5,4,8,2,6];
    rot(14,:) = [8,4,6,2,7,3,5,1];
    rot(15,:) = [8,6,7,5,4,2,3,1];
    rot(16,:) = [5,7,6,8,1,3,2,4];
    rot(17,:) = [1,3,5,7,2,4,6,8];
    rot(18,:) = [1,5,2,6,3,7,4,8];
    rot(19,:) = [6,8,2,4,5,7,1,3];
    rot(20,:) = [7,3,8,4,5,1,6,2];
    rot(21,:) = [6,2,5,1,8,4,7,3];
    rot(22,:) = [4,2,8,6,3,1,7,5];
    rot(23,:) = [4,8,3,7,2,6,1,5];
    rot(24,:) = [7,5,3,1,8,6,4,2];

    % transform f, unroll
    ft = zeros(1,8);
    for i=1:8
        ft(i) = f(rot(index,i));
    end
    if posSym_ == 1
       fres = -1 * ft;
    else 
       fres = ft;
    end
end


% the standard MC
function s_mc(cstate_,coord,f,f1,mc,i0)
   % define large tables
   % edge intersection pattern for MC case
e_pattern = uint16([ ...
    0, 265, 515, 778, 2060, 2309, 2575, 2822, ...
    1030, 1295, 1541, 1804, 3082, 3331, 3593, 3840, ...
    400, 153, 915, 666, 2460, 2197, 2975, 2710, ...
    1430, 1183, 1941, 1692, 3482, 3219, 3993, 3728, ...
    560, 825, 51, 314, 2620, 2869, 2111, 2358, ...
    1590, 1855, 1077, 1340, 3642, 3891, 3129, 3376, ...
    928, 681, 419, 170, 2988, 2725, 2479, 2214, ...
    1958, 1711, 1445, 1196, 4010, 3747, 3497, 3232, ...
    2240, 2505, 2755, 3018, 204, 453, 719, 966, ...
    3270, 3535, 3781, 4044, 1226, 1475, 1737, 1984, ...
    2384, 2137, 2899, 2650, 348, 85, 863, 598, ...
    3414, 3167, 3925, 3676, 1370, 1107, 1881, 1616, ...
    2800, 3065, 2291, 2554, 764, 1013, 255, 502, ...
    3830, 4095, 3317, 3580, 1786, 2035, 1273, 1520, ...
    2912, 2665, 2403, 2154, 876, 613, 367, 102, ...
    3942, 3695, 3429, 3180, 1898, 1635, 1385, 1120, ...
    1120, 1385, 1635, 1898, 3180, 3429, 3695, 3942, ...
    102, 367, 613, 876, 2154, 2403, 2665, 2912, ...
    1520, 1273, 2035, 1786, 3580, 3317, 4095, 3830, ...
    502, 255, 1013, 764, 2554, 2291, 3065, 2800, ...
    1616, 1881, 1107, 1370, 3676, 3925, 3167, 3414, ...
    598, 863, 85, 348, 2650, 2899, 2137, 2384, ...
    1984, 1737, 1475, 1226, 4044, 3781, 3535, 3270, ...
    966, 719, 453, 204, 3018, 2755, 2505, 2240, ...
    3232, 3497, 3747, 4010, 1196, 1445, 1711, 1958, ...
    2214, 2479, 2725, 2988, 170, 419, 681, 928, ...
    3376, 3129, 3891, 3642, 1340, 1077, 1855, 1590, ...
    2358, 2111, 2869, 2620, 314, 51, 825, 560, ...
    3728, 3993, 3219, 3482, 1692, 1941, 1183, 1430, ...
    2710, 2975, 2197, 2460, 666, 915, 153, 400, ...
    3840, 3593, 3331, 3082, 1804, 1541, 1295, 1030,...
    2822, 2575, 2309, 2060, 778, 515, 265, 0]);


% list of triangles for Marching Cubes case t, position at t*16 + tri
t_pattern = int32([ ...
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 0 <-> mc: 0, class rep: 0
    0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 1 <-> mc: 1, class rep: 1
    0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 2 <-> mc: 2, class rep: 1
    1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 3 <-> mc: 3, class rep: 3
    3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 4 <-> mc: 8, class rep: 1
    0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 5 <-> mc: 9, class rep: 3
    1, 0, 9, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 6 <-> mc: 10, class rep: 6
    1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 7 <-> mc: 11, class rep: 7
    1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 8 <-> mc: 4, class rep: 1
    0, 3, 8, 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 9 <-> mc: 5, class rep: 6
    9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 10 <-> mc: 6, class rep: 3
    2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 11 <-> mc: 7, class rep: 7
    3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 12 <-> mc: 12, class rep: 3
    0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 13 <-> mc: 13, class rep: 7
    3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 14 <-> mc: 14, class rep: 7
    9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 15 <-> mc: 15, class rep: 15
    4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 16 <-> mc: 16, class rep: 1
    4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 17 <-> mc: 17, class rep: 3
    0, 9, 1, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 18 <-> mc: 18, class rep: 6
    4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 19 <-> mc: 19, class rep: 7
    8, 7, 4, 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 20 <-> mc: 24, class rep: 6
    11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 21 <-> mc: 25, class rep: 7
    9, 1, 0, 8, 7, 4, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 22 <-> mc: 26, class rep: 22
    4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1, ... %%quitte: 23 <-> mc: 27, class rep: 23
    1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 24 <-> mc: 20, class rep: 24
    3, 7, 4, 3, 4, 0, 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 25 <-> mc: 21, class rep: 25
    9, 10, 2, 9, 2, 0, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 26 <-> mc: 22, class rep: 25
    2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1, ... %%quitte: 27 <-> mc: 23, class rep: 27
    3, 1, 10, 3, 10, 11, 7, 4, 8, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 28 <-> mc: 28, class rep: 25
    1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1, ... %%quitte: 29 <-> mc: 29, class rep: 29
    4, 8, 7, 9, 11, 0, 9, 10, 11, 11, 3, 0, -1, -1, -1, -1, ... %%quitte: 30 <-> mc: 30, class rep: 30
    4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 31 <-> mc: 31, class rep: 7
    9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 32 <-> mc: 32, class rep: 1
    9, 4, 5, 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 33 <-> mc: 33, class rep: 6
    0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 34 <-> mc: 34, class rep: 3
    8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 35 <-> mc: 35, class rep: 7
    9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 36 <-> mc: 40, class rep: 24
    0, 2, 11, 0, 11, 8, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 37 <-> mc: 41, class rep: 25
    0, 4, 5, 0, 5, 1, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 38 <-> mc: 42, class rep: 25
    2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1, ... %%quitte: 39 <-> mc: 43, class rep: 29
    1, 10, 2, 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 40 <-> mc: 36, class rep: 6
    3, 8, 0, 1, 10, 2, 4, 5, 9, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 41 <-> mc: 37, class rep: 22
    5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 42 <-> mc: 38, class rep: 7
    2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1, ... %%quitte: 43 <-> mc: 39, class rep: 23
    10, 11, 3, 10, 3, 1, 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 44 <-> mc: 44, class rep: 25
    4, 5, 9, 0, 1, 8, 8, 1, 10, 8, 10, 11, -1, -1, -1, -1, ... %%quitte: 45 <-> mc: 45, class rep: 30
    5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1, ... %%quitte: 46 <-> mc: 46, class rep: 27
    5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 47 <-> mc: 47, class rep: 7
    9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 48 <-> mc: 48, class rep: 3
    9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 49 <-> mc: 49, class rep: 7
    0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 50 <-> mc: 50, class rep: 7
    1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 51 <-> mc: 51, class rep: 15
    7, 5, 9, 7, 9, 8, 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 52 <-> mc: 56, class rep: 25
    9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1, ... %%quitte: 53 <-> mc: 57, class rep: 27
    2, 11, 3, 0, 8, 1, 1, 8, 7, 1, 7, 5, -1, -1, -1, -1, ... %%quitte: 54 <-> mc: 58, class rep: 30
    11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 55 <-> mc: 59, class rep: 7
    9, 8, 7, 9, 7, 5, 10, 2, 1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 56 <-> mc: 52, class rep: 25
    10, 2, 1, 9, 0, 5, 5, 0, 3, 5, 3, 7, -1, -1, -1, -1, ... %%quitte: 57 <-> mc: 53, class rep: 30
    8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1, ... %%quitte: 58 <-> mc: 54, class rep: 29
    2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 59 <-> mc: 55, class rep: 7
    9, 8, 5, 8, 7, 5, 10, 3, 1, 10, 11, 3, -1, -1, -1, -1, ... %%quitte: 60 <-> mc: 60, class rep: 60
    5, 0, 7, 5, 9, 0, 7, 0, 11, 1, 10, 0, 11, 0, 10, -1, ... %%quitte: 61 <-> mc: 61, class rep: 25
    11, 0, 10, 11, 3, 0, 10, 0, 5, 8, 7, 0, 5, 0, 7, -1, ... %%quitte: 62 <-> mc: 62, class rep: 25
    11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 63 <-> mc: 63, class rep: 3
    7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 64 <-> mc: 128, class rep: 1
    3, 8, 0, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 65 <-> mc: 129, class rep: 6
    0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 66 <-> mc: 130, class rep: 24
    8, 9, 1, 8, 1, 3, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 67 <-> mc: 131, class rep: 25
    7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 68 <-> mc: 136, class rep: 3
    7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 69 <-> mc: 137, class rep: 7
    2, 6, 7, 2, 7, 3, 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 70 <-> mc: 138, class rep: 25
    1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1, ... %%quitte: 71 <-> mc: 139, class rep: 27
    10, 2, 1, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 72 <-> mc: 132, class rep: 6
    1, 10, 2, 3, 8, 0, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 73 <-> mc: 133, class rep: 22
    2, 0, 9, 2, 9, 10, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 74 <-> mc: 134, class rep: 25
    6, 7, 11, 2, 3, 10, 10, 3, 8, 10, 8, 9, -1, -1, -1, -1, ... %%quitte: 75 <-> mc: 135, class rep: 30
    10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 76 <-> mc: 140, class rep: 7
    10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1, ... %%quitte: 77 <-> mc: 141, class rep: 23
    0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1, ... %%quitte: 78 <-> mc: 142, class rep: 29
    7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 79 <-> mc: 143, class rep: 7
    6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 80 <-> mc: 144, class rep: 3
    3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 81 <-> mc: 145, class rep: 7
    8, 11, 6, 8, 6, 4, 9, 1, 0, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 82 <-> mc: 146, class rep: 25
    9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1, ... %%quitte: 83 <-> mc: 147, class rep: 29
    8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 84 <-> mc: 152, class rep: 7
    0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 85 <-> mc: 153, class rep: 15
    1, 0, 9, 2, 4, 3, 2, 6, 4, 4, 8, 3, -1, -1, -1, -1, ... %%quitte: 86 <-> mc: 154, class rep: 30
    1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 87 <-> mc: 155, class rep: 7
    6, 4, 8, 6, 8, 11, 2, 1, 10, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 88 <-> mc: 148, class rep: 25
    1, 10, 2, 3, 11, 0, 0, 11, 6, 0, 6, 4, -1, -1, -1, -1, ... %%quitte: 89 <-> mc: 149, class rep: 30
    4, 8, 11, 4, 11, 6, 0, 9, 2, 2, 9, 10, -1, -1, -1, -1, ... %%quitte: 90 <-> mc: 150, class rep: 60
    10, 3, 9, 10, 2, 3, 9, 3, 4, 11, 6, 3, 4, 3, 6, -1, ... %%quitte: 91 <-> mc: 151, class rep: 25
    8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1, ... %%quitte: 92 <-> mc: 156, class rep: 27
    10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 93 <-> mc: 157, class rep: 7
    4, 3, 6, 4, 8, 3, 6, 3, 10, 0, 9, 3, 10, 3, 9, -1, ... %%quitte: 94 <-> mc: 158, class rep: 25
    10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 95 <-> mc: 159, class rep: 3
    4, 5, 9, 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 96 <-> mc: 160, class rep: 6
    0, 3, 8, 4, 5, 9, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 97 <-> mc: 161, class rep: 22
    5, 1, 0, 5, 0, 4, 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 98 <-> mc: 162, class rep: 25
    11, 6, 7, 8, 4, 3, 3, 4, 5, 3, 5, 1, -1, -1, -1, -1, ... %%quitte: 99 <-> mc: 163, class rep: 30
    7, 3, 2, 7, 2, 6, 5, 9, 4, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 100 <-> mc: 168, class rep: 25
    9, 4, 5, 0, 6, 8, 0, 2, 6, 6, 7, 8, -1, -1, -1, -1, ... %%quitte: 101 <-> mc: 169, class rep: 30
    3, 2, 6, 3, 6, 7, 1, 0, 5, 5, 0, 4, -1, -1, -1, -1, ... %%quitte: 102 <-> mc: 170, class rep: 60
    6, 8, 2, 6, 7, 8, 2, 8, 1, 4, 5, 8, 1, 8, 5, -1, ... %%quitte: 103 <-> mc: 171, class rep: 25
    9, 4, 5, 10, 2, 1, 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 104 <-> mc: 164, class rep: 22
    6, 7, 11, 1, 10, 2, 0, 3, 8, 4, 5, 9, -1, -1, -1, -1, ... %%quitte: 105 <-> mc: 165, class rep: 105
    7, 11, 6, 5, 10, 4, 4, 10, 2, 4, 2, 0, -1, -1, -1, -1, ... %%quitte: 106 <-> mc: 166, class rep: 30
    3, 8, 4, 3, 4, 5, 3, 5, 2, 10, 2, 5, 11, 6, 7, -1, ... %%quitte: 107 <-> mc: 167, class rep: 22
    9, 4, 5, 10, 6, 1, 1, 6, 7, 1, 7, 3, -1, -1, -1, -1, ... %%quitte: 108 <-> mc: 172, class rep: 30
    1, 10, 6, 1, 6, 7, 1, 7, 0, 8, 0, 7, 9, 4, 5, -1, ... %%quitte: 109 <-> mc: 173, class rep: 22
    4, 10, 0, 4, 5, 10, 0, 10, 3, 6, 7, 10, 3, 10, 7, -1, ... %%quitte: 110 <-> mc: 174, class rep: 25
    7, 10, 6, 7, 8, 10, 5, 10, 4, 4, 10, 8, -1, -1, -1, -1, ... %%quitte: 111 <-> mc: 175, class rep: 6
    6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 112 <-> mc: 176, class rep: 7
    3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1, ... %%quitte: 113 <-> mc: 177, class rep: 23
    0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1, ... %%quitte: 114 <-> mc: 178, class rep: 27
    6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 115 <-> mc: 179, class rep: 7
    5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1, ... %%quitte: 116 <-> mc: 184, class rep: 29
    9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 117 <-> mc: 185, class rep: 7
    1, 8, 5, 1, 0, 8, 5, 8, 6, 3, 2, 8, 6, 8, 2, -1, ... %%quitte: 118 <-> mc: 186, class rep: 25
    1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 119 <-> mc: 187, class rep: 3
    1, 10, 2, 9, 11, 5, 9, 8, 11, 11, 6, 5, -1, -1, -1, -1, ... %%quitte: 120 <-> mc: 180, class rep: 30
    0, 3, 11, 0, 11, 6, 0, 6, 9, 5, 9, 6, 1, 10, 2, -1, ... %%quitte: 121 <-> mc: 181, class rep: 22
    11, 5, 8, 11, 6, 5, 8, 5, 0, 10, 2, 5, 0, 5, 2, -1, ... %%quitte: 122 <-> mc: 182, class rep: 25
    6, 3, 11, 6, 5, 3, 2, 3, 10, 10, 3, 5, -1, -1, -1, -1, ... %%quitte: 123 <-> mc: 183, class rep: 6
    1, 6, 3, 1, 10, 6, 3, 6, 8, 5, 9, 6, 8, 6, 9, -1, ... %%quitte: 124 <-> mc: 188, class rep: 25
    10, 0, 1, 10, 6, 0, 9, 0, 5, 5, 0, 6, -1, -1, -1, -1, ... %%quitte: 125 <-> mc: 189, class rep: 6
    0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 126 <-> mc: 190, class rep: 24
    10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 127 <-> mc: 191, class rep: 1
    10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 128 <-> mc: 64, class rep: 1
    0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 129 <-> mc: 65, class rep: 24
    9, 1, 0, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 130 <-> mc: 66, class rep: 6
    1, 3, 8, 1, 8, 9, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 131 <-> mc: 67, class rep: 25
    2, 11, 3, 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 132 <-> mc: 72, class rep: 6
    11, 8, 0, 11, 0, 2, 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 133 <-> mc: 73, class rep: 25
    0, 9, 1, 2, 11, 3, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 134 <-> mc: 74, class rep: 22
    5, 6, 10, 1, 2, 9, 9, 2, 11, 9, 11, 8, -1, -1, -1, -1, ... %%quitte: 135 <-> mc: 75, class rep: 30
    1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 136 <-> mc: 68, class rep: 3
    1, 5, 6, 1, 6, 2, 3, 8, 0, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 137 <-> mc: 69, class rep: 25
    9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 138 <-> mc: 70, class rep: 7
    5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1, ... %%quitte: 139 <-> mc: 71, class rep: 29
    6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 140 <-> mc: 76, class rep: 7
    0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1, ... %%quitte: 141 <-> mc: 77, class rep: 27
    3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1, ... %%quitte: 142 <-> mc: 78, class rep: 23
    6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 143 <-> mc: 79, class rep: 7
    5, 6, 10, 4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 144 <-> mc: 80, class rep: 6
    4, 0, 3, 4, 3, 7, 6, 10, 5, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 145 <-> mc: 81, class rep: 25
    1, 0, 9, 5, 6, 10, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 146 <-> mc: 82, class rep: 22
    10, 5, 6, 1, 7, 9, 1, 3, 7, 7, 4, 9, -1, -1, -1, -1, ... %%quitte: 147 <-> mc: 83, class rep: 30
    3, 2, 11, 7, 4, 8, 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 148 <-> mc: 88, class rep: 22
    5, 6, 10, 4, 2, 7, 4, 0, 2, 2, 11, 7, -1, -1, -1, -1, ... %%quitte: 149 <-> mc: 89, class rep: 30
    0, 9, 1, 4, 8, 7, 2, 11, 3, 5, 6, 10, -1, -1, -1, -1, ... %%quitte: 150 <-> mc: 90, class rep: 105
    9, 1, 2, 9, 2, 11, 9, 11, 4, 7, 4, 11, 5, 6, 10, -1, ... %%quitte: 151 <-> mc: 91, class rep: 22
    6, 2, 1, 6, 1, 5, 4, 8, 7, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 152 <-> mc: 84, class rep: 25
    1, 5, 2, 5, 6, 2, 3, 4, 0, 3, 7, 4, -1, -1, -1, -1, ... %%quitte: 153 <-> mc: 85, class rep: 60
    8, 7, 4, 9, 5, 0, 0, 5, 6, 0, 6, 2, -1, -1, -1, -1, ... %%quitte: 154 <-> mc: 86, class rep: 30
    7, 9, 3, 7, 4, 9, 3, 9, 2, 5, 6, 9, 2, 9, 6, -1, ... %%quitte: 155 <-> mc: 87, class rep: 25
    8, 7, 4, 3, 5, 11, 3, 1, 5, 5, 6, 11, -1, -1, -1, -1, ... %%quitte: 156 <-> mc: 92, class rep: 30
    5, 11, 1, 5, 6, 11, 1, 11, 0, 7, 4, 11, 0, 11, 4, -1, ... %%quitte: 157 <-> mc: 93, class rep: 25
    0, 9, 5, 0, 5, 6, 0, 6, 3, 11, 3, 6, 8, 7, 4, -1, ... %%quitte: 158 <-> mc: 94, class rep: 22
    6, 9, 5, 6, 11, 9, 4, 9, 7, 7, 9, 11, -1, -1, -1, -1, ... %%quitte: 159 <-> mc: 95, class rep: 6
    10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 160 <-> mc: 96, class rep: 3
    4, 6, 10, 4, 10, 9, 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 161 <-> mc: 97, class rep: 25
    10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 162 <-> mc: 98, class rep: 7
    8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1, ... %%quitte: 163 <-> mc: 99, class rep: 27
    10, 9, 4, 10, 4, 6, 11, 3, 2, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 164 <-> mc: 104, class rep: 25
    0, 2, 8, 2, 11, 8, 4, 10, 9, 4, 6, 10, -1, -1, -1, -1, ... %%quitte: 165 <-> mc: 105, class rep: 60
    3, 2, 11, 0, 6, 1, 0, 4, 6, 6, 10, 1, -1, -1, -1, -1, ... %%quitte: 166 <-> mc: 106, class rep: 30
    6, 1, 4, 6, 10, 1, 4, 1, 8, 2, 11, 1, 8, 1, 11, -1, ... %%quitte: 167 <-> mc: 107, class rep: 25
    1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 168 <-> mc: 100, class rep: 7
    3, 8, 0, 1, 9, 2, 2, 9, 4, 2, 4, 6, -1, -1, -1, -1, ... %%quitte: 169 <-> mc: 101, class rep: 30
    0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 170 <-> mc: 102, class rep: 15
    8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 171 <-> mc: 103, class rep: 7
    9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1, ... %%quitte: 172 <-> mc: 108, class rep: 29
    8, 1, 11, 8, 0, 1, 11, 1, 6, 9, 4, 1, 6, 1, 4, -1, ... %%quitte: 173 <-> mc: 109, class rep: 25
    3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 174 <-> mc: 110, class rep: 7
    6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 175 <-> mc: 111, class rep: 3
    7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 176 <-> mc: 112, class rep: 7
    0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1, ... %%quitte: 177 <-> mc: 113, class rep: 29
    10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1, ... %%quitte: 178 <-> mc: 114, class rep: 23
    10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 179 <-> mc: 115, class rep: 7
    2, 11, 3, 10, 8, 6, 10, 9, 8, 8, 7, 6, -1, -1, -1, -1, ... %%quitte: 180 <-> mc: 120, class rep: 30
    2, 7, 0, 2, 11, 7, 0, 7, 9, 6, 10, 7, 9, 7, 10, -1, ... %%quitte: 181 <-> mc: 121, class rep: 25
    1, 0, 8, 1, 8, 7, 1, 7, 10, 6, 10, 7, 2, 11, 3, -1, ... %%quitte: 182 <-> mc: 122, class rep: 22
    11, 1, 2, 11, 7, 1, 10, 1, 6, 6, 1, 7, -1, -1, -1, -1, ... %%quitte: 183 <-> mc: 123, class rep: 6
    1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1, ... %%quitte: 184 <-> mc: 116, class rep: 27
    2, 9, 6, 2, 1, 9, 6, 9, 7, 0, 3, 9, 7, 9, 3, -1, ... %%quitte: 185 <-> mc: 117, class rep: 25
    7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 186 <-> mc: 118, class rep: 7
    7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 187 <-> mc: 119, class rep: 3
    8, 6, 9, 8, 7, 6, 9, 6, 1, 11, 3, 6, 1, 6, 3, -1, ... %%quitte: 188 <-> mc: 124, class rep: 25
    0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 189 <-> mc: 125, class rep: 24
    7, 0, 8, 7, 6, 0, 3, 0, 11, 11, 0, 6, -1, -1, -1, -1, ... %%quitte: 190 <-> mc: 126, class rep: 6
    7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 191 <-> mc: 127, class rep: 1
    11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 192 <-> mc: 192, class rep: 3
    11, 10, 5, 11, 5, 7, 8, 0, 3, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 193 <-> mc: 193, class rep: 25
    5, 7, 11, 5, 11, 10, 1, 0, 9, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 194 <-> mc: 194, class rep: 25
    10, 5, 7, 10, 7, 11, 9, 1, 8, 8, 1, 3, -1, -1, -1, -1, ... %%quitte: 195 <-> mc: 195, class rep: 60
    2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 196 <-> mc: 200, class rep: 7
    8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1, ... %%quitte: 197 <-> mc: 201, class rep: 29
    9, 1, 0, 5, 3, 10, 5, 7, 3, 3, 2, 10, -1, -1, -1, -1, ... %%quitte: 198 <-> mc: 202, class rep: 30
    9, 2, 8, 9, 1, 2, 8, 2, 7, 10, 5, 2, 7, 2, 5, -1, ... %%quitte: 199 <-> mc: 203, class rep: 25
    11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 200 <-> mc: 196, class rep: 7
    0, 3, 8, 1, 7, 2, 1, 5, 7, 7, 11, 2, -1, -1, -1, -1, ... %%quitte: 201 <-> mc: 197, class rep: 30
    9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1, ... %%quitte: 202 <-> mc: 198, class rep: 27
    7, 2, 5, 7, 11, 2, 5, 2, 9, 3, 8, 2, 9, 2, 8, -1, ... %%quitte: 203 <-> mc: 199, class rep: 25
    1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 204 <-> mc: 204, class rep: 15
    0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 205 <-> mc: 205, class rep: 7
    9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 206 <-> mc: 206, class rep: 7
    9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 207 <-> mc: 207, class rep: 3
    5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 208 <-> mc: 208, class rep: 7
    5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1, ... %%quitte: 209 <-> mc: 209, class rep: 27
    0, 9, 1, 8, 10, 4, 8, 11, 10, 10, 5, 4, -1, -1, -1, -1, ... %%quitte: 210 <-> mc: 210, class rep: 30
    10, 4, 11, 10, 5, 4, 11, 4, 3, 9, 1, 4, 3, 4, 1, -1, ... %%quitte: 211 <-> mc: 211, class rep: 25
    2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1, ... %%quitte: 212 <-> mc: 216, class rep: 23
    5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 213 <-> mc: 217, class rep: 7
    3, 2, 10, 3, 10, 5, 3, 5, 8, 4, 8, 5, 0, 9, 1, -1, ... %%quitte: 214 <-> mc: 218, class rep: 22
    5, 2, 10, 5, 4, 2, 1, 2, 9, 9, 2, 4, -1, -1, -1, -1, ... %%quitte: 215 <-> mc: 219, class rep: 6
    2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1, ... %%quitte: 216 <-> mc: 212, class rep: 29
    0, 11, 4, 0, 3, 11, 4, 11, 5, 2, 1, 11, 5, 11, 1, -1, ... %%quitte: 217 <-> mc: 213, class rep: 25
    0, 5, 2, 0, 9, 5, 2, 5, 11, 4, 8, 5, 11, 5, 8, -1, ... %%quitte: 218 <-> mc: 214, class rep: 25
    9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 219 <-> mc: 215, class rep: 24
    8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 220 <-> mc: 220, class rep: 7
    0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 221 <-> mc: 221, class rep: 3
    8, 5, 4, 8, 3, 5, 9, 5, 0, 0, 5, 3, -1, -1, -1, -1, ... %%quitte: 222 <-> mc: 222, class rep: 6
    9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 223 <-> mc: 223, class rep: 1
    4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 224 <-> mc: 224, class rep: 7
    0, 3, 8, 4, 7, 9, 9, 7, 11, 9, 11, 10, -1, -1, -1, -1, ... %%quitte: 225 <-> mc: 225, class rep: 30
    1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1, ... %%quitte: 226 <-> mc: 226, class rep: 29
    3, 4, 1, 3, 8, 4, 1, 4, 10, 7, 11, 4, 10, 4, 11, -1, ... %%quitte: 227 <-> mc: 227, class rep: 25
    2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1, ... %%quitte: 228 <-> mc: 232, class rep: 27
    9, 7, 10, 9, 4, 7, 10, 7, 2, 8, 0, 7, 2, 7, 0, -1, ... %%quitte: 229 <-> mc: 233, class rep: 25
    3, 10, 7, 3, 2, 10, 7, 10, 4, 1, 0, 10, 4, 10, 0, -1, ... %%quitte: 230 <-> mc: 234, class rep: 25
    1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 231 <-> mc: 235, class rep: 24
    4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1, ... %%quitte: 232 <-> mc: 228, class rep: 23
    9, 4, 7, 9, 7, 11, 9, 11, 1, 2, 1, 11, 0, 3, 8, -1, ... %%quitte: 233 <-> mc: 229, class rep: 22
    11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 234 <-> mc: 230, class rep: 7
    11, 4, 7, 11, 2, 4, 8, 4, 3, 3, 4, 2, -1, -1, -1, -1, ... %%quitte: 235 <-> mc: 231, class rep: 6
    4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 236 <-> mc: 236, class rep: 7
    4, 1, 9, 4, 7, 1, 0, 1, 8, 8, 1, 7, -1, -1, -1, -1, ... %%quitte: 237 <-> mc: 237, class rep: 6
    4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 238 <-> mc: 238, class rep: 3
    4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 239 <-> mc: 239, class rep: 1
    9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 240 <-> mc: 240, class rep: 15
    3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 241 <-> mc: 241, class rep: 7
    0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 242 <-> mc: 242, class rep: 7
    3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 243 <-> mc: 243, class rep: 3
    2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 244 <-> mc: 248, class rep: 7
    9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 245 <-> mc: 249, class rep: 3
    2, 8, 3, 2, 10, 8, 0, 8, 1, 1, 8, 10, -1, -1, -1, -1, ... %%quitte: 246 <-> mc: 250, class rep: 6
    1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 247 <-> mc: 251, class rep: 1
    1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 248 <-> mc: 244, class rep: 7
    3, 9, 0, 3, 11, 9, 1, 9, 2, 2, 9, 11, -1, -1, -1, -1, ... %%quitte: 249 <-> mc: 245, class rep: 6
    0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 250 <-> mc: 246, class rep: 3
    3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 251 <-> mc: 247, class rep: 1
    1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 252 <-> mc: 252, class rep: 3
    0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 253 <-> mc: 253, class rep: 1
    0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, ... %%quitte: 254 <-> mc: 254, class rep: 1
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1  %%quitte: 255 <-> mc: 255, class rep: 0
    ]);

    % clear plot
    clearMCTriangles(cstate_);

    c_pos = 2;
    color  = getColor(c_pos);
    
    % edge vertices table
    l_edges = [[1,2]; [2,4]; [3,4]; [1,3]; [5,6]; [6,8]; [7,8]; [5,7]; [1,5]; [2,6]; [4,8]; [3,7]];
    vertices = [[0,0,0]; [1,0,0]; [0,1,0]; [1,1,0]; [0,0,1]; [1,0,1]; [0,1,1]; [1,1,1]];
    % implement marching cubes
    % loop over the edges
    mc_bits = [1,2,4,8,16,32,64,128,256,512,1024,2048];
    vt = zeros(12,3);
    % compute case
    mc = 0;
    for i = 1:8
        if f(i) <= i0
            mc = bitset(mc,i);
        end
    end
    for e = 1:12
        if bitand(e_pattern(mc+1),mc_bits(e))
            v1 = l_edges(e,1);
            v2 = l_edges(e,2);
            l = (i0 - f(v1))/(f(v2) - f(v1));
            vt(e,:) = vertices(v1,:) + l*(vertices(v2,:) - vertices(v1,:));
        end
    end % loop over the edges
    % connect edges to build triangles
    % there might be a maximum of 5 triangles
    t_handles = [];
    tris = zeros(5,3);
    t_count = 0;
    for t = 1:3:16
        t_index = 16 * mc + t;
        if t_pattern(t_index) == -1
            break;
        end
        t_count = t_count + 1;
        % t_pattern start with index 0!
        e1 = t_pattern(t_index)   + 1;
        e2 = t_pattern(t_index+1) + 1;
        e3 = t_pattern(t_index+2) + 1;
        tris(t_count,:) = [e1,e2,e3];
    end % loop over triangles for this mc case

    % collect and plot
    if t_count > 0
        t_indices = zeros(t_count,3);
        v_coords  = zeros(t_count,3);
        for t = 1:t_count
            t_indices(t,:) = tris(t,:);
        end
        %t_hdl = patch('Faces',t_indices,'Vertices',vt,'FaceColor',[0.8,0.8,0.8]); %,'EdgeColor','b');
        t_hdl = patch('Faces',t_indices,'Vertices',vt,'FaceColor',color); %,'EdgeColor','b');
        t_handles = [t_handles,t_hdl];
    end
    
    % process second cell
    c_pos = c_pos + 1;
    color = getColor(c_pos);
    mc = 0;
    for i = 1:8
        if f1(i) <= i0
            mc = bitset(mc,i);
        end
    end
    % shift vertices
    coordsh = [0,0,0];
    if coord == 1
        coordsh = [1,0,0];
    elseif coord == 2
        coordsh = [0,1,0];
    elseif coord == 3
        coordsh = [0,0,1];
    end
    vt = zeros(12,3);
    for e = 1:12
        if bitand(e_pattern(mc+1),mc_bits(e))
            v1 = l_edges(e,1);
            v2 = l_edges(e,2);
            l = (i0 - f1(v1))/(f1(v2) - f1(v1));
            vt(e,:) = vertices(v1,:) + l*(vertices(v2,:) - vertices(v1,:));
            vt(e,:) = vt(e,:) + coordsh;
        end
    end % loop over the edges
    % connect edges to build triangles
    % there might be a maximum of 5 triangles
    tris = zeros(5,3);
    t_count = 0;
    for t = 1:3:16
        t_index = 16 * mc + t;
        if t_pattern(t_index) == -1
            break;
        end
        t_count = t_count + 1;
        % t_pattern start with index 0!
        e1 = t_pattern(t_index)   + 1;
        e2 = t_pattern(t_index+1) + 1;
        e3 = t_pattern(t_index+2) + 1;
        tris(t_count,:) = [e1,e2,e3];
    end % loop over triangles for this mc case

    % collect and plot
    if t_count > 0
        t_indices = zeros(t_count,3);
        v_coords  = zeros(t_count,3);
        for t = 1:t_count
            t_indices(t,:) = tris(t,:);
        end
        %t_hdl = patch('Faces',t_indices,'Vertices',vt,'FaceColor',[0.8,0.8,0.8]); %,'EdgeColor','b');
        t_hdl = patch('Faces',t_indices,'Vertices',vt,'FaceColor',color); %,'EdgeColor','b');
        t_handles = [t_handles,t_hdl];
    end
    
    % save handler to graphics object
    cstate_('MCTrisHandles') = t_handles;
    for t = 1:length(t_handles)
        set(t_handles(t),'FaceLighting','none');
    end
end


%
% Plot Asymptotes for asymptotic decider
%
function plotAsymptotes(cstate_,coord,f,f1,mc,i0)
    % keep pointer to plots
    l_handles = [];
    % process all six faces
    face = 1;
    [u0,v0,idx] = getAsymptotes(face,f);
    if (0 < v0 && v0 < 1)
        uc = [0,1];
        vc = [v0,v0];
        wc = [0,0];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < u0 && u0 < 1)
        uc = [u0,u0];
        vc = [0,1];
        wc = [0,0];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 2;
    [u0,v0,idx] = getAsymptotes(face,f);
    if (0 < u0 && u0 < 1)
       uc = [u0,u0];
        vc = [0,1];
        wc = [1,1];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < v0 && v0 < 1)
        uc = [0,1];
        vc = [v0,v0];
        wc = [1,1];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 3;
    [u0,w0,idx] = getAsymptotes(face,f);
    if (0 < u0 && u0 < 1)
        uc = [u0,u0];
        vc = [0,0];
        wc = [0,1];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [0,1];
        vc = [0,0];
        wc = [w0,w0];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 4;
    [u0,w0,idx] = getAsymptotes(face,f);
    if (0 < u0 && u0 < 1)
        uc = [u0,u0];
        vc = [1,1];
        wc = [0,1];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [0,1];
        vc = [1,1];
        wc = [w0,w0];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 5;
    [v0,w0,idx] = getAsymptotes(face,f);
    if (0 < v0 && v0 < 1)
        uc = [0,0];
        vc = [v0,v0];
        wc = [0,1];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [0,0];
        vc = [0,1];
        wc = [w0,w0];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 6;
    [v0,w0,idx] = getAsymptotes(face,f);
    if (0 < v0 && v0 < 1)
        uc = [1,1];
        vc = [v0,v0];
        wc = [0,1];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [1,1];
        vc = [0,1];
        wc = [w0,w0];
        l_ = plot3(uc,vc,wc,'c','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    
    % plot second cell
    face = 1;
    xsh = 0;
    ysh = 0;
    zsh = 0;
    if coord == 1
        xsh = 1;
    elseif coord == 2
        ysh = 1;
    elseif coord == 3
        zsh = 1;
    end
    [u0,v0,idx] = getAsymptotes(face,f1);
    if (0 < v0 && v0 < 1)
        uc = [0,1] + xsh;
        vc = [v0,v0] + ysh;
        wc = [0,0] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < u0 && u0 < 1)
        uc = [u0,u0] + xsh;
        vc = [0,1] + ysh;
        wc = [0,0] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 2;
    [u0,v0,idx] = getAsymptotes(face,f1);
    if (0 < u0 && u0 < 1)
       uc = [u0,u0] + xsh;
        vc = [0,1] + ysh;
        wc = [1,1] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < v0 && v0 < 1)
        uc = [0,1] + xsh;
        vc = [v0,v0] + ysh;
        wc = [1,1] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 3;
    [u0,w0,idx] = getAsymptotes(face,f1);
    if (0 < u0 && u0 < 1)
        uc = [u0,u0] + xsh;
        vc = [0,0] + ysh;
        wc = [0,1] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [0,1] + xsh;
        vc = [0,0] + ysh; 
        wc = [w0,w0] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 4;
    [u0,w0,idx] = getAsymptotes(face,f1);
    if (0 < u0 && u0 < 1)
        uc = [u0,u0] + xsh;
        vc = [1,1] + ysh;
        wc = [0,1] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [0,1] + xsh;
        vc = [1,1] + ysh;
        wc = [w0,w0] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 5;
    [v0,w0,idx] = getAsymptotes(face,f1);
    if (0 < v0 && v0 < 1)
        uc = [0,0] + xsh;
        vc = [v0,v0] + ysh;
        wc = [0,1] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [0,0] + xsh;
        vc = [0,1] + ysh;
        wc = [w0,w0] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    %%
    face = 6;
    [v0,w0,idx] = getAsymptotes(face,f1);
    if (0 < v0 && v0 < 1)
        uc = [1,1] + xsh;
        vc = [v0,v0] + ysh;
        wc = [0,1] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    if (0 < w0 && w0 < 1)
        uc = [1,1] + xsh;
        vc = [0,1] + ysh;
        wc = [w0,w0] + zsh;
        l_ = plot3(uc,vc,wc,'y','LineWidth',2);
        l_handles = [l_handles,l_];
    end
    
    % save pointers to plots in global variable
    cstate_('Asymptotes') = l_handles;
    
end


function plotHyperbolicArcs(cstate_,f,mc,i0)
    % u handles
    u_handles = [];
    v_handles = [];
    w_handles = [];
    % there are two opposite faces
    % face 5, 6
    f_color = {'m','c'};
    u_face = [0,1];
    v_face = [0,1];
    w_face = [0,1];
    for fnr = 1:2
        [vu,wu,eu,au] = uSurface(f,u_face(fnr));
        sz_v = 500;
        sal = 0.0000000000000001;
        dvl = (vu-sal) / (sz_v - 1);
        vl = 0;
        dvr = (1-vu-sal)/(sz_v - 1);
        vr = vu + sal;
        vcl = [];
        wcl = [];
        ucl = [];
        ucr = [];
        vcr = [];
        wcr = [];
        w1 = 0;
        v1 = vu + (i0-au)/(eu*(w1-wu));
        w2 = 1;
        v2 = vu + (i0-au)/(eu*(w2-wu));
        if v2 < vu
            tmp = v2;
            v2 = v1;
            v1 = tmp;
            tmp = w2;
            w2 = w1;
            w1 = tmp;
        end
        c_count = 1;
        for i=1:sz_v
            w = wu + (i0 - au) / (eu*(vl - vu));
            if w >= 0 && w <= 1
                vcl(c_count) = vl;
                wcl(c_count) = w;
                ucl(c_count) = u_face(fnr);
                c_count = c_count + 1;
            end
            vl = vl + dvl;
        end
        vcl(c_count) = v1;
        wcl(c_count) = w1;
        ucl(c_count) = u_face(fnr);
        l_ = plot3(ucl,vcl,wcl,f_color{fnr},'LineWidth',2);
        set(l_,'Visible','off');
        u_handles = [u_handles,l_];

        c_count = 1;
        vcr(c_count) = v2;
        wcr(c_count) = w2;
        ucr(c_count) = u_face(fnr);
        c_count = c_count + 1;
        for i=1:sz_v
            w = wu + (i0 - au) / (eu*(vr-vu));
            if w >= 0 && w <= 1
                vcr(c_count) = vr;
                wcr(c_count) = w;
                ucr(c_count) = u_face(fnr);
                c_count = c_count + 1;
            end
            vr = vr + dvr;
        end
        l_ = plot3(ucr,vcr,wcr,f_color{fnr},'LineWidth',2);
        set(l_,'Visible','off');
        u_handles = [u_handles,l_];
    end % face 5, 6

    % face 3,4
    for fnr = 1:2
        [uv,wv,ev,av] = vSurface(f,v_face(fnr));
        sz_u = 500;
        sal = 0.0000000000000001;
        dul = (uv-sal) / (sz_u - 1);
        ul = 0;
        dur = (1-uv-sal)/(sz_u - 1);
        ur = uv + sal;
        ucl = [];
        vcl = [];
        wcl = [];
        ucr = [];
        vcr = [];
        wcr = [];
        w1 = 0;
        u1 = uv + (i0-av)/(ev*(w1-wv));
        w2 = 1;
        u2 = uv + (i0-av)/(ev*(w2-wv));
        if u2 < uv
            tmp = u2;
            u2 = u1;
            u1 = tmp;
            tmp = w2;
            w2 = w1;
            w1 = tmp;
        end
        c_count = 1;
        for i=1:sz_u
            w = wv + (i0 - av) / (ev*(ul - uv));
            if w >= 0 && w <= 1
                ucl(c_count) = ul;
                vcl(c_count) = v_face(fnr);
                wcl(c_count) = w;
                c_count = c_count + 1;
            end
            ul = ul + dul;
        end
        ucl(c_count) = u1;
        vcl(c_count) = v_face(fnr);
        wcl(c_count) = w1;

        l_ = plot3(ucl,vcl,wcl,f_color{fnr},'LineWidth',2);
        set(l_,'Visible','off');
        v_handles = [v_handles,l_];

        c_count = 1;
        ucr(c_count) = u2
        vcr(c_count) = v_face(fnr);
        wcr(c_count) = w2;

        c_count = c_count + 1;
        for i=1:sz_u
            w = wv + (i0 - av) / (ev*(ur-uv));
            if w >= 0 && w <= 1
                ucr(c_count) = ur;
                vcr(c_count) = v_face(fnr);
                wcr(c_count) = w;

                c_count = c_count + 1;
            end
            ur = ur + dur;
        end
        l_ = plot3(ucr,vcr,wcr,f_color{fnr},'LineWidth',2);
        set(l_,'Visible','off');
        v_handles = [v_handles,l_];
    end % face 5, 6

    % face 1, 2
    for fnr = 1:2
        [uw,vw,ew,aw] = wSurface(f,w_face(fnr));
        sz_u = 500;
        sal = 0.0000000000000001;
        dul = (uw-sal) / (sz_u - 1);
        ul = 0;
        dur = (1-uw-sal)/(sz_u - 1);
        ur = uw + sal;
        vcl = [];
        wcl = [];
        ucl = [];
        ucr = [];
        vcr = [];
        wcr = [];
        v1 = 0;
        u1 = uw + (i0-aw)/(ew*(v1-vw));
        v2 = 1;
        u2 = uw + (i0-aw)/(ew*(v2-vw));
        if u2 < uw
            tmp = u2;
            u2 = u1;
            u1 = tmp;
            tmp = v2;
            v2 = v1;
            v1 = tmp;
        end
        c_count = 1;
        for i=1:sz_u
            v = vw + (i0 - aw) / (ew*(ul - uw));
            if v >= 0 && v <= 1
                ucl(c_count) = ul;
                vcl(c_count) = v;
                wcl(c_count) = w_face(fnr);
                c_count = c_count + 1;
            end
            ul = ul + dul;
        end
        ucl(c_count) = u1;
        vcl(c_count) = v1;
        wcl(c_count) = w_face(fnr);
        l_ = plot3(ucl,vcl,wcl,f_color{fnr},'LineWidth',2);
        set(l_,'Visible','off');
        w_handles = [w_handles,l_];

        c_count = 1;
        ucr(c_count) = u2;
        vcr(c_count) = v2;
        wcr(c_count) = w_face(fnr);
        c_count = c_count + 1;
        for i=1:sz_u
            v = vw + (i0 - aw) / (ew*(ur-uw));
            if v >= 0 && v <= 1
                ucr(c_count) = ur;
                vcr(c_count) = v;
                wcr(c_count) = w_face(fnr);
                c_count = c_count + 1;
            end
            ur = ur + dur;
        end
        l_ = plot3(ucr,vcr,wcr,f_color{fnr},'LineWidth',2);
        set(l_,'Visible','off');
        w_handles = [w_handles,l_];
    end % face 5, 6

    % cstate_
    cstate_('HyperbolicArcsU') = u_handles;
    cstate_('HyperbolicArcsV') = v_handles;
    cstate_('HyperbolicArcsW') = w_handles;

end

function setHyperbolicArcsVisibility(cstate_, vis_flag, face)

    switch face
        case 'u_face'
            u_handles = cstate_('HyperbolicArcsU');
            for uh = 1:length(u_handles)
                if vis_flag
                    set(u_handles(uh),'Visible','on');
                else
                    set(u_handles(uh),'Visible','off');
                end
            end
        case 'v_face'
            v_handles = cstate_('HyperbolicArcsV');
            for vh = 1:length(v_handles)
                if vis_flag
                    set(v_handles(vh),'Visible','on');
                else
                    set(v_handles(vh),'Visible','off');
                end
            end
        case 'w_face'
            w_handles = cstate_('HyperbolicArcsW');
            for wh = 1:length(w_handles)
                if vis_flag
                    set(w_handles(wh),'Visible','on');
                else
                    set(w_handles(wh),'Visible','off');
                end
            end
    end
end


function plotInnerHexagon(cstate_,F,mc,i0)
    [ic,cvt,fs,ui,vi,wi] = quadraticEquationsAdvanced(F,i0);
    ict = ic(1) + ic(2) + ic(3);
    if ict == 6
        lt_handles = [];
        uc = [cvt(1:6,1);cvt(1,1)];
        vc = [cvt(1:6,2);cvt(1,2)];
        wc = [cvt(1:6,3);cvt(1,3)];
        l_ = plot3(uc,vc,wc,'black','LineWidth',3);
        lt_handles = [lt_handles,l_];
        cstate_('InnerHexaHandles') = lt_handles;
    else
        cstate_('InnerHexaHandles') = [];
    end
end


% compute vertex at face for singular case
function fvt = getFaceVertex(face,f)
    [u0,v0,idx] = getAsymptotes(face,f);
    switch(face)
        case 1
            fvt = [u0,v0,0];
        case 2
            fvt = [u0,v0,1];
        case 3
            fvt = [u0,0,v0];
        case 4
            fvt = [u0,1,v0];
        case 5
            fvt = [0,u0,v0];
        case 6
            fvt = [1,u0,v0];
    end
end

% Compute contours
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [nrv,vt,nrc,cns,contours,face, fvt,s_flag, a_cnt] = computeOrientedContours(f,mc,isov)
    % Use the fact that the isosurface does not start or end in the cell.
    % It must cut the cell at the edges building closed hyperbolic arcs.
    % Triangulatiom
    %   i) compute intersection of level set with edges
    %  ii) compute closed polygons, we call these polygon contours
    % iii) compute the solution to the three quadratic equations
    %  iv) classify and triangulate contours

    % cell topology
    % vertices accroding to the coordinate axes, x-axes fastes
    % v1 = 000, v2 = 001, v3 = 010, v4 = 011
    % v5 = 100, v6 = 101, v7 = 110, v8 = 111
    % edges:
    % e1 = (1,2), e2 = (2,4), e3 = (3,4), e4 = (1,3)
    % e5 = (5,6), e6 = (6,8), e7 = (7,8), e8 = (5,7)
    % e9 = (1,5), e10 = (2,6), e11 = (4,8), e12 = (3,7)
    % faces:
    % f1 = (1,2,3,4), f2 = (5,6,7,8),
    % f3 = (1,2,5,6), f4 = (3,4,7,8),
    % f5 = (1,3,5,7), f6 = (2,4,6,8)

    % compute offsets along axis,
    % these are local coordinates used for interpolating vertices later on
    edges = [[1 2 -1  0  0]; ...
        [2 4  1 -1  0]; ...
        [3 4 -1  1  0]; ...
        [1 3  0 -1  0]; ...
        [5 6 -1  0  1]; ...
        [6 8  1 -1  1]; ...
        [7 8 -1  1  1]; ...
        [5 7  0 -1  1]; ...
        [1 5  0  0 -1]; ...
        [2 6  1  0 -1]; ...
        [4 8  1  1 -1]; ...
        [3 7  0  1 -1]];
    
    fvt = -ones(6,3);
    cnt_flag = false(1,4);
    vt = zeros(12,4); % 4. coordinate means if vertex is set
    nrv = 0;
    for i = 1:12
        v1 = edges(i,1);
        v2 = edges(i,2);
        vals = [f(v1),f(v2)];
        vals = sort(vals);
        if vals(1) < isov && isov <= vals(2)
            % there is an intersection
            nrv = nrv + 1;
            iv = (isov - f(v1)) / (f(v2) - f(v1));
            if (edges(i,3) < 0)
                vt(i,1) = iv;
            else
                vt(i,1) = edges(i,3);
            end
            if (edges(i,4) < 0)
                vt(i,2) = iv;
            else
                vt(i,2) = edges(i,4);
            end
            if (edges(i,5) < 0)
                vt(i,3) = iv;
            else
                vt(i,3) = edges(i,5);
            end
            vt(i,4) = 1;
        else
            vt(i,4) = -1;
        end
    end
    % build contours, porcess facewise
    % faces
    fedge = [[1 2 3 4];[5 6 7 8];[1 10 5 9];[3 11 7 12];[4 12 8 9];[2 11 6 10]];
    % code face vertices
   	face_v = [[1,3,2,4];[6,8,5,7];[5,1,6,2];[3,7,4,8];[5,7,1,3];[2,4,6,8]];
    face_e = [[4,3,2,1];[6,7,8,5];[9,1,10,5];[12,7,11,3];[8,12,4,9];[2,11,6,10]];
    segments = zeros(12,2);
    face_index = zeros(1,12);
    nrs = 1;
    ambiguous = zeros(1,6);
    ambcase   = zeros(1,6);
    s_amb = false(1,12);
    for i = 1:6
        f_case = 0;
        if f(face_v(i,1)) >= isov
            f_case = bitset(f_case,1);
        end
        if f(face_v(i,2)) >= isov
            f_case = bitset(f_case,2);
        end
        if f(face_v(i,3)) >= isov
            f_case = bitset(f_case,3);
        end
        if f(face_v(i,4)) >= isov
            f_case = bitset(f_case,4);
        end
        % triangulate face
        e1 = face_e(i,1);
        e2 = face_e(i,2);
        e3 = face_e(i,3);
        e4 = face_e(i,4);
        f1 = f(face_v(i,1));
        f2 = f(face_v(i,2));
        f3 = f(face_v(i,3));
        f4 = f(face_v(i,4));
        switch f_case
            case 1
                segments(e1,1) = e4; % outgoing to segment e4
                segments(e4,2) = e1; % incoming from segment e1
            case 2
                segments(e2,1) = e1;
                segments(e1,2) = e2;
            case 3
                segments(e2,1) = e4;
                segments(e4,2) = e2;
            case 4
                segments(e4,1) = e3;
                segments(e3,2) = e4;
            case 5
                segments(e1,1) = e3;
                segments(e3,2) = e1;
            case 6
                alp = aDecider(f1,f2,f3,f4);
                small = 0.00000000000001;
                if abs(alp - isov) < small
                    alp = isov;
                end
                if alp >= isov
                    segments(e4,1) = e1;
                    segments(e1,2) = e4;
                    segments(e2,1) = e3;
                    segments(e3,2) = e2;
                    if alp == isov
                        s_amb(e4) = true;
                        s_amb(e2) = true;
                        fvt(i,:) = getFaceVertex(i,f);
                        face_index(e4) = i;
                        face_index(e2) = i; %for this edge the face is ambiguous and singular
                    end
                else
                    segments(e4,1) = e3;
                    segments(e3,2) = e4;
                    segments(e2,1) = e1;
                    segments(e1,2) = e2;
                end
            case 7
                segments(e2,1) = e3;
                segments(e3,2) = e2;
            case 8
                segments(e3,1) = e2;
                segments(e2,2) = e3;
            case 9
                alp = aDecider(f1,f2,f3,f4);
                small = 0.00000000000001;
                if abs(alp - isov) < small
                    alp = isov;
                end
                if alp >= isov
                    segments(e1,1) = e2;
                    segments(e2,2) = e1;
                    segments(e3,1) = e4;
                    segments(e4,2) = e3;
                    if alp == isov
                        s_amb(e1) = true;
                        s_amb(e3) = true;
                        fvt(i,:) = getFaceVertex(i,f);
                        face_index(e1) = i;
                        face_index(e3) = i;
                    end
                else
                    segments(e1,1) = e4;
                    segments(e4,2) = e1;
                    segments(e3,1) = e2;
                    segments(e2,2) = e3;
                end
            case 10
                segments(e3,1) = e1;
                segments(e1,2) = e3;
            case 11
                segments(e3,1) = e4;
                segments(e4,2) = e3;
            case 12
                segments(e4,1) = e2;
                segments(e2,2) = e4;
            case 13
                segments(e1,1) = e2;
                segments(e2,2) = e1;
            case 14
                segments(e4,1) = e1;
                segments(e1,2) = e4;
        end
    end % loop over faces

    % compute closed contours
    [contours, s_flag,a_cnt,face] = connectOrietedContours(nrv,vt,s_amb,face_index,segments);
    % count nr. of countours and size
    cns = zeros(1,4);
    nrc = 0;
    for i = 1:4
        pos = 0;
        for j = 1:12
            if contours(i,j) > 0
                pos = j;
            end
        end
        cns(i) = pos;
        if pos > 0
            nrc = nrc+1;
        end
    end
end

% connect vertices into a closed contour
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [contours, s_amb,a_cnt,face] = connectOrietedContours(nrv,vt,s_flag,face_idx,segment)
    %
    %
    % there might be at most 4 contours with at most 12 vertices
    tcont = zeros(4,13);
    face = -ones(4,12);
    % each vertex appears at most once in only one contour
    % a contour has at least three vertices, i.e. it is a triangle or a larger
    % polygon
    cnt = 0;
    a_cnt = false(1,4);
    s_amb = false(4,12);
    for e = 1:12
        eTo = segment(e,1);
        eIn = segment(e,2);
        if eTo > 0
            cnt = cnt + 1;
            % connect contour
            eStart = e;
            pos = 1;
            tcont(cnt,pos) = eStart;
            if (s_flag(eStart))
                s_amb(cnt,pos) = true;
                a_cnt(cnt) = true;
                face(cnt,pos) = face_idx(eStart);
            end
            while eTo ~= eStart
                pos = pos + 1;
                tcont(cnt,pos) = eTo;
                if (s_flag(eTo))
                    s_amb(cnt,pos) = true;
                    a_cnt(cnt) = true;
                    face(cnt,pos) = face_idx(eTo);
                end
                eIn = eTo;
                eTo = segment(eIn,1);
                segment(eIn,1) = -1;
            end
        end
    end
    contours = tcont;
end


function [ct,cvt,fs,ui,vi,wi] = quadraticEquations(F,i0)
    % init
    ct = zeros(1,3);
    cvt = zeros(6,3); 
    ui = zeros(1,2);
    vi = zeros(1,2);
    wi = zeros(1,2);
    fs = zeros(3,2);
    % face pair 1,3
    f1 = F(1);
    f2 = F(2);
    f3 = F(3);
    f4 = F(4);
    h1 = F(5);
    h2 = F(6);
    h3 = F(7);
    h4 = F(8);
    
    % solutions
    % quadratic equation
    a = (f1-f2)*(-h3+h4+h1-h2)-(h1-h2)*(-f3+f4+f1-f2);
    b = (i0-f1)*(-h3+h4+h1-h2)+(f1-f2)*(h3-h1)-(i0-h1)*(-f3+f4+f1-f2)-(h1-h2)*(f3-f1);
    c = (i0-f1)*(h3-h1)-(i0-h1)*(f3-f1);
    d = b^2 - 4*a*c;
    if d > 0
        d = sqrt(d);
        % compute ui
        ui1 = (-b-d) / (2*a);
        ui2 = (-b+d) / (2*a);
        g1 = f1*(1-ui1) + f2*ui1;
        g2 = f3*(1-ui1) + f4*ui1;
        vi1 = (i0 - g1)/(g2-g1);
        g1 = f1*(1-ui2) + f2*ui2;
        g2 = f3*(1-ui2) + f4*ui2;  
        vi2 = (i0 - g1)/(g2-g1);
        ui(1) = ui1;
        ui(2) = ui2;
        vi(1) = vi1;
        vi(2) = vi2;
        
        % compute wi
        f1 = F(1);
        f2 = F(2);
        f3 = F(5);
        f4 = F(6);
        g1 = f1*(1-ui1) + f2*ui1;
        g2 = f3*(1-ui1) + f4*ui1;
        %if abs(g1 - g2) > 0.000001
            wi(1) = (i0 - g1)/(g2-g1);
%         else 
%             f1 = F(1);
%             f2 = F(3);
%             f3 = F(5);
%             f4 = F(7);
%             g1 = f1*(1-vi2) + f2*vi2;
%             g2 = f3*(1-vi2) + f4*vi2;
%             wi(1) = (i0 - g1)/(g2-g1);
%         end
        
        f1 = F(1);
        f2 = F(2);
        f3 = F(5);
        f4 = F(6);
        g1 = f1*(1-ui2) + f2*ui2;
        g2 = f3*(1-ui2) + f4*ui2;
        wi(2) = (i0 - g1)/(g2-g1);
        %if abs(g1 - g2) > 0.00000000001
            wi(2) = (i0 - g1)/(g2-g1);
%         else 
%             f1 = F(3);
%             f2 = F(1);
%             f3 = F(7);
%             f4 = F(5);
%             g1 = f1*(1-vi1) + f2*vi1;
%             g2 = f3*(1-vi1) + f4*vi1;
%             wi(2) = (i0 - g1)/(g2-g1);
%         end
        
        % numerical error
        small = 0.00000000000001;
        for i = 1:2
            if ui(i) < 0 && abs(ui(i)) < small
                ui(i) = 0;
            end
            if ui(i) > 1 && abs(ui(i) - 1) < small
                ui(i) = 1;
            end
            if vi(i) < 0 && abs(vi(i)) < small
                vi(i) = 0;
            end
            if vi(i) > 1 && abs(vi(i) - 1) < small
                vi(i) = 1;
            end
            if wi(i) < 0 && abs(wi(i)) < small
                wi(i) = 0;
            end
            if wi(i) > 1 && abs(wi(i) - 1) < small
                wi(i) = 1;
            end
        end
        
        % create 6 points
        cvt(1,:) = [ui(1),vi(1),wi(1)];
        cvt(2,:) = [ui(1),vi(1),wi(2)];
        cvt(3,:) = [ui(2),vi(1),wi(2)];
        cvt(4,:) = [ui(2),vi(2),wi(2)];
        cvt(5,:) = [ui(2),vi(2),wi(1)];
        cvt(6,:) = [ui(1),vi(2),wi(1)];
        
        face = 1;
        [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
        e1 = f0+f3-f1-f2;
        alpha_face1 = (f0*f3 - f1*f2)/e1;
        e2 = h0+h3-h1-h2;
        alpha_face2 = (h0*h3 - h1*h2)/e2;
        face = 3;
        [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
        e1 = f0+f3-f1-f2;
        alpha_face3 = (f0*f3 - f1*f2)/e1;
        e2 = h0+h3-h1-h2;
        alpha_face4 = (h0*h3 - h1*h2)/e2;
        face = 5;
        [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
        e1 = f0+f3-f1-f2;
        alpha_face5 = (f0*f3 - f1*f2)/e1;
        e2 = h0+h3-h1-h2;
        alpha_face6 = (h0*h3 - h1*h2)/e2;
        
%         if (alpha_face1 == i0) 
%             if wi(1) < wi(2)
%                 wi(1) = 0;
%             else
%                 wi(2) = 0;
%             end
%         end
%         if (alpha_face2 == i0) 
%             if wi(1) > ui(2)
%                 wi(1) = 1;
%             else
%                 wi(2) = 1;
%             end
%         end
%         if (alpha_face3 == i0) 
%             if vi(1) < vi(2)
%                 vi(1) = 0;
%             else
%                 vi(2) = 0;
%             end
%         end
%         if (alpha_face4 == i0) 
%             if vi(1) > vi(2)
%                 vi(1) = 1;
%             else
%                 vi(2) = 1;
%             end
%         end
%         if (alpha_face5 == i0) 
%             if ui(1) < ui(2)
%                 ui(1) = 0;
%             else
%                 ui(2) = 0;
%             end
%         end
%         if (alpha_face6 == i0) 
%             if ui(1) > ui(2)
%                 ui(1) = 1;
%             else
%                 ui(2) = 1;
%             end
%         end
        
        % count nr. of solutions
        if (ui(1) >= 0 && ui(1) <= 1) && (vi(1) >= 0 && vi(1) <= 1)
            ct(1) = ct(1) + 1;
        end
        if (ui(2) >= 0 && ui(2) <= 1) && (vi(2) >= 0 && vi(2) <= 1)
            ct(1) = ct(1) + 1;
        end
        if (ui(1) >= 0 && ui(1) <= 1) && (wi(1) >= 0 && wi(1) <= 1)
            ct(2) = ct(2) + 1;
        end
        if (ui(2) >= 0 && ui(2) <= 1) && (wi(2) >= 0 && wi(2) <= 1)
            ct(2) = ct(2) + 1;
        end
        if (vi(1) >= 0 && vi(1) <= 1) && (wi(2) >= 0 && wi(2) <= 1)
            ct(3) = ct(3) + 1;
        end
        if (vi(2) >= 0 && vi(2) <= 1) && (wi(1) >= 0 && wi(1) <= 1)
            ct(3) = ct(3) + 1;
        end
        if ui(1) >= 0 && ui(1) <= 1
            fs(1,1) = 1;
        end
        if ui(2) >= 0 && ui(2) <= 1
            fs(1,2) = 1;
        end
        if vi(1) >= 0 && vi(1) <= 1
            fs(2,1) = 1;
        end
        if vi(2) >= 0 && vi(2) <= 1
            fs(2,2) = 1;
        end
        if wi(1) >= 0 && wi(1) <= 1
            fs(3,1) = 1;
        end
        if wi(2) >= 0 && wi(2) <= 1
            fs(3,2) = 1;
        end
    else
        disp('complex solution');
        %l_ = [];
    end
    % debugging purpose
    %disp(['ui1 = ',num2str(ui(1)), ', ui2 = ',num2str(ui(2))]);
    %disp(['vi1 = ',num2str(vi(1)), ', vi2 = ',num2str(vi(2))]);
    %disp(['wi1 = ',num2str(wi(1)), ', wi2 = ',num2str(wi(2))]);
end


function [ct,cvt,fs,ui,vi,wi] = quadraticEquationsSingular(F,i0)
    % init
    ct = zeros(1,3);
    cvt = zeros(6,3); 
    ui = zeros(1,2);
    vi = zeros(1,2);
    wi = zeros(1,2);
    fs = zeros(3,2);
    % face pair 1,3
    f1 = F(1);
    f2 = F(2);
    f3 = F(3);
    f4 = F(4);
    h1 = F(5);
    h2 = F(6);
    h3 = F(7);
    h4 = F(8);
    
    % solutions
    % quadratic equation
    a = (f1-f2)*(-h3+h4+h1-h2)-(h1-h2)*(-f3+f4+f1-f2);
    b = (i0-f1)*(-h3+h4+h1-h2)+(f1-f2)*(h3-h1)-(i0-h1)*(-f3+f4+f1-f2)-(h1-h2)*(f3-f1);
    c = (i0-f1)*(h3-h1)-(i0-h1)*(f3-f1);
    d = b^2 - 4*a*c;
    if d > 0
        d = sqrt(d);
        % compute ui
        ui1 = (-b-d) / (2*a);
        ui2 = (-b+d) / (2*a);
        g1 = f1*(1-ui1) + f2*ui1;
        g2 = f3*(1-ui1) + f4*ui1;
        vi1 = (i0 - g1)/(g2-g1);
        g1 = f1*(1-ui2) + f2*ui2;
        g2 = f3*(1-ui2) + f4*ui2;  
        vi2 = (i0 - g1)/(g2-g1);
        ui(1) = ui1;
        ui(2) = ui2;
        vi(1) = vi1;
        vi(2) = vi2;
        
        % compute wi
        f1 = F(1);
        f2 = F(2);
        f3 = F(5);
        f4 = F(6);
        g1 = f1*(1-ui1) + f2*ui1;
        g2 = f3*(1-ui1) + f4*ui1;
        if abs(g1 - g2) > 0.000001
            wi(1) = (i0 - g1)/(g2-g1);
        else 
            f1 = F(1);
            f2 = F(3);
            f3 = F(5);
            f4 = F(7);
            g1 = f1*(1-vi2) + f2*vi2;
            g2 = f3*(1-vi2) + f4*vi2;
            wi(1) = (i0 - g1)/(g2-g1);
        end
        
        f1 = F(1);
        f2 = F(2);
        f3 = F(5);
        f4 = F(6);
        g1 = f1*(1-ui2) + f2*ui2;
        g2 = f3*(1-ui2) + f4*ui2;
        if abs(g1 - g2) > 0.000001
            wi(2) = (i0 - g1)/(g2-g1);
        else 
            f1 = F(1);
            f2 = F(3);
            f3 = F(5);
            f4 = F(7);
            g1 = f1*(1-vi1) + f2*vi1;
            g2 = f3*(1-vi1) + f4*vi1;
            wi(2) = (i0 - g1)/(g2-g1);
        end
        
        % create 6 points
        cvt(1,:) = [ui(1),vi(1),wi(1)];
        cvt(2,:) = [ui(1),vi(1),wi(2)];
        cvt(3,:) = [ui(2),vi(1),wi(2)];
        cvt(4,:) = [ui(2),vi(2),wi(2)];
        cvt(5,:) = [ui(2),vi(2),wi(1)];
        cvt(6,:) = [ui(1),vi(2),wi(1)];
        
        face = 1;
        [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
        e1 = f0+f3-f1-f2;
        alpha_face1 = (f0*f3 - f1*f2)/e1;
        e2 = h0+h3-h1-h2;
        alpha_face2 = (h0*h3 - h1*h2)/e2;
        face = 3;
        [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
        e1 = f0+f3-f1-f2;
        alpha_face3 = (f0*f3 - f1*f2)/e1;
        e2 = h0+h3-h1-h2;
        alpha_face4 = (h0*h3 - h1*h2)/e2;
        face = 5;
        [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
        e1 = f0+f3-f1-f2;
        alpha_face5 = (f0*f3 - f1*f2)/e1;
        e2 = h0+h3-h1-h2;
        alpha_face6 = (h0*h3 - h1*h2)/e2;
        
        if (alpha_face1 == i0) 
            if wi(1) < wi(2)
                wi(1) = 0;
            else
                wi(2) = 0;
            end
        end
        if (alpha_face2 == i0) 
            if wi(1) > ui(2)
                wi(1) = 1;
            else
                wi(2) = 1;
            end
        end
        if (alpha_face3 == i0) 
            if vi(1) < vi(2)
                vi(1) = 0;
            else
                vi(2) = 0;
            end
        end
        if (alpha_face4 == i0) 
            if vi(1) > vi(2)
                vi(1) = 1;
            else
                vi(2) = 1;
            end
        end
        if (alpha_face5 == i0) 
            if ui(1) < ui(2)
                ui(1) = 0;
            else
                ui(2) = 0;
            end
        end
        if (alpha_face6 == i0) 
            if ui(1) > ui(2)
                ui(1) = 1;
            else
                ui(2) = 1;
            end
        end
        
        % count nr. of solutions
        if (ui(1) > 0 && ui(1) < 1) && (vi(1) > 0 && vi(1) < 1)
            ct(1) = ct(1) + 1;
        end
        if (ui(2) > 0 && ui(2) < 1) && (vi(2) > 0 && vi(2) < 1)
            ct(1) = ct(1) + 1;
        end
        if (ui(1) > 0 && ui(1) < 1) && (wi(1) > 0 && wi(1) < 1)
            ct(2) = ct(2) + 1;
        end
        if (ui(2) > 0 && ui(2) < 1) && (wi(2) > 0 && wi(2) < 1)
            ct(2) = ct(2) + 1;
        end
        if (vi(1) > 0 && vi(1) < 1) && (wi(2) > 0 && wi(2) < 1)
            ct(3) = ct(3) + 1;
        end
        if (vi(2) > 0 && vi(2) < 1) && (wi(1) > 0 && wi(1) < 1)
            ct(3) = ct(3) + 1;
        end
        if ui(1) > 0 && ui(1) < 1
            fs(1,1) = 1;
        end
        if ui(2) > 0 && ui(2) < 1
            fs(1,2) = 1;
        end
        if vi(1) > 0 && vi(1) < 1
            fs(2,1) = 1;
        end
        if vi(2) > 0 && vi(2) < 1
            fs(2,2) = 1;
        end
        if wi(1) > 0 && wi(1) < 1
            fs(3,1) = 1;
        end
        if wi(2) > 0 && wi(2) < 1
            fs(3,2) = 1;
        end
    else
        disp('complex solution');
        %l_ = [];
    end
end


% a different version of quadratic equations
function [ct,cvt,fs,ui,vi,wi] = quadraticEquationsAdvanced(F,i0)
    % init
    ct = zeros(1,3);
    cvt = zeros(6,3); 
    ui = zeros(1,2);
    vi = zeros(1,2);
    wi = zeros(1,2);
    fs = zeros(3,2);
    % face pair 1,2, w=0 and w=1
    face = 1;
    [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
    
    e1 = f0+f3-f1-f2;
    v1 = (f0-f1)/e1;
    u1 = (f0-f2)/e1;
    a1 = (f0*f3 - f1*f2)/e1;
    e2 = h0+h3-h1-h2;
    v2 = (h0-h1)/e2;
    u2 = (h0-h2)/e2;
    a2 = (h0*h3 - h1*h2)/e2;
    
    a = (v2-v1)*e1*e2;
    b = -a*(u1+u2) + (i0-a2)*e1 - (i0-a1)*e2;
    c = a*u1*u2 - (i0-a2)*e1*u1 + (i0-a1)*e2*u2;
    d = b*b - 4*a*c;
    if (d >= 0) 
        d = sqrt(d);
        ui(2) = (-b+d) / (2*a);
        ui(1) = (-b-d) / (2*a);
    end    
    % remember alpha at faces
    alpha_face1 = a1;
    alpha_face2 = a2;
    
    % faces 3,4, v=0, v=1
    face = 3;
    [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
    
    e1 = f0+f3-f1-f2;
    v1 = (f0-f1)/e1;
    u1 = (f0-f2)/e1;
    a1 = (f0*f3 - f1*f2)/e1;
    e2 = h0+h3-h1-h2;
    v2 = (h0-h1)/e2;
    u2 = (h0-h2)/e2;
    a2 = (h0*h3 - h1*h2)/e2;
    a = (u2-u1)*e1*e2;
    b = -a*(v1+v2) + (i0-a2)*e1 - (i0-a1)*e2;
    c = a*v1*v2 - (i0-a2)*e1*v1 + (i0-a1)*e2*v2;
    d = b*b - 4*a*c;
    if (d >= 0) 
       d = sqrt(d);
       wi(1) = (-b+d) / (2*a);
       wi(2) = (-b-d) / (2*a);
    end
    % remember alpha at faces
    alpha_face3 = a1;
    alpha_face4 = a2;
    
    % faces 5,6, v=0, v=1
    face = 5;
    [f0, f1, f2, f3, h0, h1, h2, h3] = setBilinearform(F,face);
    
    e1 = f0+f3-f1-f2;
    v1 = (f0-f1)/e1;
    u1 = (f0-f2)/e1;
    a1 = (f0*f3 - f1*f2)/e1;
    e2 = h0+h3-h1-h2;
    v2 = (h0-h1)/e2;
    u2 = (h0-h2)/e2;
    a2 = (h0*h3 - h1*h2)/e2;
    a = (v2-v1)*e1*e2;
    b = -a*(u1+u2) + (i0-a2)*e1 - (i0-a1)*e2;
    c = a*u1*u2 - (i0-a2)*e1*u1 + (i0-a1)*e2*u2;
    d = b*b - 4*a*c;
    if (d >= 0) 
       d = sqrt(d);
       vi(1) = (-b+d) / (2*a);
       vi(2) = (-b-d) / (2*a);
    end
    
    % numerical error
    small = 0.00000000000001;
    for i = 1:2
        if ui(i) < 0 && abs(ui(i)) < small
            ui(i) = 0;
        end
        if ui(i) > 1 && abs(ui(i) - 1) < small
            ui(i) = 1;
        end
        if vi(i) < 0 && abs(vi(i)) < small
            vi(i) = 0;
        end
        if vi(i) > 1 && abs(vi(i) - 1) < small
            vi(i) = 1;
        end
        if wi(i) < 0 && abs(wi(i)) < small
            wi(i) = 0;
        end
        if wi(i) > 1 && abs(wi(i) - 1) < small
            wi(i) = 1;
        end
    end
    
    % remember alpha at faces
%     alpha_face5 = a1;
%     alpha_face6 = a2;
%         
%     % check if singular case
%     if (alpha_face1 == i0) 
%         if wi(1) < wi(2)
%             wi(1) = 0;
%         else
%             wi(2) = 0;
%         end
%     end
%     if (alpha_face2 == i0) 
%         if wi(1) > ui(2)
%             wi(1) = 1;
%         else
%             wi(2) = 1;
%         end
%     end
%     if (alpha_face3 == i0) 
%         if vi(1) < vi(2)
%             vi(1) = 0;
%         else
%             vi(2) = 0;
%         end
%     end
%     if (alpha_face4 == i0) 
%         if vi(1) > vi(2)
%             vi(1) = 1;
%         else
%             vi(2) = 1;
%         end
%     end
%     if (alpha_face5 == i0) 
%         if ui(1) < ui(2)
%             ui(1) = 0;
%         else
%             ui(2) = 0;
%         end
%     end
%     if (alpha_face6 == i0) 
%         if ui(1) > ui(2)
%             ui(1) = 1;
%         else
%             ui(2) = 1;
%         end
%     end
    
    % create 6 points
    cvt(1,:) = [ui(1),vi(1),wi(1)];
    cvt(2,:) = [ui(1),vi(1),wi(2)];
    cvt(3,:) = [ui(2),vi(1),wi(2)];
    cvt(4,:) = [ui(2),vi(2),wi(2)];
    cvt(5,:) = [ui(2),vi(2),wi(1)];
    cvt(6,:) = [ui(1),vi(2),wi(1)];

    % count nr. of solutions
    if (ui(1) >= 0 && ui(1) <= 1) && (vi(1) >= 0 && vi(1) <= 1)
        ct(1) = ct(1) + 1;
    end
    if (ui(2) >= 0 && ui(2) <= 1) && (vi(2) >= 0 && vi(2) <= 1)
        ct(1) = ct(1) + 1;
    end
    if (ui(1) >= 0 && ui(1) <= 1) && (wi(1) >= 0 && wi(1) <= 1)
        ct(2) = ct(2) + 1;
    end
    if (ui(2) >= 0 && ui(2) <= 1) && (wi(2) >= 0 && wi(2) <= 1)
        ct(2) = ct(2) + 1;
    end
    if (vi(1) >= 0 && vi(1) <= 1) && (wi(2) >= 0 && wi(2) <= 1)
        ct(3) = ct(3) + 1;
    end
    if (vi(2) >= 0 && vi(2) <= 1) && (wi(1) >= 0 && wi(1) <= 1)
        ct(3) = ct(3) + 1;
    end
    if ui(1) > 0 && ui(1) < 1
        fs(1,1) = 1;
    end
    if ui(2) > 0 && ui(2) < 1
        fs(1,2) = 1;
    end
    if vi(1) > 0 && vi(1) < 1
        fs(2,1) = 1;
    end
    if vi(2) > 0 && vi(2) < 1
        fs(2,2) = 1;
    end
    if wi(1) > 0 && wi(1) < 1
        fs(3,1) = 1;
    end
    if wi(2) > 0 && wi(2) < 1
        fs(3,2) = 1;
    end
   
    % debugging purpose
    %disp(['ui1 = ',num2str(ui(1)), ', ui2 = ',num2str(ui(2))]);
    %disp(['vi1 = ',num2str(vi(1)), ', vi2 = ',num2str(vi(2))]);
    %disp(['wi1 = ',num2str(wi(1)), ', wi2 = ',num2str(wi(2))]);
end


% compute the distance between two numbers in a periodic ring
function d = distanceInnerHexagon(d1,d2)
    r = d1 - d2;
    if r < 0
        r = -r;
    end
    
    if r > 2
        d = 6 - r;
    else
        d = r;
    end    
end

function d = distanceInnerPolygon(psz, d1, d2)
    r = d1 - d2;
    if r < 0
        r = -r;
    end
    
    if r > 2
        d = psz - r;
    else
        d = r;
    end    
end

% compute the midpoint between to integers in a ring of integers mod 6
function mp = midpointIntegRing(d1,d2)
    dmax = max(d1,d2);
    dmin = min(d1,d2);
    if mod(dmax+2,6) == dmin
        mp = mod(dmax+1,6);
    else
        mp = (dmax + dmin) / 2;
    end
end
% Triangulate contours. For Tunnel and controu with 12 vertices
% use all 6 vertices of inner hexagon to produce better triangulations
% ==%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nrT,tris,verts] = triangulateVts(i0,i1,d1,d2,vt,contour,cvt,flag)
    %dst = distanceInnerHexagon(d1-1,d2-1);
    rt_size = size(cvt,1);
    dst = distanceInnerPolygon(size(cvt,1),d1-1,d2-1);
    if dst == 0
        if flag == true
            nrT = 0;
            tris = zeros(1,3);
            verts = zeros(3,3);
        else
            % there is only one triangle
            nrT = 1;
            Trs = zeros(1,3);
            Trs(:) = [i0,i1,d1];
            verts = zeros(3,3);
            verts(1,:) = vt(contour(i0),1:3);
            verts(2,:) = vt(contour(i1),1:3);
            verts(3,:) = cvt(d1,:);
            tris = [1,2,3];
        end
    elseif dst == 2
        % there are three triangles
        if flag == true
            % eliminate tirangle at border
            nrT = 2;
            Trs = zeros(2,3);
            dm = midpointIntegRing(d1-1,d2-1) + 1;
            Trs(1,:) = [i0,dm,d1];
            Trs(2,:) = [i1,d2,dm];
            verts = zeros(5,3);
            verts(1,:) = vt(contour(i0),1:3);
            verts(2,:) = vt(contour(i1),1:3);
            verts(3,:) = cvt(d2,:);
            verts(4,:) = cvt(d1,:);
            verts(5,:) = cvt(dm,:);
            tris = [[1,5,4];[2,3,5]];
        else 
            nrT = 3;
            Trs = zeros(3,3);
            dm = midpointIntegRing(d1-1,d2-1) + 1;
            Trs(1,:) = [i0,i1,dm];
            Trs(2,:) = [i0,dm,d1];
            Trs(3,:) = [i1,d2,dm];
            verts = zeros(5,3);
            verts(1,:) = vt(contour(i0),1:3);
            verts(2,:) = vt(contour(i1),1:3);
            verts(3,:) = cvt(d2,:);
            verts(4,:) = cvt(d1,:);
            verts(5,:) = cvt(dm,:);
            tris = [[1,2,5];[1,5,4];[2,3,5]];
        end
    elseif dst == 1 % points are consecutive
        % there are only two triangles
        if flag == true
            if rt_size == 6
                nrT = 2;
                Trs = zeros(2,3);
                tris = zeros(2,3);
                Trs(1,:) = [i0,d2,d1];
                tris(1,:) = [1,4,3];
                Trs(2,:) = [i1,d2,d1];
                tris(2,:) = [2,4,3];
                verts = zeros(4,3);
                verts(1,:) = vt(contour(i0),1:3);
                verts(2,:) = vt(contour(i1),1:3);
                verts(3,:) = cvt(d1,:);
                verts(4,:) = cvt(d2,:);
            else 
                nrT = 1;
                Trs = zeros(1,3);
                tris = zeros(1,3);
                % triangulate along the shortest diagonal
                l1 = norm(vt(contour(i0),1:3)-cvt(d2,:));
                l2 = norm(vt(contour(i1),1:3)-cvt(d1,:));
                if l1 < l2
                    Trs(1,:) = [i0,d2,d1];
                    tris = [1,4,3];
                else
                    Trs(1,:) = [i1,d2,d1];
                    tris = [2,4,3];
                end
                verts = zeros(4,3);
                verts(1,:) = vt(contour(i0),1:3);
                verts(2,:) = vt(contour(i1),1:3);
                verts(3,:) = cvt(d1,:);
                verts(4,:) = cvt(d2,:);
            end
        else 
            nrT = 2;
            Trs = zeros(2,3);
            tris = zeros(2,3);
            % triangulate along the shortest diagonal
            l1 = norm(vt(contour(i0),1:3)-cvt(d2,:));
            l2 = norm(vt(contour(i1),1:3)-cvt(d1,:));
            if l1 < l2
                Trs(1,:) = [i0,i1,d2];
                Trs(2,:) = [i0,d2,d1];
                tris = [[1,2,4];[1,4,3]];
            else
                Trs(1,:) = [i0,i1,d1];
                Trs(2,:) = [i1,d2,d1];
                tris = [[1,2,3];[2,4,3]];
            end
            verts = zeros(4,3);
            verts(1,:) = vt(contour(i0),1:3);
            verts(2,:) = vt(contour(i1),1:3);
            verts(3,:) = cvt(d1,:);
            verts(4,:) = cvt(d2,:);
        end
%         t_hdl = patch('Faces',tris,'Vertices',verts,'FaceColor',[0.8,0.8,0.8]);
%         t_handles = [t_handles,t_hdl];
    end

end


