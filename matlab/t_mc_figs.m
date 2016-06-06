
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function t_mc( )
    %
    cstate_ = containers.Map;
    cstate_('LevelsetHandle') = 0;
    cstate_('LineHandles') = [];
    cstate_('TrisHandles') = [];
    %f = figure('Visible','off');
    close all;
    figure('Visible','off');
    hold on;
    % plot coordinate axis
    vm = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    %patch('Vertices',vm,'Faces',fm, 'FaceVertexCData',hsv(6),'FaceColor','none');
    %ha = axes('Position',[0.07,0.07,0.75,0.75]);

    h_bbox = patch('Vertices',vm,'Faces',fm,'FaceColor','none');
    daspect([1 1 1]);
    % set camera
    %view(3);
    view(25,50);
    axis tight;
    axis square;
    camlight;
    lighting gouraud;
    xlabel('u');
    ylabel('v');
    zlabel('w');
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);
    set(gca,'Position',[0.06 0.15 0.7 0.7]);



    % read current state values
    fig_names = {{'fig1.fig','fig1_mesh.fig'},...
        {'fig2.fig','fig2_mesh.fig'},...
        {'fig3.fig','fig3_mesh.fig'},...
        {'fig4.fig','fig4_mesh.fig'},...
        {'fig5.fig','fig5_mesh.fig'},...
        {'fig6.fig','fig6_mesh.fig'},...
        {'fig7.fig','fig7_mesh.fig'},...
        {'fig8.fig','fig8_mesh.fig'},...
        {'fig9.fig','fig9_mesh.fig'},...
        {'fig10.fig','fig10_mesh.fig'},...
        {'fig11.fig','fig11_mesh.fig'},...
        {'fig12.fig','fig12_mesh.fig'}};

    vals = [0.0387,9.8588,9.9994608191478135,0.0097,9.99946,-1.8061,-1.7660,-0.6,1.10918,0.0708,-1.3064,0.0293];
    mc_lt   = [6,25, 25,30,60,22,24,22,60,30,105,105];
    for i = 1:12
        fig_ = i+1;
        % default case
        f = zeros(1,8);
        switch fig_
            case 2 % contour with 6 vertices, MC 3, bit patter 6
                f = [-0.960492426392903 0.793207329559554 0.916735525189067 -0.422761282626275 -0.934993247757551 -0.850129305868777 -0.0367116785741896 -0.656740699156587];
            case 3 % contour with 7 vertices, type A, MC 6, bit patter 25
                f =  [10.2967816247556,9.45145192686147,9.54753711271687,10.6482067822841,9.81494966341055,9.31168538578250,9.80950580411527,10.7451536262220];
            case 4 % contour with 7 vertices, type B, MC 6, bit pattern 25
                f = [9.9998593195995547,9.9993381282115549,9.9979160205452544,9.9986053863704142,9.9999374908631235,9.999424800002032,9.9983922749132219,9.999579324965488];
            case 5 % contour with 8 vertices, MC 12, bit patter 30
                f = [0.454797708726920,	0.801330575352402,	0.649991492712356,	-0.973974554763863,	-0.134171007607162,	-0.0844698148589140,	-0.826313795402046,	0.433391503783462];
            case 6 % contour with 8 vertices, case b, MC 10, bit patter 60
                f = [9.9985934885536665,9.9998695572230147,9.9999045831713928,9.999316745478131,9.9986117521866866,9.9998754368055813,9.9999031760062458,9.9992041920402936];
            case 7 % contour with 9 vertices, MC 7, bit patter 22
                f = [-15.6504952739285,2.90290077342601,24.5454566157887,-24.5274127623786,21.6741877710053,-4.49696327433901,-19.7891575872492,-15.5588482753161];
            case 8 % MC 4 with tunnel, bit pattern 24
                f = [-7.70146936482581,-3.21868369245987,-5.44023748418735,15.6051950593180,12.7611835388515,-4.46952393442309,-11.7240576326183,-9.23038948829007];
            case 9 % MC 7 with tunnel, bit patter 22
                f = [-3.42744283804455,0.621278122151001,4.48110777981235,-1.95551129669134,2.30448107596369,-1.04182240925489,-3.51087814405650,-6.44976786808517];
            case 10 % MC 10 with tunnel, bit patter 60
                f = [-0.100000000000000,-6.11000000000000,2,10.2000000000000,10.8000000000000,1.80000000000000,-8.20000000000000,-0.180000000000000];
            case 11 % MC 12 with tunnel, bit patter 30
                f=[-3.37811990337124,0.473258332744286,2.54344310345736,7.87658724379480,4.38700713005133,-1.49950251870885,-4.21025867362045,-1.00233824192217];
            case 12 % MC 13 with tunnel, MC 3, bit patter 6
                f = [2.74742516087490,-3.39187542578189,-12.5297639669456,0.431517989649243,-6.92460546400188,2.52228314017858,14.6950568276448,-10.0732624062474];
            case 13 % MC 13 with contour with 12 vertices
                f = [0.546912886195662,	-0.421103532406922,	-0.643375084081520,	0.855507421818445,	-0.260686312588506,	0.206413666735986,	0.237274227130530,	-0.183297728364877];
                %f = [4.69314856679690,-25.6666841753773,-19.3439361061026,11.2891663083649,-5.72871085708909,12.8485897893816,14.4616618309557,-3.61934839891487];
                %f = [2.56647391270132,-7.87546704094998,-24.0314386830922,0.877608326864389,-27.8666241843413,21.9109258856636,14.6592692141074,-17.3567518307032];

        end


        % set a nice color
        c_pos = randi(25);
        color = getColor(c_pos);

        % the isovalue and MC case
        i0 = vals(i);
        mc = mc_lt(i);

        % compute level set within unit cell
        levelset_ = 1;
        % generate isosurface
        if levelset_ == 1
            cc =0:0.02:1;
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

            p_handle = patch(isosurface(x,y,z,fd,i0));
            cstate_('LevelsetHandle') = p_handle;
            set(p_handle,'FaceColor',color,'EdgeColor','none');
            %set(p_handle,'FaceColor',colors(c_pos,:,:,:),'EdgeColor','none');
            % plot straight lines joining opposite faces
            plotStraightLines(cstate_,f,i0);
            savefig(fig_names{i}{1});
            clearLevelSetPlot(cstate_);
        end

        % plot triangulation
        celltris_ = 1;
        if celltris_
            clearTriangles(cstate_);
            triangulateCell(cstate_,f,mc,i0)
            savefig(fig_names{i}{2});
            clearTriangles(cstate_);
        end
	end
end


% Trilinear interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = trilinear(f,u,v,w)
    val = (1-w)*((1-v)*(f(1)*(1-u)+f(2)*u)+v*(f(3)*(1-u)+f(4)*u))+ w*((1-v)*(f(5)*(1-u)+f(6)*u)+v*(f(7)*(1-u)+f(8)*u));
end

%
% Compute the bilinear coefficients for a pair of opposite faces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Plot straight lines lying on the level set joining
% opposite faces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotStraightLines(cstate_,F,i0)
    ui = zeros(3,2);
    vi = zeros(3,2);
    ic = zeros(1,3);
    idx = 1;
    for f = 1:2:6
        % init
        [f1, f2, f3, f4, h1, h2, h3, h4] = setBilinearform(F,f);
        a = (f1-f2)*(-h3+h4+h1-h2)-(h1-h2)*(-f3+f4+f1-f2);
        b = (i0-f1)*(-h3+h4+h1-h2)+(f1-f2)*(h3-h1)-(i0-h1)*(-f3+f4+f1-f2)-(h1-h2)*(f3-f1);
        c = (i0-f1)*(h3-h1)-(i0-h1)*(f3-f1);
        d = b^2 - 4*a*c;
        if d > 0
            pos = 1;
            d = sqrt(d);
            % compute ui
            ui1 = (-b-d) / (2*a);
            ui2 = (-b+d) / (2*a);
            us1 = min(ui1,ui2);
            us2 = max(ui1,ui2);
            if 0 < us1 && us1 < 1
                g1 = f1*(1-us1) + f2*us1;
                g2 = f3*(1-us1) + f4*us1;
                vs1 = (i0 - g1)/(g2-g1);
                if 0 < vs1 && vs1 < 1
                    ui(idx,pos) = us1;
                    vi(idx,pos) = vs1;
                    ic(idx) = ic(idx) + 1;
                    pos = pos + 1;
                end
            end
            if 0 < us2 && us2 < 1
                g1 = f1*(1-us2) + f2*us2;
                g2 = f3*(1-us2) + f4*us2;
                vs2 = (i0 - g1)/(g2-g1);
                if 0 < vs2 && vs2 < 1
                    ui(idx,pos) = us2;
                    vi(idx,pos) = vs2;
                    ic(idx) = ic(idx) + 1;
                end
            end
        end
        idx  = idx + 1;
    end

    % check if there is a tunnel
    ict = ic(1) + ic(2) + ic(3);

    % plot straight lines joining opposite faces
    lt_handles = [];
    for i = 1:ic(1)
        cu1 = [ui(1,i),ui(1,i)];
        cv1 = [vi(1,i),vi(1,i)];
        cw1 = [0, 1];
        l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
        lt_handles = [lt_handles,l_];
    end
    for i = 1:ic(2)
        cu1 = [ui(2,i),ui(2,i)];
        cw1 = [vi(2,i),vi(2,i)];
        cv1 = [0, 1];
        l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
        lt_handles = [lt_handles,l_];
    end
    for i = 1:ic(3)
        cv1 = [ui(3,i),ui(3,i)];
        cw1 = [vi(3,i),vi(3,i)];
        cu1 = [0, 1];
        l_ = plot3(cu1,cv1,cw1,'b','LineWidth',2);
        lt_handles = [lt_handles,l_];
    end
    cstate_('LineHandles') = lt_handles;
end % plotStaightLines(cstate_,F,i0)


% Delete graphic objects from figure
%
function clearLevelSetPlot(cstate_)
    p_handle = cstate_('LevelsetHandle');
	if p_handle ~= 0
        delete(p_handle);
        cstate_('LevelsetHandle') = 0;
    end
    lt_handles = cstate_('LineHandles');
    for l =1:length(lt_handles)
        delete(lt_handles(l));
    end
    lt_handles = [];
    cstate_('LinesHandles') = lt_handles;
end

function clearTriangles(cstate_)
    t_handles = cstate_('TrisHandles');
    for t = 1:length(t_handles)
        delete(t_handles(t));
    end
    cstate_('TrisHandles') = [];
end

%
%
function c_ = getColor(c_index)
    c_table = [[0,117,220];[153,63,0];[116,0,132];[103,100,120];...
               [0,92,49];[43,206,72];[255,204,153];[138,128,158];...
               [148,255,181];[143,124,0];[157,204,0];[194,0,136];...
               [30,81,158];[255,164,5];[255,168,187];[66,102,0];...
               [255,0,16];[94,241,242];[0,153,143];[224,255,102];...
               [116,10,255];[153,0,0];[255,255,128];[255,255,0];...
               [255,80,5]];
    if c_index > 26
        c_index = 1;
    end
    c_ = c_table(c_index,:)/255;
end

%
%
% Compute contours of isosurface at faces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triangulateCell(cstate_,f,mc,isov)
    [nrv,vt,nrc,cns,contours] = computeContorus(f,mc,isov);
    triangulateContours(cstate_,f,mc,isov,nrv,vt,nrc,cns,contours);
end

% Compute contours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% check if for this configuration there is a tunnel
% boolean flag tells if two contours build a tunnel
% nr is the number of countours building a tunnel, it must be two
% ltc contains the id of the contours which build a tunnel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triangulateContours(cstate_,F,mc,i0,nrv,vt,nrc,cns,contours)
    % find out which contours build up a single piece of surface
    nr = 0;
    ltc = 0;
    face = 0;

    % handles to geometry
    t_handles = [];

    % compute all 6 hyperbola intersections
    ui = zeros(3,2);
    vi = zeros(3,2);
    ic = zeros(1,3);
    idx = 1;
    for f = 1:2:6
        % init
        [f1, f2, f3, f4, h1, h2, h3, h4] = setBilinearform(F,f);
        a = (f1-f2)*(-h3+h4+h1-h2)-(h1-h2)*(-f3+f4+f1-f2);
        b = (i0-f1)*(-h3+h4+h1-h2)+(f1-f2)*(h3-h1)-(i0-h1)*(-f3+f4+f1-f2)-(h1-h2)*(f3-f1);
        c = (i0-f1)*(h3-h1)-(i0-h1)*(f3-f1);
        d = b^2 - 4*a*c;
        if d > 0
            pos = 1;
            d = sqrt(d);
            % compute ui
            ui1 = (-b-d) / (2*a);
            ui2 = (-b+d) / (2*a);
            us1 = min(ui1,ui2);
            us2 = max(ui1,ui2);
            if 0 < us1 && us1 < 1
                g1 = f1*(1-us1) + f2*us1;
                g2 = f3*(1-us1) + f4*us1;
                vs1 = (i0 - g1)/(g2-g1);
                if 0 < vs1 && vs1 < 1
                    ui(idx,pos) = us1;
                    vi(idx,pos) = vs1;
                    ic(idx) = ic(idx) + 1;
                    pos = pos + 1;
                end
            end
            if 0 < us2 && us2 < 1
                g1 = f1*(1-us2) + f2*us2;
                g2 = f3*(1-us2) + f4*us2;
                vs2 = (i0 - g1)/(g2-g1);
                if 0 < vs2 && vs2 < 1
                    ui(idx,pos) = us2;
                    vi(idx,pos) = vs2;
                    ic(idx) = ic(idx) + 1;
                end
            end
        end
        idx  = idx + 1;
    end

    % check if there is a tunnel
    ict = ic(1) + ic(2) + ic(3);


    % triangulate tunnel
    if ict == 6
        % if there are three contours, one does not belong to the tunnel
        c1 = 1;
        c2 = 2;
        c3 = 3;
        if nrc == 3
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
                    if ui(1,1) > umax || ui(1,2) < umin
                        c3 = c; % c3 does not belong to the tunnel
                    else
                        c2 = c;
                    end
                else
                    c1 = c;
                end
            end
        end
        % create 6 vertices
        cvt = zeros(6,3);
        % face 1,2 and face 3,4 common coordinate is u
        u = ui(1,1);
        v = vi(1,1);
        w = vi(2,1);
        p1 = 1;
        p2 = 1;
        p3 = 1;
        cvt(1,:) = [u,v,w];
        % face 3,4 and face 5,6 common coord is w
        if abs(w - vi(3,2)) < 0.00005
            p3 = 2;
        end

        % connect alternating pi
        % get new v coordinate from face 4,5
        v = ui(3,p3);
        cvt(2,:) = [u,v,w];
        % get new u coordinate from face 0,1
        p1 = mod(p1,2) + 1;
        u = ui(1,p1);
        cvt(3,:) = [u,v,w];
        % get new w coordinate from face 2,3
        p2 = mod(p2,2) + 1;
        w = vi(2,p2);
        cvt(4,:) = [u,v,w];
        % get new v coordinate from face 4,5
        p3 = mod(p3,2) + 1;
        v = ui(3,p3);
        cvt(5,:) = [u,v,w];
        % get nuew u coordinate from face 0,1
        p1 = mod(p1,2) + 1;
        u = ui(1,p1);
        cvt(6,:) = [u,v,w];

        % fussion verties
        fvert = zeros(6,3);
        fvert(1,:) = (cvt(1,:) + cvt(2,:))/2;
        fvert(2,:) = (cvt(1,:) + cvt(2,:))/2;
        fvert(3,:) = (cvt(3,:) + cvt(4,:))/2;
        fvert(4,:) = (cvt(3,:) + cvt(4,:))/2;
        fvert(5,:) = (cvt(5,:) + cvt(6,:))/2;
        fvert(6,:) = (cvt(5,:) + cvt(6,:))/2;

        % generate triangle points
        tp = zeros(3,3);
        tp(1,:) = (cvt(1,:) + cvt(2,:))/2;
        tp(2,:) = (cvt(3,:) + cvt(4,:))/2;
        tp(3,:) = (cvt(5,:) + cvt(6,:))/2;


        flag_special_case = false;
        if (cns(1) > 6 || cns(2) > 6 ||cns(3) > 6)
            flag_special_case = true;
            %disp('case with one contour');
        end

        % collect vertex pairs
        cid = [c1,c2];
        for cq = 1:2
            xc = cid(cq);
            gconne = zeros(1,12);
            for i = 1:cns(xc)
                index = -1;
                d = 3;
                for j = 1:6
                    val = norm(vt(contours(xc,i),1:3)-cvt(j,:));
                    if d > val
                        index = j;
                        d = val;
                    end
                end
                gconne(i) = index;
            end
            % merged triangles
            for i = 1:cns(xc)
                i0 = i;
                i1 = mod(i,cns(xc)) + 1;
                d1 = gconne(i0);
                d2 = gconne(i1);
                if d1 == d2 || int32(d1)/int32(2) == int32(d2)/int32(2)
                    verts = [vt(contours(xc,i0),1:3);vt(contours(xc,i1),1:3);fvert(d1,:)];
                    triTunnel = zeros(1,3);
                    triTunnel(1,:) = [1,2,3];
                    t_hdl = patch('Faces',triTunnel,'Vertices',verts,'FaceColor',[0.8,0.8,0.8]); %,'EdgeColor','b');
                    t_handles = [t_handles,t_hdl];
                else
                    tris = zeros(2,3);
                    tris(1,:) = [1,2,4];
                    tris(2,:) = [2,3,4];
                    verts = zeros(4,3);
                    verts(1,:) = vt(contours(xc,i0),1:3);
                    verts(2,:) = vt(contours(xc,i1),1:3);
                    verts(3,:) = fvert(d2,:);
                    verts(4,:) = fvert(d1,:);
                    t_hdl = patch('Faces',tris,'Vertices',verts,'FaceColor',[0.8,0.8,0.8]);
                    t_handles = [t_handles,t_hdl];
                end
            end
            % triangulate contour with 12 vertices
            if cns(1) == 12
                tri12 = zeros(1,3);
                tri12(1,:) = [1,2,3];
                vt12 = zeros(3,3);
                vt12(1,:) = fvert(1,:);
                vt12(2,:) = fvert(3,:);
                vt12(3,:) = fvert(5,:);
                t_hdl = patch('Faces',tri12,'Vertices',vt12,'FaceColor',[0.8,0.8,0.8]); %,'EdgeColor','b');
                t_handles = [t_handles,t_hdl];
            end
        end
        % if anny triangulate third contour
        if (nrc == 3)
            trOout = zeros(1,3);
            trOout(1,:) = [1,2,3];
            vtOut = zeros(3,3);
            vtOut(1,:) = vt(contours(c3,1),1:3);
            vtOut(2,:) = vt(contours(c3,2),1:3);
            vtOut(3,:) = vt(contours(c3,3),1:3);
            t_hdl = patch('Faces',trOout,'Vertices',vtOut,'FaceColor',[0.8,0.8,0.8]);
            t_handles = [t_handles,t_hdl];
        end
    else % there is no tunnel
        if ((ict == 2 && ic(1) == 2) || (ict == 2 && ic(2) == 2)  ||  (ict == 2 && ic(3) == 2) || ict < 2)
            % there is no saddle point, this is a simple polygon
            for s = 1:nrc
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
            end
        else % there are some where a saddle point
            % plot bary center
            ucoord = 0;
            vcoord = 0;
            wcoord = 0;
            if (ic(1) == 2) % common coordinate is w
                ucoord = ui(2,1);
                vcoord = ui(3,1);
                wcoord = vi(2,1);
            elseif (ic(2) == 2) % common coordinate is v
                ucoord = ui(1,1);
                vcoord = vi(1,1);
                wcoord = vi(3,1);
            elseif (ic(3) == 2) % common coordinate is u
                ucoord = ui(1,1);
                vcoord = vi(1,1);
                wcoord = vi(2,1);
            else
                ucoord = ui(1,1) + ui(1,2) + ui(2,1) + ui(2,2);
                vcoord = vi(1,1) + vi(1,2) + ui(3,1) + ui(3,2);
                wcoord = vi(2,1) + vi(2,2) + vi(3,1) + vi(3,2);
                ucoord = ucoord/(ic(1) + ic(2));
                vcoord = vcoord/(ic(1) + ic(3));
                wcoord = wcoord/(ic(2) + ic(3));
            end
            %
            % At this point the code is not so clear
            % we are triangulating also non ambiguous cases
            % which are handled different in the c-code
            % this is just a demonstration software
            % star triangulation
            for s = 1:nrc
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
                    vtStar(cns(s)+1,:) =  [ucoord,vcoord,wcoord];

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
        % store handles to manage geometry objects

    end

    % keep track of graphic obejcts
    cstate_('TrisHandles') = t_handles;
end % triangulateContours


%
% Compute asymptotes at faces
% The functions accept arrays as input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uw,vw,ew,aw] = wSurface(f,w)
    f0 = f(1);
    f1 = f(2);
    f2 = f(3);
    f3 = f(4);
    f4 = f(5);
    f5 = f(6);
    f6 = f(7);
    f7 = f(8);
    ew = (f0+f3-f1-f2)*(1-w) + (f4+f7-f5-f6)*w;
    uw = ((f0-f2)*(1-w)+(f4-f6)*w) ./ ew;
    vw = ((f0-f1)*(1-w)+(f4-f5)*w) ./ ew;
    aw = ((f0*(1-w)+f4*w).*(f3*(1-w)+f7*w) - (f1*(1-w)+f5*w).*(f2*(1-w)+f6*w)) ./ ew;
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
    eu = (f0+f6-f2-f4)*(1-u) + (f1+f7-f3-f5)*u;
    vu = ((f0-f4)*(1-u)+(f1-f5)*u) ./ eu;
    wu = ((f0-f2)*(1-u)+(f1-f3)*u) ./ eu;
    au = ((f0*(1-u)+f1*u).*(f6*(1-u)+f7*u) - (f2*(1-u)+f3*u).*(f4*(1-u)+f5*u)) ./ eu;
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
    ev = (f0+f5-f1-f4)*(1-v) + (f2+f7-f3-f6)*v;
    uv = ((f0-f4)*(1-v)+(f2-f6)*v) ./ ev;
    wv = ((f0-f1)*(1-v)+(f2-f3)*v) ./ ev;
    av = ((f0*(1-v)+f2*v).*(f5*(1-v)+f7*v) - (f1*(1-v)+f3*v).*(f4*(1-v)+f6*v)) ./ ev;
end

% computes u0 and v0 where the asymptotes intersect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = transform(index,f)
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
    ft(1:8) = f(1:8);
    for i=1:8
        f(i) = ft(rot(index,i));
    end

end
