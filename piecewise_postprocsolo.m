function piecewise_postprocsolo(t,z)
global gv

% global variable
L         =gv.L;
dist_base =gv.dist_base;
dist_m1   =gv.dist_m1;
dist_m2   =gv.dist_m2;
dist_m3   =gv.dist_m3;
dist_m4   =gv.dist_m4;
R         =gv.R;
g         =gv.g;
eta       =gv.eta;
nsez      =gv.nsez;
npie      =gv.npie;
nsol      =gv.nsol;

%-------------------------------------------------------------------------
% pre-processing

% initialization
error_mean      =zeros(nsol,1);
error_tip       =zeros(nsol,1);

% loading test 
% load('C:\Users\Federico.Renda\Documents\MATLAB\Piecewise\Dynamics\Multi Section\2016_11_07 - EXP\TestFede dinamico\ricostruite\scena 2\Scena2_ric2d');
% load('C:\Users\Federico.Renda\Documents\MATLAB\Piecewise\Dynamics\Multi Section\2016_11_07 - EXP\TestFede dinamico\ricostruite\scena 13\Scena13_ric3d');
load('C:\Users\Federico.Renda\Documents\MATLAB\Piecewise\Dynamics\Multi Section\2016_11_07 - EXP\TestFede dinamico\ricostruite\scena 16\Scena16_ric2d');

scena           =points2d;                      % 100 fps
% scena           =Scena13_ric3d;               

g_scena         =[0 1 0 0; 1 0 0 0; 0 0 -1 0];  % from scena to e1-e2-e3

mkdir('.\LAST RUN\');

% get the solution
tau       =zeros(nsol,npie);
xci       =zeros(nsol,npie);
k         =zeros(nsol,npie);
q         =zeros(nsol,npie);
p         =zeros(nsol,npie);
r         =zeros(nsol,npie);
for ii=1:npie
    tau(:,ii) =z(:,6*(ii-1)+1);
    xci(:,ii) =z(:,6*(ii-1)+2);
    k(:,ii)   =z(:,6*(ii-1)+3);
    q(:,ii)   =z(:,6*(ii-1)+4);
    p(:,ii)   =z(:,6*(ii-1)+5);
    r(:,ii)   =z(:,6*(ii-1)+6);
end
xcidot    =z(:,6*npie+1:12*npie);

%-------------------------------------------------------------------------
% save risults

save('.\LAST RUN\postproc','t','z')
save('.\LAST RUN\long strain','t','q');
save('.\LAST RUN\tras n strain','t','p');
save('.\LAST RUN\tras b strain','t','r');
save('.\LAST RUN\torsione','t','tau');
save('.\LAST RUN\curavtura su n','t','xci');
save('.\LAST RUN\curvatura su b','t','k');
save('.\LAST RUN\strain velocity','t','xcidot');

save('.\LAST RUN\Rototraslation','t','g');
save('.\LAST RUN\velocity','t','eta');

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% plots

for ii=1:npie
    % deformations

    figure
    plot(t,tau(:,ii))
    grid on
    title(strcat('torsion of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('tau [1/m]')
    auxstr  =strcat('.\LAST RUN\torsione',num2str(ii),'.png');
    print('-dpng',auxstr)

    figure
    plot(t,xci(:,ii))
    grid on
    title(strcat('curvature on y of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('xci [1/m]')
    auxstr  =strcat('.\LAST RUN\curvature',num2str(ii),'_on_y.png');
    print('-dpng',auxstr)

    figure
    plot(t,k(:,ii))
    grid on
    title(strcat('curvature on z of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('k [1/m]')
    auxstr  =strcat('.\LAST RUN\curvature',num2str(ii),'_on_z.png');
    print('-dpng',auxstr)

    figure
    plot(t,q(:,ii))
    grid on
    title(strcat('longitudinal strain of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('q [-]')
    auxstr  =strcat('.\LAST RUN\longitudinal_strain',num2str(ii),'.png');
    print('-dpng',auxstr)

    figure
    plot(t,p(:,ii))
    grid on
    title(strcat('tras y strain of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('p [-]')
    auxstr  =strcat('.\LAST RUN\tras_y_strain',num2str(ii),'.png');
    print('-dpng',auxstr)

    figure
    plot(t,r(:,ii))
    grid on
    title(strcat('tras z strain of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('r [-]')
    auxstr  =strcat('.\LAST RUN\tras_z_strain',num2str(ii),'.png');
    print('-dpng',auxstr)

end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Video Poseidrone

mov               =VideoWriter(strcat('.\LAST RUN\Section Dynamics'),'MPEG-4');
mov.FrameRate     =10^2;                                    % fps
open(mov)
scrsz             =get(0,'ScreenSize');
figure('color','w','Position',[scrsz(3)/24 2*scrsz(4)/48 11*scrsz(3)/12 9*scrsz(4)/10])

% figure invariants
ang               =linspace(0,2*pi,180);

for ii=1:nsol                               % per ogni istante
    
        
    % get rototraslation
    g_now              =g(4*(ii-1)+1:4*(ii-1)+4,:);
    
    % calcolo posizione marker simulati
    delta_base          =[0;-dist_base;0;1];
    delta_m1            =[0;-dist_m1;0;1];
    delta_m2            =[0;-dist_m2;0;1];
    delta_m3            =[0;-dist_m3;0;1];
    delta_m4            =[0;-dist_m4;0;1];
    
    markersim_base      =g_now(:,1:4)*delta_base;
    markersim_base      =markersim_base(1:3);
    markersim_1         =g_now(:,4*(nsez-1)+1:4*(nsez-1)+4)*delta_m1;
    markersim_1         =markersim_1(1:3);
    markersim_2         =g_now(:,4*nsez+4*(nsez-1)+1:4*nsez+4*(nsez-1)+4)*delta_m2;
    markersim_2         =markersim_2(1:3);
    markersim_3         =g_now(:,4*nsez*2+4*(nsez-1)+1:4*nsez*2+4*(nsez-1)+4)*delta_m3;
    markersim_3         =markersim_3(1:3);
    markersim_4         =g_now(:,4*nsez*3+4*(nsez-1)+1:4*nsez*3+4*(nsez-1)+4)*delta_m4;
    markersim_4         =markersim_4(1:3);
    
    % get real marker
    marker_base         =g_scena*[scena(ii,[1 2]) 0 1]';              % caso 2d
    marker_base         =marker_base(1:3);
    marker_1            =g_scena*[scena(ii,[3 4]) 0 1]';
    marker_1            =marker_1(1:3);
    marker_2            =g_scena*[scena(ii,[5 6]) 0 1]';
    marker_2            =marker_2(1:3);
    marker_3            =g_scena*[scena(ii,[7 8]) 0 1]';
    marker_3            =marker_3(1:3);
    marker_4            =g_scena*[scena(ii,[9 10]) 0 1]';
    marker_4            =marker_4(1:3);
%     marker_base         =g_scena*[scena(ii,[1 2 3]) 1]';                 % caso 3d
%     marker_base         =marker_base(1:3);
%     marker_1            =g_scena*[scena(ii,[4 5 6]) 1]';
%     marker_1            =marker_1(1:3);
%     marker_2            =g_scena*[scena(ii,[7 8 9]) 1]';
%     marker_2            =marker_2(1:3);
%     marker_3            =g_scena*[scena(ii,[10 11 12]) 1]';
%     marker_3            =marker_3(1:3);
%     marker_4            =g_scena*[scena(ii,[13 14 15]) 1]';
%     marker_4            =marker_4(1:3);

    % calcolo errore
    er_base             =sqrt((marker_base(1)-markersim_base(1))^2 ...
                         +(marker_base(2)-markersim_base(2))^2 ...
                         +(marker_base(3)-markersim_base(3))^2);
    er_1                =sqrt((marker_1(1)-markersim_1(1))^2 ...
                         +(marker_1(2)-markersim_1(2))^2 ...
                         +(marker_1(3)-markersim_1(3))^2);
    er_2                =sqrt((marker_2(1)-markersim_2(1))^2 ...
                         +(marker_2(2)-markersim_2(2))^2 ...
                         +(marker_2(3)-markersim_2(3))^2);
    er_3                =sqrt((marker_3(1)-markersim_3(1))^2 ...
                         +(marker_3(2)-markersim_3(2))^2 ...
                         +(marker_3(3)-markersim_3(3))^2);
    er_4                =sqrt((marker_4(1)-markersim_4(1))^2 ...
                         +(marker_4(2)-markersim_4(2))^2 ...
                         +(marker_4(3)-markersim_4(3))^2);                 
    error_mean(ii)      =100*(er_base + er_1 + er_2 + er_3 + er_4)/(4*sum(L));
    error_tip(ii)       =100*er_4/sum(L);
    
    clf
    
    % set the graph options
    set(gca,'CameraPosition',[-0.5*npie*L(1) npie*L(1) -npie*L(1)],...
        'CameraTarget',[0 0 0],...
        'CameraUpVector',[1 0 0])
    axis equal
    grid on
    hold on
    xlabel('E1 [m]')
    ylabel('E2 [m]')
    zlabel('E3 [m]')
    title(strcat('t= ',num2str(t(ii))))
    axis ([-npie*L(1) npie*L(1) -L(1) 1.5*npie*L(1) -npie*L(1) npie*L(1)])

    % disegno la sezione
    for zz=1:npie
        sez     =[zeros(1,180) 0;R(zz)*sin(ang) 0;R(zz)*cos(ang) 0;ones(1,180) 1];
        for jj=1:nsez
            sez_qui  =g_now(:,4*nsez*(zz-1)+4*(jj-1)+1:4*nsez*(zz-1)+4*(jj-1)+4)*sez;
            plot3(sez_qui(1,:),sez_qui(2,:),sez_qui(3,:),'Color',[1-mod(zz,2),0,mod(zz,2)])
        end
        drawnow
    end
    
    % disegno marker simulati
    line([markersim_base(1) markersim_1(1)],[markersim_base(2) markersim_1(2)],[markersim_base(3) markersim_1(3)],'Color','k')
    line([markersim_1(1) markersim_2(1)],[markersim_1(2) markersim_2(2)],[markersim_1(3) markersim_2(3)],'Color','k')
    line([markersim_2(1) markersim_3(1)],[markersim_2(2) markersim_3(2)],[markersim_2(3) markersim_3(3)],'Color','k')
    line([markersim_3(1) markersim_4(1)],[markersim_3(2) markersim_4(2)],[markersim_3(3) markersim_4(3)],'Color','k')
    plot3(markersim_base(1),markersim_base(2),markersim_base(3),'Color','k','Marker','d')
    plot3(markersim_1(1),markersim_1(2),markersim_1(3),'Color','k','Marker','d')
    plot3(markersim_2(1),markersim_2(2),markersim_2(3),'Color','k','Marker','d')
    plot3(markersim_3(1),markersim_3(2),markersim_3(3),'Color','k','Marker','d')
    plot3(markersim_4(1),markersim_4(2),markersim_4(3),'Color','k','Marker','d')
    
    % disegno marker reali
    line([marker_base(1) marker_1(1)],[marker_base(2) marker_1(2)],[marker_base(3) marker_1(3)],'Color','b')
    line([marker_1(1) marker_2(1)],[marker_1(2) marker_2(2)],[marker_1(3) marker_2(3)],'Color','b')
    line([marker_2(1) marker_3(1)],[marker_2(2) marker_3(2)],[marker_2(3) marker_3(3)],'Color','b')
    line([marker_3(1) marker_4(1)],[marker_3(2) marker_4(2)],[marker_3(3) marker_4(3)],'Color','b')
    plot3(marker_base(1),marker_base(2),marker_base(3),'Color','b','Marker','d')
    plot3(marker_1(1),marker_1(2),marker_1(3),'Color','b','Marker','d')
    plot3(marker_2(1),marker_2(2),marker_2(3),'Color','b','Marker','d')
    plot3(marker_3(1),marker_3(2),marker_3(3),'Color','b','Marker','d')
    plot3(marker_4(1),marker_4(2),marker_4(3),'Color','b','Marker','d')
    
    % force drawing
    drawnow
    
    % for movie
    F   = getframe(gcf);
    writeVideo(mov,F);
end

close(mov);

% save and plot error
save('.\LAST RUN\error_mean','error_mean');
save('.\LAST RUN\error_tip','error_tip');

figure
plot(t(1:end-1),error_mean(1:end-1),'b',t(1:end-1),error_tip(1:end-1),'r')
h = legend('Mean Error','Tip Error',1);
set(h,'Interpreter','none','Location','NorthWest')
% title('Bending')
% title('Fetching')
title('Reaching')
ylabel('e [%]')
xlabel('t [s]')
% auxstr  =strcat('.\LAST RUN\','error_bending.png');
% auxstr  =strcat('.\LAST RUN\','error_fetching.png');
auxstr  =strcat('.\LAST RUN\','error_reaching.png');
print('-dpng',auxstr)

% eof