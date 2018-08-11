function pos = piecewise_driver(inp,tim,ini_cond)
format long


global gv

npie        =4;    

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------


% load('C:\Users\Federico.Renda\Documents\MATLAB\Piecewise\Dynamics\Multi Section\2016_11_07 - EXP\TestFede dinamico\ricostruite\scena 2\t_vs_T13')
% load('C:\Users\Federico.Renda\Documents\MATLAB\Piecewise\Dynamics\Multi Section\2016_11_07 - EXP\TestFede dinamico\ricostruite\scena 13\t_vs_T13_T21_T32')
%load('C:\Users\Federico.Renda\Documents\MATLAB\Piecewise\Dynamics\Multi Section\2016_11_07 - EXP\TestFede dinamico\ricostruite\scena 16\t_vs_T11_T13_T31_T33')
t_tens=0:0.01:tim;
T11      = zeros(length(t_tens),1);
T13      = zeros(length(t_tens),1);
T21      =zeros(length(t_tens),1);
T31      =inp(1,:)';
T32      =inp(2,:)';
T33      = inp(3,:)';

% posizione marker
% dist_base   =3.28e-2;               % [m] distanza del marker alla base;
dist_base   =1.51e-2;               % [m] distanza del marker alla base x scena 16;
dist_m1     =1.54e-2;               % [m] distanza del marker 1;
dist_m2     =1.15e-2;               % [m] distanza del marker 2;
dist_m3     =0.4e-2;                % [m] distanza del marker 3;
dist_m4     =0.32e-2;               % [m] distanza del marker 4;

%-------------------------------------------------------------------------
% beginning of input section

% Geometrical input del braccio siliconico
E        =10*110e3;                         % [Pa] modulo di young 110e3
eta      =2*3e2;                           % [Pa*s] viscosità 3e3
Poi      =0.5;                           % [-] modulo di Poisson 0.5
G        =E/(2*(1+Poi));                 % [Pa] modulo di taglio
R1       =13.8e-3;                       % [m] Raggio sezione 1
R2       =11.1e-3;                       % [m] Raggio sezione 2
R3       =8.2e-3;                        % [m] Raggio sezione 3
R4       =5.4e-3;                        % [m] Raggio sezione 4
R        =[R1 R2 R3 R4];
d1       =3e-3;                          % [m] braccio cavi punta
d2       =6e-3;                          % [m] braccio cavi centrale
d3       =9e-3;                          % [m] braccio cavi base
L1       =98e-3;                         % [m] Lunghezza sezione 1
L2       =105e-3;                        % [m] Lunghezza sezione 2
L3       =108e-3;                        % [m] Lunghezza sezione 3
L4       =107e-3;                        % [m] Lunghezza sezione 4
L        =[L1 L2 L3 L4];
nsez     =floor(L3*2e2+1);               % una sezione per mezzo centimetro floor(L*2e2+1)
X1       =linspace(0,L1,nsez);           % [m] curvilinear abscissa sezione 1
X2       =linspace(0,L2,nsez);           % [m] curvilinear abscissa sezione 2
X3       =linspace(0,L3,nsez);           % [m] curvilinear abscissa sezione 3
X4       =linspace(0,L4,nsez);           % [m] curvilinear abscissa sezione 4
X        =[X1 X2 X3 X4];
A1       =pi*R1^2;                       % [m^2] sezione 1
A2       =pi*R2^2;                       % [m^2] sezione 2
A3       =pi*R3^2;                       % [m^2] sezione 3
A4       =pi*R4^2;                       % [m^2] sezione 4
J1       =pi*R1^4/4;                     % [m^4] sezione 1
J2       =pi*R2^4/4;                     % [m^4] sezione 2
J3       =pi*R3^4/4;                     % [m^4] sezione 3
J4       =pi*R4^4/4;                     % [m^4] sezione 4
I1       =pi*R1^4/2;                     % [m^4] sezione 1
I2       =pi*R2^4/2;                     % [m^4] sezione 2
I3       =pi*R3^4/2;                     % [m^4] sezione 3
I4       =pi*R4^4/2;                     % [m^4] sezione 4

%-------------------------------------------------------------------------
% initial configuration t=0

xci_star    =[0;0;0;1;0;0];

%-------------------------------------------------------------------------
% dinamic parameters

ro_water    =100;%1022;                                   % [Kg/m^3] densità
ro_arm      =1080;                                   % [Kg/m^3] densità nominale
Clx         =0.01;                                      % [-] coeff. viscosità longitudinale
Cly         =2.5;                                   % [-] coeff. viscosità trasversale
Clz         =2.5;                                   % [-] coeff. viscosità trasversale
Bly         =1.5;                                    % [-] coeff. massa aggiunta trasversale
Blz         =1.5;                                    % [-] coeff. massa aggiunta trasversale
Gra         =[0;0;0;-9.81;0;0];                      % [m/s^2] vettore gravità [0;0;0;-9.81;0;0]
D1          =ro_water*diag([0 0 0 0.5*pi*R1*Clx R1*Cly R1*Clz]);    % drag coef matrix section 1
D2          =ro_water*diag([0 0 0 0.5*pi*R2*Clx R2*Cly R2*Clz]);    % drag coef matrix section 2
D3          =ro_water*diag([0 0 0 0.5*pi*R3*Clx R3*Cly R3*Clz]);    % drag coef matrix section 3
D4          =ro_water*diag([0 0 0 0.5*pi*R4*Clx R4*Cly R4*Clz]);    % drag coef matrix section 4


% D1          =100*diag([0 0 0 0.5*pi*R1*Clx R1*Cly R1*Clz]);    % drag coef matrix section 1
% D2          =100*diag([0 0 0 0.5*pi*R2*Clx R2*Cly R2*Clz]);    % drag coef matrix section 2
% D3          =100*diag([0 0 0 0.5*pi*R3*Clx R3*Cly R3*Clz]);    % drag coef matrix section 3
% D4          =100*diag([0 0 0 0.5*pi*R4*Clx R4*Cly R4*Clz]);    % drag coef matrix section 4

D           =[D1 D2 D3 D4];
Ag1         =ro_water*diag([0 0 0 0 A1*Bly A1*Blz]); % massa da aggiungere section 1
Ag2         =ro_water*diag([0 0 0 0 A2*Bly A2*Blz]); % massa da aggiungere section 2
Ag3         =ro_water*diag([0 0 0 0 A3*Bly A3*Blz]); % massa da aggiungere section 3
Ag4         =ro_water*diag([0 0 0 0 A4*Bly A4*Blz]); % massa da aggiungere section 4


% Ag1         =1022*diag([0 0 0 0 A1*Bly A1*Blz]); % massa da aggiungere section 1
% Ag2         =1022*diag([0 0 0 0 A2*Bly A2*Blz]); % massa da aggiungere section 2
% Ag3         =1022*diag([0 0 0 0 A3*Bly A3*Blz]); % massa da aggiungere section 3
% Ag4         =1022*diag([0 0 0 0 A4*Bly A4*Blz]); % massa da aggiungere section 4


Eps1        =diag([G*I1 E*J1 E*J1 E*A1 G*A1 G*A1]);  % stifness matrix section 1
Eps2        =diag([G*I2 E*J2 E*J2 E*A2 G*A2 G*A2]);  % stifness matrix section 2
Eps3        =diag([G*I3 E*J3 E*J3 E*A3 G*A3 G*A3]);  % stifness matrix section 3
Eps4        =diag([G*I4 E*J4 E*J4 E*A4 G*A4 G*A4]);  % stifness matrix section 4
Eps         =[Eps1 Eps2 Eps3 Eps4];
Ipsi1       =eta*diag([I1 3*J1 3*J1 3*A1 A1 A1]);    % viscosity matrix section 1
Ipsi2       =eta*diag([I2 3*J2 3*J2 3*A2 A2 A2]);    % viscosity matrix section 2
Ipsi3       =eta*diag([I3 3*J3 3*J3 3*A3 A3 A3]);    % viscosity matrix section 3
Ipsi4       =eta*diag([I4 3*J4 3*J4 3*A4 A4 A4]);    % viscosity matrix section 4
Ipsi        =[Ipsi1 Ipsi2 Ipsi3 Ipsi4];
M1          =ro_arm*diag([I1 J1 J1 A1 A1 A1]);       % inertia matrix section 1
M2          =ro_arm*diag([I2 J2 J2 A2 A2 A2]);       % inertia matrix section 2
M3          =ro_arm*diag([I3 J3 J3 A3 A3 A3]);       % inertia matrix section 3
M4          =ro_arm*diag([I4 J4 J4 A4 A4 A4]);       % inertia matrix section 4
M           =[M1 M2 M3 M4];
Ma1         =M1+Ag1;                                 % added inertia matrix section 1
Ma2         =M2+Ag2;                                 % added inertia matrix section 2
Ma3         =M3+Ag3;                                 % added inertia matrix section 3
Ma4         =M4+Ag4;                                 % added inertia matrix section 4
Ma          =[Ma1 Ma2 Ma3 Ma4];

%-------------------------------------------------------------------------
% numerical setting
time        =floor(t_tens(end));     % [s]
nsol        =time*10^2+1;            % una soluzione ogni centisecondo
tspan       =linspace(0,time,nsol);  % [s] time
                  % numero di pezzi
dX1         =L1/(nsez-1);            % delta X section 1
dX2         =L2/(nsez-1);            % delta X section 2
dX3         =L3/(nsez-1);            % delta X section 3
dX4         =L4/(nsez-1);            % delta X section 4
dX          =[dX1 dX2 dX3 dX4];

%-------------------------------------------------------------------------
% actuation load (body coordinate)

tact        =1;                         % [s] tempo torque in dir z o y
trel        =9.0;                       % [s] tempo rilassamento
Fax         =[-T31-T32-T33 -T21 -T13-T11...
    zeros(length(t_tens),1)]; % [N] contraction load
Fay         =[zeros(length(t_tens),1)...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)]; % [N] lateral y load
Faz         =[zeros(length(t_tens),1)...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)]; % [N] lateral z load
Famx        =[zeros(length(t_tens),1)...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)]; % [Nm] torsion torque
Famy        =[-T32*d3...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)...
    zeros(length(t_tens),1)]; % [Nm] bending torque
Famz        =[(T31-T33)*d3 T21*d2 (T11-T13)*d1...
    zeros(length(t_tens),1)]; % [Nm] bending torque
%-------------------------------------------------------------------------
% external tip load (base (X=0) coordinate)

Fpx         =[0 0 0 0];              % [N] contraction load
Fpy         =[0 0 0 0];              % [N] lateral y load
Fpz         =[0 0 0 0];              % [N] lateral z load
Fpmx        =[0 0 0 0];              % [Nm] torsion torque
Fpmy        =[0 0 0 0];              % [Nm] bending torque
Fpmz        =[0 0 0 0];              % [Nm] bending torque

%-------------------------------------------------------------------------
% osservabili

g           =zeros(4*nsol,4*nsez*npie);
eta         =zeros(6*nsol,nsez*npie);
nstep       =1;

% global variable
gv.ro_arm      =ro_arm;
gv.ro_water    =ro_water;
gv.Gra         =Gra;
gv.L           =L;
gv.dist_base   =dist_base;
gv.dist_m1     =dist_m1;
gv.dist_m2     =dist_m2;
gv.dist_m3     =dist_m3;
gv.dist_m4     =dist_m4;
gv.X           =X;
gv.R           =R;
gv.xci_star    =xci_star;
gv.D           =D;
gv.Eps         =Eps;
gv.Ipsi        =Ipsi;
gv.M           =M;
gv.Ma          =Ma;
gv.nsol        =nsol;
gv.nsez        =nsez;
gv.npie        =npie;
gv.dX          =dX;
gv.t_tens      =t_tens;
gv.time        =time;
gv.tspan       =tspan;
gv.tact        =tact;
gv.trel        =trel;
gv.Fax         =Fax;
gv.Fay         =Fay;
gv.Faz         =Faz;
gv.Famx        =Famx;
gv.Famy        =Famy;
gv.Famz        =Famz;
gv.Fpx         =Fpx;
gv.Fpy         =Fpy;
gv.Fpz         =Fpz;
gv.Fpmx        =Fpmx;
gv.Fpmy        =Fpmy;
gv.Fpmz        =Fpmz;

% osservabili
gv.g           =g;
gv.eta         =eta;
gv.nstep       =nstep;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% solution initialization
% sol=[xci*npie xci_dot*npie]

disp('Time-advancing')
myopt          =odeset('RelTol',1e-4,'OutputFcn',@piecewise_observables);

%-------------------------------------------------------------------------
% condizioni iniziali temporali

xci_0          =[0;0;0;1;0;0];
xcidot_0       =[0;0;0;0;0;0];

% ini_cond        =[repmat(xci_0',[1,npie]) repmat(xcidot_0',[1,npie])];

% integrate
[t,z]          =ode45(@piecewise_derivatives,tspan,ini_cond(1,1:12*npie),myopt);%24


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% postproc


nsol=size(z,1);



%piecewise_postproc(t,z)

ang               =linspace(0,2*pi,180);
g=gv.g;

for ii=1:nsol
    
    g_now              =g(4*(ii-1)+1:4*(ii-1)+4,:);
    
    
    for zz=1:npie
        sez     =[zeros(1,180) 0;R(zz)*sin(ang) 0;R(zz)*cos(ang) 0;ones(1,180) 1];
        %sez     =g_now(:,4*(jj-1)+1:4*(jj-1)+4)*sez;
        for jj=1:nsez
            sez_qui  =g_now(:,4*nsez*(zz-1)+4*(jj-1)+1:4*nsez*(zz-1)+4*(jj-1)+4)*sez;
            % plot3(sez_qui(1,:),sez_qui(2,:),sez_qui(3,:),'Color',[1-mod(zz,2),0,mod(zz,2)])
            pos(:,jj,ii,zz)=[sez_qui(1:3,181)];
        end
        
        
        %        plot3(sez(1,:),sez(2,:),sez(3,:),'Color',[1,0,0])
        %
        %        drawnow
    end
end
% piecewise_postproc(t,z)
%piecewise_postprocsolo(t,z)

end

% end