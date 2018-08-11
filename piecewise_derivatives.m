function z_punto = piecewise_derivatives(t,z)
global gv

t

% global variable
ro_arm      =gv.ro_arm;
ro_water    =gv.ro_water;
L           =gv.L;
D           =gv.D;
Eps         =gv.Eps;
Ipsi        =gv.Ipsi;
M           =gv.M;
Ma          =gv.Ma;
xci_star    =gv.xci_star;
Gra         =gv.Gra;
dX          =gv.dX;
X           =gv.X;
t_tens      =gv.t_tens;
nsez        =gv.nsez;
npie        =gv.npie;
tact        =gv.tact;
trel        =gv.trel;
Fax         =gv.Fax;
Fay         =gv.Fay;
Faz         =gv.Faz;
Famx        =gv.Famx;
Famy        =gv.Famy;
Famz        =gv.Famz;
Fpx         =gv.Fpx;
Fpy         =gv.Fpy;
Fpz         =gv.Fpz;
Fpmx        =gv.Fpmx;
Fpmy        =gv.Fpmy;
Fpmz        =gv.Fpmz;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% actual solution xci xcidot

Xci              =z(1:6*npie);
Xcidot           =z(6*npie+1:12*npie);             

% inizializzo coefficienti dinamica
genMagM          =zeros(6*npie,6*npie);
genCo1M          =zeros(6*npie,6*npie);
genCo2M          =zeros(6*npie,6*npie);
genDraM          =zeros(6*npie,6*npie);
genForV          =zeros(6*npie,1);
genGraV          =zeros(6*npie,6);
genTipV          =zeros(6*npie,1);

% inizializzo la cnematica precedente
% g_r              =diag([-1 -1 1 1]);         % a testa in giu`
g_r              =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % cantilever
Jaco_prec        =diag([1 1 1 1 1 1 zeros(1, 6*(npie-1))]);
g_prec           =diag([1 1 1 1]);
eta_prec         =zeros(6,1);
adetan_prec      =zeros(6*npie,6*npie);

%-------------------------------------------------------------------------
% calcolo le componenti dei coefficienti dinamica

% mass e coriolis 1 del primo pezzo
xci1             =Xci(1:6,:);
xcidot1          =Xcidot(1:6,:);
k1               =xci1(1:3);
theta1           =sqrt(k1'*k1); 

LMasX            =zeros(6,6*nsez);
LRMagX           =zeros(6,6*nsez);
LRCo1X           =zeros(6,6*nsez);
LRDraX           =zeros(6,6*nsez);
LMas_prec        =zeros(6,6);
LRMag_prec       =zeros(6,6);
LRDra_prec       =zeros(6,6);
LRCo1_prec       =zeros(6,6);
for ii=1:nsez
    coAdjg1_here                     =piecewise_coAdjoint(X(ii),theta1,xci1);
    invAdjg1_here                    =piecewise_invAdjoint(X(ii),theta1,xci1);
    intdAdjg1_here                   =piecewise_ADJ(X(ii),theta1,xci1);
        
    % mass
    Mas_here                         =coAdjg1_here*M(:,1:6)*invAdjg1_here;
    LMas_here                        =intdAdjg1_here'*Mas_here;
    trapz                            =dX(1)*(LMas_prec+LMas_here)/2;
    LMasX(:,6*(ii-1)+1:6*nsez)       =LMasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    LMas_prec                        =LMas_here;
    Mag_here                         =coAdjg1_here*Ma(:,1:6)*invAdjg1_here;
    LRMag_here                       =intdAdjg1_here'*Mag_here*intdAdjg1_here;
    trapz                            =dX(1)*(LRMag_prec+LRMag_here)/2;
    LRMagX(:,6*(ii-1)+1:6*nsez)      =LRMagX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    LRMag_prec                       =LRMag_here;
        
    % coriolises 1
    eta_here                         =eta_prec+intdAdjg1_here*xcidot1;
    LRCo1_here                       =intdAdjg1_here'*dinamico_coadj(eta_here)*Mag_here*intdAdjg1_here;
    trapz                            =dX(1)*(LRCo1_prec+LRCo1_here)/2;
    LRCo1X(:,6*(ii-1)+1:6*nsez)      =LRCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    LRCo1_prec                       =LRCo1_here;
    
    % drag
    Dra_here                         =coAdjg1_here*D(:,1:6)*invAdjg1_here*norm(eta_here(4:6));
    LRDra_here                       =intdAdjg1_here'*Dra_here*intdAdjg1_here;
    trapz                            =dX(1)*(LRDra_prec+LRDra_here)/2;
    LRDraX(:,6*(ii-1)+1:6*nsez)      =LRDraX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    LRDra_prec                       =LRDra_here;
    
end
LMasX            =LMasX-repmat(LMasX(:,1:6),[1,nsez]);
LRMagX           =LRMagX-repmat(LRMagX(:,1:6),[1,nsez]);
LRCo1X           =LRCo1X-repmat(LRCo1X(:,1:6),[1,nsez]);
LRDraX           =LRDraX-repmat(LRDraX(:,1:6),[1,nsez]);
LMas             =LMasX(:,6*(nsez-1)+1:6*nsez);
LRMag            =LRMagX(:,6*(nsez-1)+1:6*nsez);
LRCo1            =LRCo1X(:,6*(nsez-1)+1:6*nsez);
LRDra            =LRDraX(:,6*(nsez-1)+1:6*nsez);

% Actuation load internal load and tip load del primo pezzo
%if t<=tact                                      % virata
  %  Fa1          =[Famx(:,1) Famy(:,1) Famz(:,1) Fax(:,1) Fay(:,1) Faz(:,1)]*t/tact;
 %   Fa1          =interp1(t_tens,Fa1,t)';
%end
%if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
    Fa1          =[Famx(:,1) Famy(:,1) Famz(:,1) Fax(:,1) Fay(:,1) Faz(:,1)];
    Fa1          =interp1(t_tens,Fa1,t)';
%end

Fi1              =Eps(:,1:6)*(xci1-xci_star)+Ipsi(:,1:6)*xcidot1;
Fp1              =[Fpmx(1);Fpmy(1);Fpmz(1);Fpx(1);Fpy(1);Fpz(1)];
    
% Actuation load successiva
if npie ~= 1
    %if (t<=tact)                                      % virata
     %   Fa1_suc     =[Famx(:,2) Famy(:,2) Famz(:,2) Fax(:,2) Fay(:,2) Faz(:,2)]*t/tact;
      %  Fa1_suc     =interp1(t_tens,Fa1_suc,t)';
    %end
    %if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
        Fa1_suc     =[Famx(:,2) Famy(:,2) Famz(:,2) Fax(:,2) Fay(:,2) Faz(:,2)];
        Fa1_suc     =interp1(t_tens,Fa1_suc,t)';
    %end
else
    Fa1_suc      =[0;0;0;0;0;0];
end

% aggiorno coefficienti dinamica
invAdjg1_last    =piecewise_invAdjoint(X(nsez),theta1,xci1);
invAdjg1R_last   =blkdiag(invAdjg1_last(1:3,1:3),invAdjg1_last(4:6,4:6));
intdAdjg1_last   =piecewise_ADJ(X(nsez),theta1,xci1);
MagB             =blkdiag(LRMag,zeros(6*(npie-1),6*(npie-1)));
genMagM          =genMagM+Jaco_prec'*MagB*Jaco_prec;
Co1B             =blkdiag(LRCo1,zeros(6*(npie-1),6*(npie-1)));
genCo1M          =genCo1M+Jaco_prec'*Co1B*Jaco_prec;
GraB             =[LMas;zeros(6*(npie-1),6)];
genGraV          =genGraV+Jaco_prec'*GraB*dinamico_Adjoint(g_prec^-1);
DraB             =blkdiag(LRDra,zeros(6*(npie-1),6*(npie-1)));
genDraM          =genDraM+Jaco_prec'*DraB*Jaco_prec;
% ForB             =[(invAdjg1_last*intdAdjg1_last)'*(Fa1-Fa1_suc);zeros(6*(npie-1),1)]; % pneumatico
% genForV          =genForV+Jaco_prec'*ForB-[L*Fi1;zeros(6*(npie-1),1)];                 % pneumatico
genForV          =genForV+[L(1)*(Fa1-Fi1);zeros(6*(npie-1),1)];  % cavo tip2base
TipB             =[(invAdjg1_last*intdAdjg1_last)'*(invAdjg1R_last*Fp1);zeros(6*(npie-1),1)];
genTipV          =genTipV+Jaco_prec'*TipB;

%----------------------------------------------------------------------
% recursive factors
if npie ~= 1
    Jaco_prec       =blkdiag(invAdjg1_last*intdAdjg1_last,zeros(6*(npie-1),6*(npie-1)))*Jaco_prec+...
                     blkdiag(zeros(6,6),diag([1 1 1 1 1 1]),zeros(6*(npie-2),6*(npie-2)));
    g_prec          =g_prec*piecewise_expmap(X(nsez),theta1,xci1);
    eta_prec        =invAdjg1_last*(eta_prec+intdAdjg1_last*xcidot1);
end

%--------------------------------------------------------------------------
% masses, coriolises 1, coriolises 2 dal secondo pezzo in poi

for jj=2:npie                   
    xcin            =Xci(6*(jj-1)+1:6*(jj-1)+6,:);
    xcidotn         =Xcidot(6*(jj-1)+1:6*(jj-1)+6,:);
    kn              =xcin(1:3);
    thetan          =sqrt(kn'*kn);
    
    MasX            =zeros(6,6*nsez);
    LMasX           =zeros(6,6*nsez);
    MagX            =zeros(6,6*nsez);
    LMagX           =zeros(6,6*nsez);
    RMagX           =zeros(6,6*nsez);
    LRMagX          =zeros(6,6*nsez);
    Co1X            =zeros(6,6*nsez);
    LCo1X           =zeros(6,6*nsez);
    RCo1X           =zeros(6,6*nsez);
    LRCo1X          =zeros(6,6*nsez);
    Co2X            =zeros(6,6*nsez);
    LCo2X           =zeros(6,6*nsez);
    DraX            =zeros(6,6*nsez);
    LDraX           =zeros(6,6*nsez);
    RDraX           =zeros(6,6*nsez);
    LRDraX          =zeros(6,6*nsez);
    Mas_prec        =zeros(6,6);
    LMas_prec       =zeros(6,6);
    Mag_prec        =zeros(6,6);
    LMag_prec       =zeros(6,6);
    RMag_prec       =zeros(6,6);
    LRMag_prec      =zeros(6,6);
    Co1_prec        =zeros(6,6);
    LCo1_prec       =zeros(6,6);
    RCo1_prec       =zeros(6,6);
    LRCo1_prec      =zeros(6,6);
    Co2_prec        =zeros(6,6);
    LCo2_prec       =zeros(6,6);
    Dra_prec        =zeros(6,6);
    LDra_prec       =zeros(6,6);
    RDra_prec       =zeros(6,6);
    LRDra_prec      =zeros(6,6);
    for ii=1:nsez
        coAdjgn_here                     =piecewise_coAdjoint(X((jj-1)*nsez+ii),thetan,xcin);
        invAdjgn_here                    =piecewise_invAdjoint(X((jj-1)*nsez+ii),thetan,xcin);
        intdAdjgn_here                   =piecewise_ADJ(X((jj-1)*nsez+ii),thetan,xcin);
        
        % masses
        Mas_here                         =coAdjgn_here*M(:,(jj-1)*6+1:(jj-1)*6+6)*invAdjgn_here;
        trapz                            =dX(jj)*(Mas_prec+Mas_here)/2;
        MasX(:,6*(ii-1)+1:6*nsez)        =MasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Mas_prec                         =Mas_here;
        LMas_here                        =intdAdjgn_here'*Mas_here;
        trapz                            =dX(jj)*(LMas_prec+LMas_here)/2;
        LMasX(:,6*(ii-1)+1:6*nsez)       =LMasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LMas_prec                        =LMas_here;
        Mag_here                         =coAdjgn_here*Ma(:,(jj-1)*6+1:(jj-1)*6+6)*invAdjgn_here;
        trapz                            =dX(jj)*(Mag_prec+Mag_here)/2;
        MagX(:,6*(ii-1)+1:6*nsez)        =MagX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Mag_prec                         =Mag_here;
        LMag_here                        =intdAdjgn_here'*Mag_here;
        trapz                            =dX(jj)*(LMag_prec+LMag_here)/2;
        LMagX(:,6*(ii-1)+1:6*nsez)       =LMagX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LMag_prec                        =LMag_here;
        RMag_here                        =Mag_here*intdAdjgn_here;
        trapz                            =dX(jj)*(RMag_prec+RMag_here)/2;
        RMagX(:,6*(ii-1)+1:6*nsez)       =RMagX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        RMag_prec                        =RMag_here;
        LRMag_here                       =intdAdjgn_here'*Mag_here*intdAdjgn_here;
        trapz                            =dX(jj)*(LRMag_prec+LRMag_here)/2;
        LRMagX(:,6*(ii-1)+1:6*nsez)      =LRMagX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LRMag_prec                       =LRMag_here;
        
        % coriolises 1
        eta_here                         =eta_prec+intdAdjgn_here*xcidotn;
        Co1_here                         =dinamico_coadj(eta_here)*Mag_here;
        trapz                            =dX(jj)*(Co1_prec+Co1_here)/2;
        Co1X(:,6*(ii-1)+1:6*nsez)        =Co1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Co1_prec                         =Co1_here;
        LCo1_here                        =intdAdjgn_here'*Co1_here;
        trapz                            =dX(jj)*(LCo1_prec+LCo1_here)/2;
        LCo1X(:,6*(ii-1)+1:6*nsez)       =LCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LCo1_prec                        =LCo1_here;
        RCo1_here                        =Co1_here*intdAdjgn_here;
        trapz                            =dX(jj)*(RCo1_prec+RCo1_here)/2;
        RCo1X(:,6*(ii-1)+1:6*nsez)       =RCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        RCo1_prec                        =RCo1_here;
        LRCo1_here                       =intdAdjgn_here'*Co1_here*intdAdjgn_here;
        trapz                            =dX(jj)*(LRCo1_prec+LRCo1_here)/2;
        LRCo1X(:,6*(ii-1)+1:6*nsez)      =LRCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LRCo1_prec                       =LRCo1_here;
        
        % coriolises 2
        Co2_here                         =Mag_here*dinamico_adj(intdAdjgn_here*xcidotn);
        trapz                            =dX(jj)*(Co2_prec+Co2_here)/2;
        Co2X(:,6*(ii-1)+1:6*nsez)        =Co2X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Co2_prec                         =Co2_here;
        LCo2_here                        =intdAdjgn_here'*Co2_here;
        trapz                            =dX(jj)*(LCo2_prec+LCo2_here)/2;
        LCo2X(:,6*(ii-1)+1:6*nsez)       =LCo2X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LCo2_prec                        =LCo2_here;
        
        % drag
        Dra_here                         =coAdjgn_here*D(:,(jj-1)*6+1:(jj-1)*6+6)*invAdjgn_here*norm(eta_here(4:6));
        trapz                            =dX(jj)*(Dra_prec+Dra_here)/2;
        DraX(:,6*(ii-1)+1:6*nsez)        =DraX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Dra_prec                         =Dra_here;
        LDra_here                        =intdAdjgn_here'*Dra_here;
        trapz                            =dX(jj)*(LDra_prec+LDra_here)/2;
        LDraX(:,6*(ii-1)+1:6*nsez)       =LDraX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LDra_prec                        =LDra_here;
        RDra_here                        =Dra_here*intdAdjgn_here;
        trapz                            =dX(jj)*(RDra_prec+RDra_here)/2;
        RDraX(:,6*(ii-1)+1:6*nsez)       =RDraX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        RDra_prec                        =RDra_here;
        LRDra_here                       =intdAdjgn_here'*Dra_here*intdAdjgn_here;
        trapz                            =dX(jj)*(LRDra_prec+LRDra_here)/2;
        LRDraX(:,6*(ii-1)+1:6*nsez)      =LRDraX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LRDra_prec                       =LRDra_here;
    end
    MasX            =MasX-repmat(MasX(:,1:6),[1,nsez]);
    LMasX           =LMasX-repmat(LMasX(:,1:6),[1,nsez]);
    MagX            =MagX-repmat(MagX(:,1:6),[1,nsez]);
    LMagX           =LMagX-repmat(LMagX(:,1:6),[1,nsez]);
    RMagX           =RMagX-repmat(RMagX(:,1:6),[1,nsez]);
    LRMagX          =LRMagX-repmat(LRMagX(:,1:6),[1,nsez]);
    Co1X            =Co1X-repmat(Co1X(:,1:6),[1,nsez]);
    LCo1X           =LCo1X-repmat(LCo1X(:,1:6),[1,nsez]);
    RCo1X           =RCo1X-repmat(RCo1X(:,1:6),[1,nsez]);
    LRCo1X          =LRCo1X-repmat(LRCo1X(:,1:6),[1,nsez]);
    Co2X            =Co2X-repmat(Co2X(:,1:6),[1,nsez]);
    LCo2X           =LCo2X-repmat(Co2X(:,1:6),[1,nsez]);
    DraX            =DraX-repmat(DraX(:,1:6),[1,nsez]);
    LDraX           =LDraX-repmat(LDraX(:,1:6),[1,nsez]);
    RDraX           =RDraX-repmat(RDraX(:,1:6),[1,nsez]);
    LRDraX          =LRDraX-repmat(LRDraX(:,1:6),[1,nsez]);
    Mas             =MasX(:,6*(nsez-1)+1:6*nsez);
    LMas            =LMasX(:,6*(nsez-1)+1:6*nsez);
    Mag             =MagX(:,6*(nsez-1)+1:6*nsez);
    LMag            =LMagX(:,6*(nsez-1)+1:6*nsez);
    RMag            =RMagX(:,6*(nsez-1)+1:6*nsez);
    LRMag           =LRMagX(:,6*(nsez-1)+1:6*nsez);
    Co1             =Co1X(:,6*(nsez-1)+1:6*nsez);
    LCo1            =LCo1X(:,6*(nsez-1)+1:6*nsez);
    RCo1            =RCo1X(:,6*(nsez-1)+1:6*nsez);
    LRCo1           =LRCo1X(:,6*(nsez-1)+1:6*nsez);
    Co2             =Co2X(:,6*(nsez-1)+1:6*nsez);
    LCo2            =LCo2X(:,6*(nsez-1)+1:6*nsez);
    Dra             =DraX(:,6*(nsez-1)+1:6*nsez);
    LDra            =LDraX(:,6*(nsez-1)+1:6*nsez);
    RDra            =RDraX(:,6*(nsez-1)+1:6*nsez);
    LRDra           =LRDraX(:,6*(nsez-1)+1:6*nsez);
    
    % Actuation and internal load
  %  if t<=tact                                      % virata
       % Fan         =[Famx(:,jj) Famy(:,jj) Famz(:,jj) Fax(:,jj) Fay(:,jj) Faz(:,jj)]*t/tact;
       % Fan          =interp1(t_tens,Fan,t)';
    %end
%     if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
         Fan         =[Famx(:,jj) Famy(:,jj) Famz(:,jj) Fax(:,jj) Fay(:,jj) Faz(:,jj)];
         Fan          =interp1(t_tens,Fan,t)';
%     end
    Fin             =Eps(:,(jj-1)*6+1:(jj-1)*6+6)*(xcin-xci_star)+Ipsi(:,(jj-1)*6+1:(jj-1)*6+6)*xcidotn;
    Fpn             =[Fpmx(jj);Fpmy(jj);Fpmz(jj);Fpx(jj);Fpy(jj);Fpz(jj)];
    
    % Actuation and internal load successiva
    if jj~=npie
        %if t<=tact                                      % virata
         %   Fan_suc     =[Famx(:,jj+1) Famy(:,jj+1) Famz(:,jj+1) Fax(:,jj+1) Fay(:,jj+1) Faz(:,jj+1)]*t/tact;
          %  Fan_suc     =interp1(t_tens,Fan_suc,t)';
        %end
        %if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
            Fan_suc     =[Famx(:,jj+1) Famy(:,jj+1) Famz(:,jj+1) Fax(:,jj+1) Fay(:,jj+1) Faz(:,jj+1)];
            Fan_suc     =interp1(t_tens,Fan_suc,t)';
        %end
    else
        Fan_suc     =[0;0;0;0;0;0];
    end
    
    % aggiorno coefficienti dinamica
    invAdjgn_last   =piecewise_invAdjoint(X((jj-1)*nsez+nsez),thetan,xcin);
    invAdjgnR_last  =blkdiag(invAdjgn_last(1:3,1:3),invAdjgn_last(4:6,4:6));
    invAdjgprec     =dinamico_Adjoint(g_prec^-1);
    invAdjgprecR    =blkdiag(invAdjgprec(1:3,1:3),invAdjgprec(4:6,4:6));
    intdAdjgn_last  =piecewise_ADJ(X((jj-1)*nsez+nsez),thetan,xcin);
    MagB            =blkdiag([repmat(Mag,[jj-1 jj-1]) repmat(RMag,[jj-1 1]); repmat(LMag,[1 jj-1]) LRMag],...
                     zeros(6*(npie-jj),6*(npie-jj)));
    genMagM         =genMagM+Jaco_prec'*MagB*Jaco_prec;
    Co1B            =blkdiag([repmat(Co1,[jj-1 jj-1]) repmat(RCo1,[jj-1 1]); repmat(LCo1,[1 jj-1]) LRCo1],...
                     zeros(6*(npie-jj),6*(npie-jj)));
    genCo1M         =genCo1M+Jaco_prec'*Co1B*Jaco_prec;
    Co2B            =blkdiag([repmat(Co2,[jj-1 jj-1]) zeros(6*(jj-1),6); repmat(LCo2,[1 jj-1]) zeros(6,6)],...
                     zeros(6*(npie-jj),6*(npie-jj)))+MagB*adetan_prec;
    genCo2M         =genCo2M+Jaco_prec'*Co2B*Jaco_prec;
    GraB            =[repmat(Mas,[jj-1 1]);LMas;zeros(6*(npie-jj),6)];
    genGraV         =genGraV+Jaco_prec'*GraB*dinamico_Adjoint(g_prec^-1);
    DraB            =blkdiag([repmat(Dra,[jj-1 jj-1]) repmat(RDra,[jj-1 1]); repmat(LDra,[1 jj-1]) LRDra],...
                     zeros(6*(npie-jj),6*(npie-jj)));
    genDraM         =genDraM+Jaco_prec'*DraB*Jaco_prec;
%     ForB            =[repmat(invAdjgn_last'*(Fan-Fan_suc),[jj-1 1]);(invAdjgn_last*intdAdjgn_last)'*(Fan-Fan_suc);...
%                       zeros(6*(npie-jj),1)];                                                 % pneumatico
%     genForV         =genForV+Jaco_prec'*ForB-[zeros(6*(jj-1),1);L*Fin;zeros(6*(npie-jj),1)]; % pneumatico
    genForV         =genForV+[repmat(L(jj)*Fan,[jj-1 1]);L(jj)*(Fan-Fin);zeros(6*(npie-jj),1)];       % cavo tip2base
    TipB            =[repmat(invAdjgn_last'*(invAdjgnR_last*invAdjgprecR*Fpn),[jj-1 1]);...
                     (invAdjgn_last*intdAdjgn_last)'*(invAdjgnR_last*invAdjgprecR*Fpn);zeros(6*(npie-jj),1)];
    genTipV         =genTipV+Jaco_prec'*TipB;
    
    %----------------------------------------------------------------------
    % recursive factors
    prec2prec       =[];
    for ii=1:jj-1
        prec2prec   =blkdiag(prec2prec,invAdjgn_last);
    end
    Jaco_prec       =blkdiag(prec2prec,invAdjgn_last*intdAdjgn_last,zeros(6*(npie-jj),6*(npie-jj)))*Jaco_prec+...
                     blkdiag(zeros(6*jj,6*jj),repmat(eye(6),[jj~=npie jj~=npie]),zeros(6*(npie-jj-1),6*(npie-jj-1)));
    g_prec          =g_prec*piecewise_expmap(X((jj-1)*nsez+nsez),thetan,xcin);
    ADxin           =intdAdjgn_last*xcidotn;
    eta_prec        =invAdjgn_last*(eta_prec+ADxin);
    prec2prec       =[];
    prec2prec_inv   =[];
    for ii=1:npie
        prec2prec     =blkdiag(prec2prec,invAdjgn_last);
        prec2prec_inv =blkdiag(prec2prec_inv,invAdjgn_last^-1);
        for zz=1:npie
            if (ii+zz == jj)
                adetan_prec(6*(ii-1)+1:6*(ii-1)+6,6*(zz-1)+1:6*(zz-1)+6)  =dinamico_adj(ADxin);
            end
        end
    end
    adetan_prec     =prec2prec*adetan_prec*prec2prec_inv;
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% calcolo le derivate

Xci_punto       =Xcidot;
Xcidot_punto    =genMagM^-1*(genForV+(1-ro_water/ro_arm)*genGraV*dinamico_Adjoint(g_r^-1)*Gra+genTipV-(genCo1M-genCo2M+genDraM)*Xcidot);

%----------------------------------------------------------------------
z_punto         =[Xci_punto;Xcidot_punto];
              
% eof