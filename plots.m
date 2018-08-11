kl=5000;
R1       =13.8e-3*kl;                       % [m] Raggio sezione 1
R2       =11.1e-3*kl;                       % [m] Raggio sezione 2
R3       =8.2e-3*kl;                        % [m] Raggio sezione 3
R4       =5.4e-3*kl;                        % [m] Raggio sezione 4

figure;
for i=1:100
    scatter3(target(1,:),target(2,:),target(3,:),'g','filled')
    h=scatter3(pos(1,:,i,1),pos(2,:,i,1),pos(3,:,i,1),R1,'r');
    hold on
    scatter3(target(1,:),target(2,:),target(3,:),'g','filled')
    h1=scatter3(pos(1,:,i,2),pos(2,:,i,2),pos(3,:,i,2),R2,'b');
    h2=scatter3(pos(1,:,i,3),pos(2,:,i,3),pos(3,:,i,3),R3,'r');
    h3=scatter3(pos(1,:,i,4),pos(2,:,i,4),pos(3,:,i,4),R4,'b');
    
    scatter3(target2(1,:),target2(2,:),target2(3,:),'b','filled')
    hh=scatter3(pos2(1,:,i,1),pos2(2,:,i,1),pos2(3,:,i,1),R1,'r');
    hold on
    hh1=scatter3(pos2(1,:,i,2),pos2(2,:,i,2),pos2(3,:,i,2),R2,'b');
    hh2=scatter3(pos2(1,:,i,3),pos2(2,:,i,3),pos2(3,:,i,3),R3,'r');
    
    hh3=scatter3(pos2(1,:,i,4),pos2(2,:,i,4),pos2(3,:,i,4),R4,'b');
    
    scatter3(target3(1,:),target3(2,:),target3(3,:),'r','filled')
    hhh=scatter3(pos3(1,:,i,1),pos3(2,:,i,1),pos3(3,:,i,1),R1,'r');
    hold on
    hhh1=scatter3(pos3(1,:,i,2),pos3(2,:,i,2),pos3(3,:,i,2),R2,'b');
    hhh2=scatter3(pos3(1,:,i,3),pos3(2,:,i,3),pos3(3,:,i,3),R3,'r');
    
    hhh3=scatter3(pos3(1,:,i,4),pos3(2,:,i,4),pos3(3,:,i,4),R4,'b');
    %view([0 90] );
    pause(0.5)
    hold on
        xlim([-0.4 0.2])
        ylim([-0.1 0.5])
        zlim([ -0.15 0.15])
    delete(h);
    delete(h1);
    delete(h2);
    delete(h3);
    
        delete(hh);
    delete(hh1);
    delete(hh2);
    delete(hh3);
    
        delete(hhh);
    delete(hhh1);
    delete(hhh2);
    delete(hhh3);
    %axis([-0.1 0.1 -0.5 0.5 -1 1])
    % clf
end