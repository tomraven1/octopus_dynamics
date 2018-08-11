tic
para.stp=100;
para.cc=2;
rng(200)
%inp1=zeros(6,para.stp);
for kk=1:20
    inp1=zeros(3,fix(para.stp/para.cc)+1);
    ang=1:1:para.stp;
    pick=randi(7000);
   % targetact=[-0.255308042020190;0.283693370769576;0.0800899533889185];
    targetact=endpos(10:12,pick);
    para.tar=targetact;
    
    A=2*ones(1,3*(para.stp+1));
    b=1;
    lb=zeros(1,6*(para.stp+1));
    ub=3*ones(1,3*(para.stp+1));
    f = @(inp)evalu_4sec(inp,para);
    options =  optimoptions (@fmincon,'Display','iter','MaxIter',40,'algorithm','sqp','UseParallel',true);
    %options =  optimset ('Display','iter');
    [inp,fval] = fmincon(f,inp1,[],[],[],[],lb,ub,[],options)%fminsearch%fminunc
    
    
    %[inp,fval] = ga(f,inp1,[],[],[],[],lb,ub,[],options)
    toc
    
    tim=(para.stp)/100;
    inp2(:,1)=inp(:,1);
    for i=2:para.stp
        if mod(i,para.cc)==0
            inp2(:,i)=inp(:,fix(i/para.cc)+1);
        else
            inp2(:,i)=inp2(:,i-1);
        end
    end
    inp2(:,para.stp+1:para.stp+21)=0;
    
    tim=tim+0.20;
    %
    % inp3(1,:,:)=inp2(1:4,:);
    % inp3(2,:,:)=inp2(5:8,:);
    
    %inp2(3,:)=0;
    
   pos = piecewise_driver3(inp2,tim,ini_cond);
    
    % test(:,1)=net(x(:,1));
    %  test(:,2)=net(x(:,2));
    target=targetact;
    %
    % for i=1:para.stp
    %
    %     test(:,i+2)=net([test(1:6,i);test(1:6,i+1);inp2(:,i);test(7:12,i);test(7:12,i+1)]);
    %
    % end
    
    siz=length(inp2);
    
   % inp2(3,:)=[];
    X = tonndata(inp2,true,false);
    T=para.T;
    %T = tonndata(t,true,false);
    [xc,xic,aic] = preparets(netc,X,{},T(1,1:siz));
    yc = netc(xc,xic,aic);
    test=cell2mat(yc);
    
    
    etf=length(test);
    target1=repmat(target,1,etf);
    
    [err,Ind]=min(rssq((test(10:12,:)-target1)));
    in(kk)=Ind;
    errall(kk)=err;
    
    testall(:,:,kk)=test;
    posall(:,:,:,:,kk)=pos;
    targetall(:,kk)=target;
end

save asd.mat

scatter3(test(10,1:para.stp-1),test(11,1:para.stp-1),test(12,1:para.stp-1))
hold on
scatter3(pos(1,end,1:end,4),pos(2,end,1:end,4),pos(3,end,1:end,4))
% hold on
scatter3(target(1,:),target(2,:),target(3,:),'g','filled')
% scatter3(endpos(10,:),endpos(11,:),endpos(12,:),'.')
% scatter3(endpos(10,1),endpos(11,1),endpos(12,1),'k')
toc
% for i=1:100
%
%     scatter3(inipos(1,:,i),inipos(2,:,i),inipos(3,:,i),'r')
%     pause(0.1)
%     hold on
%     if mod(i,10)==0
%         clear figure;
%     end
% end


%
% h=figure;
% axis([-0.1 0.1 -0.5 0.5 -1 1])
%
%
kl=5000;
R1       =13.8e-3*kl;                       % [m] Raggio sezione 1
R2       =11.1e-3*kl;                       % [m] Raggio sezione 2
R3       =8.2e-3*kl;                        % [m] Raggio sezione 3
R4       =5.4e-3*kl;                        % [m] Raggio sezione 4

figure;
for i=1:600
   % scatter3(target(1,:),target(2,:),target(3,:),'g','filled')
     hold on
    h=scatter3(pos(1,:,i,1),pos(2,:,i,1),pos(3,:,i,1),R1,'r');
    hold on
    h1=scatter3(pos(1,:,i,2),pos(2,:,i,2),pos(3,:,i,2),R2,'b');
    h2=scatter3(pos(1,:,i,3),pos(2,:,i,3),pos(3,:,i,3),R3,'r');
    
    h3=scatter3(pos(1,:,i,4),pos(2,:,i,4),pos(3,:,i,4),R4,'b');
    %view([0 90] );
    pause(0.2)
    hold on
         xlim([-0.3 0.1])
        ylim([-0.1 0.5])
         zlim([ -0.02 0.15])
    delete(h);
    delete(h1);
    delete(h2);
    delete(h3);
    %axis([-0.1 0.1 -0.5 0.5 -1 1])
    % clf
end
%   h=scatter3(pos(1,:,i,1),pos(2,:,i,1),pos(3,:,i,1),R1,'r');
%     hold on
%     h1=scatter3(pos(1,:,i,2),pos(2,:,i,2),pos(3,:,i,2),R2,'b');
%    h2=scatter3(pos(1,:,i,3),pos(2,:,i,3),pos(3,:,i,3),R3,'r');
%
%      h3=scatter3(pos(1,:,i,4),pos(2,:,i,4),pos(3,:,i,4),R4,'b');
%
%
%    plot(  rssq(squeeze(pos(:,end,2:end-1,4))-squeeze(pos(:,end,1:end-2,4))))
%
%       plot( (squeeze(pos(1,end,2:end-1,4))-squeeze(pos(1,end,1:end-2,4))))
%
%



for i=1:20
    if in(i)>31
        asd=cross((squeeze(posall(:,end,1:in(i)-2,4,i)-posall(:,end,2:in(i)-1,4,i)))',(squeeze(posall(:,end,2:in(i)-1,4,i)-posall(:,end-1,2:in(i)-1,4,i)))')';
        %asd=cross((squeeze(posall(:,end,1:end-2,4,i)-posall(:,end,2:end-1,4,i)))',(squeeze(posall(:,end,2:end-1,4,i)-posall(:,end-1,2:end-1,4,i)))')';
        
        % asd=cross((squeeze(testall(10:12,1:in(i)-2,i)-testall(10:12,2:in(i)-1,i)))',(squeeze(testall(10:12,2:in(i)-1,i)-testall(7:9,2:in(i)-1,i)))')';
        
        %plot(  rssq(asd))
        [juk,ind]=max(rssq(asd));
        ll=length(asd)-1;
        plot( (-ind:ll-ind) ,rssq(asd)/max(rssq(asd)),'k')
%         
       %  wtf=rssq(asd)/max(rssq(asd));
        % hop(:,i)=wtf(ind-30:ind+20);
        %/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
        hold on
        %plot(  rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i))))%/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
        hold on
    end
end



clear asd;
for jj=2:22
    
    for ij=1:4
        for i=3
            if in(i)>10
                asd(:,:,ij)=cross((squeeze(posall(:,jj,1:in(i)-2,ij,i)-posall(:,jj,2:in(i)-1,ij,i)))',(squeeze(posall(:,jj,2:in(i)-1,ij,i)-posall(:,jj-1,2:in(i)-1,ij,i)))')';
                hold on
                
                % asd=cross((squeeze(testall(10:12,1:in(i)-2,i)-testall(10:12,2:in(i)-1,i)))',(squeeze(testall(10:12,2:in(i)-1,i)-testall(7:9,2:in(i)-1,i)))')';
                
                % plot(  rssq(asd))
                % [juk,ind]=max(rssq(asd));
                %  ll=length(asd)-1;
                % plot( (-ind:ll-ind) ,rssq(asd)/max(rssq(asd)))
                
                %plot(  rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i))))%/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
                
            end
        end
    end
    ll=length(asd);
    asda=reshape(asd,[3,ll*4]);
    plot( rssq(asda)/max(rssq(asda)))
    %/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
    hold on
    
    
end

% asd=[];
% for i=1:22
%     for k=1:4
%         asd=[asd ; squeeze(inipos(:,i,:,k));];
%     end
% end


for i=1:20
    pos=posall(:,:,:,:,i);
    target=targetall(:,i);
    
    target1=repmat(target,1,121);
    
    
    [err,Ind]=min(rssq((squeeze(pos(:,end,:,4))-target1)));
    
    errall(i)=err;
    in(i)=Ind;
end
