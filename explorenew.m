time=2;
rng(1111);
nsol=time*100+1;
 load('ini_cond.mat')
%  ini(1:3,:)=unifrnd(-0.5,0.5,3,nsol);
%  ini(4:6,:)=unifrnd(-0.05,0.05,3,nsol);
%  ini=ones(6,nsol)*0.05;
% 
% ini=filter(ones(1,1000),1 , ini);
% sc
% ini=ini*10;
%global net


ini=zeros(3,nsol);
%init=abs(randn(3,1))*5;
init=ones(3,1)*10;
eps=1;
for i=2:50:nsol

    for j=1:3
            f1=rand;
    f2=rand;
%         if abs(ini(j,i-1)-init(j,1))<eps
%             init(j,1)=-abs(randn)*5;
%             ini(j,i)=ini(j,i-1);
%         else
%             ini(j,i)=ini(j,i-1)+sign((-ini(j,i-1)+init(j,1)))*rand*eps;
%         end
%         if f1>0.95
%             ini(j,i)=-abs(randn)*0;
%         end
        if i>100&&f2>0.2
            ini(j,i-49:i)=0;
        elseif i>50
           % ini(j,i-49:i)=abs(randn)*5;
             ini(j,i-49:i)=10;
        else
            ini(j,i:i+100)=0;
        end
        if f1>0.9&&i>50
           % ini(j,i-29:i)=abs(randn)*5;
             ini(j,i-29:i)=10;
        end
    end
    
end
 %ini(1:3,1:10)=0;
%ini(3,1:200)=-1;

ini(:,:)=zeros(3,nsol);
%ini(1:6,:)=zeros(6,nsol);
%ini(1,:)=zeros(1,nsol);
ini(1,1:end)=10;

inipos = piecewise_driver(ini,time,ini_cond);

save chaosairxx.mat


endpos=squeeze(inipos(:,end,:));
for i=1:nsol-1
    
    scatter3(inipos(1,:,i,4),inipos(2,:,i,4),inipos(3,:,i,4),'r')
    pause(0.1)
    hold on
    if mod(i,10)==0
        clear figure;
    end
end


% time=0.1;
% ini(:,1:2)=zeros(6,2);
% 
% 
% 
% 
% 
% inipos = piecewise_driver(ini',time,net);
% 
% 
% 
% 
% 
% for par=2:100
%        target(par-1,:)=[-0.25+0.5*sin((par-2)/50);0;0];
%     targetv(par-1,:)=[-0.25+0.5*sin((par-2)/50)+0.01;0;0]; 
% end
