
function rew= evalu_4sec(inp,para)
x=para.x;
%net=para.net;
netc=para.netc;

target=para.tar;

 
 inp2(:,1)=inp(:,1);
for i=2:para.stp
    if mod(i,para.cc)==0
        inp2(:,i)=inp(:,fix(i/para.cc)+1);
    else
        inp2(:,i)=inp2(:,i-1);
    end
end
 inp2(:,para.stp+1)=0;
 
  siz=length(inp2);
 
 X = tonndata(inp2,true,false);
 T=para.T;
%T = tonndata(t,true,false);
[xc,xic,aic] = preparets(netc,X,{},T(1,1:siz));
yc = netc(xc,xic,aic);
yc=cell2mat(yc);

%effort=max(rssq(inp2));
effort=sum(rssq(inp2))+1;
velo=rssq(diff(yc(10:12,:)',1));
jerk=diff(velo,2);
%minjerk=max(jerk);
minjerk=sum(jerk);
yc=(yc(:,end));

%yc=cell2mat(yc(end));
%yc=cell2mat(yc);
etf=length(yc);
%target=repmat(target,1,etf);

%[err,Ind]=min(rssq((yc(10:12,:)-target))*10);
err=(rssq((yc(10:12,:)-target))*10);
%err=sum(rssq(yc(4:6,:)-target(:,1:end-1))*10);
%rew=err^2+errvec;
%energ=1/((test(1,i+2)+0.25)*10+0.5*rssq(test(:,i+2)-test(:,i+1)));


rew=err;%+(effort/100);%(effort/100);%minjerk/100;%
%rew=err+(Ind/1000);%+(Ind/10000)^2;
%rew=err;
% 
% test(:,1)=net(x(:,1));
% test(:,2)=net(x(:,2));
% %target=[0.0495555303072983;0.0202573482751266;0.0914981560750420];
% %target=[-0.25;0.1;0];
% %tarvec=target-test(1:3,1);
% 
% for i=1:para.stp
%     
%     % test(:,i+2)=net([test(1:3,i);test(1:3,i+1);inp(:,i);test(4:6,i);test(4:6,i+1)]);
%     if mod(i,para.cc)==0
%         test(:,i+2)=net([test(1:6,i);test(1:6,i+1);inp(:,fix(i/para.cc)+1);test(7:12,i);test(7:12,i+1)]);
%     else
%         test(:,i+2)=net([test(1:6,i);test(1:6,i+1);zeros(6,1);test(7:12,i);test(7:12,i+1)]);
%     end
%     
%     %test(:,i+2)=net([test(1:3,i);test(1:3,i+1);ini(:,i+3)]);
%     % test(:,i+2)=net([test(1:3,i);test(1:3,i+1);ini(:,i+2);x(13:end,i+2)]);
%     % test(:,i+2)=net([test(1:3,i);test(1:3,i+1);ini(:,i+3);test(4:6,i);test(4:6,i+1);test(7:9,i);test(7:9,i+1)]);
%     %target(:,i)=[-0.25+i/(10*para.stp);i/(10*para.stp);0];
%     target(:,i)=para.tar;
%     cerr(i)=rssq(test(4:6,i+2)-target(:,i))*i^2;
% end
% %actvec=target-test(1:3,para.stp+2);
% 
% %errvec=abs(atan2d(norm(cross(tarvec,actvec)),dot(tarvec,actvec)))/50;
% err=rssq(test(4:6,i+2)-target(:,i))*10;
% %rew=err^2+errvec;
% energ=1/((test(1,i+2)+0.25)*10+0.5*rssq(test(:,i+2)-test(:,i+1)));
% rew=err+0*energ;
% %rew=sum(cerr);

end





