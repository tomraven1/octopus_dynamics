asd=size(inipos,3);
new=7000;
endpos=squeeze(inipos(:,end,1:end,:));
endpos=cat(1,endpos(:,:,1),(endpos(:,:,2)),endpos(:,:,3),endpos(:,:,4));
%endpos=cat(1,endpos(:,:,1),endpos(:,:,4));
%endpos=reshape(endpos,[],asd);
%endori=squeeze(inipos(:,end-1,:));
endori=squeeze(inipos(:,end-1,1:end,:));
endori=cat(1,endori(:,:,1),(endori(:,:,2)),endori(:,:,3),endori(:,:,4));
%endori=cat(1,endori(:,:,1),endori(:,:,4));
%endori=reshape(endori,[],asd);
endori2=squeeze(inipos(:,end-2,1:end,:));
endori2=cat(1,endori2(:,:,1),endori2(:,:,2));
%endori2=reshape(endori2,[],asd);

inp=ini;
% Fax         =[(-inp(1,:)-inp(2,:)-inp(3,:)-inp(4,:)-inp(5,:)-inp(6,:));
%     (-inp(4,:)-inp(5,:)-inp(6,:))  ] ;             % [N] contraction load
% 
% Famy        =[((inp(1,:)+inp(4,:))*0.01-(inp(2,:)+inp(5,:))*0.01*cosd(60)-(inp(3,:)+inp(6,:))*0.01*cosd(60));
%     (0.01*inp(4,:)-0.01*inp(5,:)*0.5-0.01*inp(6,:)*0.5)];              % [Nm] bending torque
% Famz        =[(-(inp(2,:)+inp(5,:))*0.01*cosd(30)+(inp(3,:)+inp(6,:))*0.01*cosd(30));
%     (-0.01*inp(5,:)*cosd(30)+0.01*inp(6,:)*cosd(30))];
% 
% FF=[Fax;Famy;Famz];
% FF=FF(:,1:new+1);

ini=ini(:,1:new);
endpos=endpos(:,1:new);
endori=endori(:,1:new);
endori2=endori2(:,1:new);
% endpos=endpos-repmat(endpos(:,1),1,asd);
% endori=endori-repmat(endori(:,1),1,asd);
% endori2=endori2-repmat(endori2(:,1),1,asd);



%ix=[endpos(:,2:end-2)-endpos(:,1:end-3);endpos(:,2:end-2) ;endpos(:,3:end-1)-endpos(:,2:end-2)];
%ix2=[endori(:,1:end-3);endori(:,2:end-2) ;endori(:,3:end-1)];
%ix3=[endori2(:,1:end-3);endori2(:,2:end-2) ;endori2(:,3:end-1)];
%x = [endpos(:,1:end-2);endpos(:,2:end-1) ;endpos(:,3:end);endori(:,1:end-2);endori(:,2:end-1) ;endori(:,3:end) ];

x=[ini(:,2:end)];
%x=[FF(:,3:end-1)];
%x=[endpos(:,1:end-3);endpos(:,2:end-2) ;ini(:,2:end-3);endori(:,1:end-3);endori(:,2:end-2) ;endori2(:,1:end-3);endori2(:,2:end-2)];
%x=[endpos(:,1:end-3);endpos(:,2:end-2) ;ini(:,2:end-3)];

%x=[ix;ix2;ix3];

%x = [inipos(:,1:end-3);  inipos(:,2:end-2)  ;inipos(:,3:end-1);inipos(:,4:end);ini(:,1:end-4)  ];

%t=ini2(:,2:end-3);
%t=[endpos(:,2:end)];
t=[endpos(:,2:end);endori(:,2:end)];%;(endpos(:,3:end)-endpos(:,2:end-1));(endori(:,3:end)-endori(:,2:end-1))];
%
% x=mapminmax(x,-1,1);
% t=mapminmax(t,-1,1);


X = tonndata(x,true,false);
T = tonndata(t,true,false);

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainbr';  % Levenberg-Marquardt backpropagation.

% Create a Nonlinear Autoregressive Network with External Input
inputDelays = 1;
feedbackDelays = 1:2;
hiddenLayerSize = [40];
%net = layrecnet(1:2,[20 15 ]);
net = narxnet(inputDelays,feedbackDelays,hiddenLayerSize,'open',trainFcn);
net.trainParam.epochs=6;
% Choose Input and Feedback Pre/Post-Processing Functions
% Settings for feedback input are automatically applied to feedback output
% For a list of all processing functions type: help nnprocess
% Customize input parameters at: net.inputs{i}.processParam
% Customize output parameters at: net.outputs{i}.processParam
net.inputs{1}.processFcns = {'removeconstantrows','mapminmax'};
net.inputs{2}.processFcns = {'removeconstantrows','mapminmax'};

% Prepare the Data for Training and Simulation
% The function PREPARETS prepares timeseries data for a particular network,
% shifting time by the minimum amount to fill input states and layer
% states. Using PREPARETS allows you to keep your original time series data
% unchanged, while easily customizing it for networks with differing
% numbers of delays, with open loop or closed loop feedback modes.
[x,xi,ai,t] = preparets(net,X,{},T);

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivide
net.divideFcn = 'divideblock';  % Divide data randomly
net.divideMode = 'time';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate', 'ploterrhist', ...
    'plotregression', 'plotresponse', 'ploterrcorr', 'plotinerrcorr'};

%net.layers{2}.transferFcn = 'satlins';
% Train the Network
[net,tr] = train(net,x,t,xi,ai);

% Test the Network
y = net(x,xi,ai);
e = gsubtract(t,y);
performance = perform(net,t,y)

% Recalculate Training, Validation and Test Performance
trainTargets = gmultiply(t,tr.trainMask);
valTargets = gmultiply(t,tr.valMask);
testTargets = gmultiply(t,tr.testMask);
trainPerformance = perform(net,trainTargets,y)
valPerformance = perform(net,valTargets,y)
testPerformance = perform(net,testTargets,y)

% View the Network
%view(net)



% Closed Loop Network
% Use this network to do multi-step prediction.
% The function CLOSELOOP replaces the feedback input with a direct
% connection from the outout layer.
netc = closeloop(net);
%[netc,xi,ai] = closeloop(netc,xi,ai);
netc.name = [net.name ' - Closed Loop'];
%view(netc)
[xc,xic,aic,tc] = preparets(netc,X,{},T);
yc = netc(xc,xic,aic);

%netc.performFcn = 'mse';


%netc.performParam.normalization = 'standard';
%netc.inputs{1}.processFcns = {'removeconstantrows','mapminmax'};
%netc.inputs{2}.processFcns = {'removeconstantrows','mapminmax'};


closedLoopPerformance = perform(net,tc,yc)
netc.trainParam.epochs=1000;
% netc.divideFcn = 'divideblock';  % Divide data randomly
% netc.divideMode = 'time';  % Divide up every sample
% netc.divideParam.trainRatio = 80/100;
% netc.divideParam.valRatio = 5/100;
% netc.divideParam.testRatio = 15/100;
 netc.trainFcn = 'trainbr';
% netc.trainParam.mu_max=10000000000000000000;


[netc,tr] = train(netc,xc,tc,xic,aic);
% netcn=netc;
% for i=1:50
% hope=bttderiv('de_dwb',netcn,xc,tc,xic,aic);
% wb=getwb(netcn);
% newwb=wb+0.0000001*hope;
% netcn=setwb(netcn,newwb);
% yc = netcn(xc,xic,aic);
% asfgs(i) = perform(netcn,tc,yc)
% end

yc = netc(xc,xic,aic);
closedLoopPerformance = perform(net,tc,yc)



asd=cell2mat(yc);
asd2=cell2mat(t);
asd3=cell2mat(y);
% asd4=cell2mat(tc);
% 
plot(asd2(4,:))
hold on
plot(asd(4,:),'r')
% hold on
% plot(asd3(4,:))

% for i=1:8
%     numTimesteps = size(x,2)-i*100;
%     knownOutputTimesteps = 1:(numTimesteps-100);
%     predictOutputTimesteps = (numTimesteps-99):numTimesteps;
%     X1 = X(:,knownOutputTimesteps);
%     T1 = T(:,knownOutputTimesteps);
%     [x1,xio,aio] = preparets(net,X1,{},T1);
%     [y1,xfo,afo] = net(x1,xio,aio);
%     % Next the the network and its final states will be converted to
%     % closed-loop form to make five predictions with only the five inputs
%     % provided.
%     x2 = X(1,predictOutputTimesteps);
%     [netc,xic,aic] = closeloop(net,xfo,afo);
%     [y2,xfc,afc] = netc(x2,xic,aic);
%     multiStepPerformance(i) = perform(net,T(1,predictOutputTimesteps),y2);
%     
%     asd(:,100*i-99:i*100)=cell2mat(T(1,predictOutputTimesteps));
%     asd2(:,100*i-99:i*100)=cell2mat(y2);
% end
% 
% plot(asd2(4,:))
% hold on
% plot(asd(4,:))
