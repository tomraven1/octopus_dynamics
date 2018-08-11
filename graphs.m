lala=length(Pos);

for i=1:20
    pos=posall(:,:,:,:,i);
    target=targetall(:,i);
    
    target1=repmat(target,1,lala);
    
    
    [err,Ind]=min(rssq((squeeze(pos(:,end,:,4))-target1)));
    
    errall(i)=err;
    in(i)=Ind;
    dist(i)=rssq(squeeze(pos(:,end,1,4))-target1(:,1));
end



for i=1:20
    if in(i)>51
        %asd=cross((squeeze(posall(:,end,1:in(i)-2,4,i)-posall(:,end,2:in(i)-1,4,i)))',(squeeze(posall(:,end,2:in(i)-1,4,i)-posall(:,end-1,2:in(i)-1,4,i)))')';
        asd=cross((squeeze(posall(:,end,1:end-2,2,i)-posall(:,end,2:end-1,2,i)))',(squeeze(posall(:,end,2:end-1,2,i)-posall(:,end-1,2:end-1,2,i)))')';
          
       
       
       %asd=cross((squeeze(Posiall(1:end-2,:,1,1,i)-Posiall(2:end-1,:,1,1,i)))',(squeeze(Posiall(2:end-1,:,1,1,i)-Posiall(2:end-1,:,3,1,i)))')';
      
           
        %asd=cross((squeeze(posact(i,:,end,1:end-2,2)-posact(i,:,end,2:end-1,2)))',(squeeze(posact(i,:,end,2:end-1,2)-posact(i,:,end-1,2:end-1,2)))')';
       
        % asd=cross((squeeze(testall(10:12,1:in(i)-2,i)-testall(10:12,2:in(i)-1,i)))',(squeeze(testall(10:12,2:in(i)-1,i)-testall(7:9,2:in(i)-1,i)))')';
        
        %plot(  rssq(asd))
        
        %elong=in(i)*500*max(rssq(asd))/rssq(posall(:,end,1,4,i)-posall(:,end,in(i),4,i));
        %asd=rssq(asd);
        %asd = imresize(asd, [1 elong], 'nearest');
        %[juk,ind]=max(rssq(asd));
        [juk,ind]=max(rssq(asd(10:end,:)'));
%         rsamp=juk*10000/dist(i);
%         rsamp=round(rsamp*10);
%         
%         aaa(i)=rsamp;
%         asd=resample(asd',rsamp,10);
%         asd=asd';
%         [juk,ind]=max(rssq(asd));
       %    [juk,ind]=max((asd));
      %  ll=length(asd)-1;
        ll=length(asd(10:end,:)')-1;
       %plot( (-ind:ll-ind) ,rssq(asd)/max(rssq(asd)),'k')
        plot( rssq(asd(10:end,:)')/max(rssq(asd(10:end,:)')),'k')
       
        %plot( (-ind:ll-ind) ,rssq(asd(10:end,:)')/max(rssq(asd(10:end,:)')),'k')
       %  plot( (-ind:ll-ind) ,rssq(asd),'k')
      %        plot( (-ind:ll-ind) ,(asd)/max((asd)),'k')
%         
%         wtf=rssq(asd)/max(rssq(asd));
%         wtf(lala-1:lala+20)=NaN;
%          hop(:,i)=wtf(ind-20:ind+20);
       
        %/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
        hold on
        %plot(  rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i))))%/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
        hold on
    end
end

hop=hop(:,any(hop));

 av= nanmean(hop');
        figure;
        plot(av)
        vari= nanvar(hop');
        figure;
        plot(vari)




% 
% 
% 
% clear asd;
% for jj=2:22
%     
%     for ij=1:4
%         for i=3
%             if in(i)>10
%                 asd(:,:,ij)=cross((squeeze(posall(:,jj,1:in(i)-2,ij,i)-posall(:,jj,2:in(i)-1,ij,i)))',(squeeze(posall(:,jj,2:in(i)-1,ij,i)-posall(:,jj-1,2:in(i)-1,ij,i)))')';
%                 hold on
%                 
%                 % asd=cross((squeeze(testall(10:12,1:in(i)-2,i)-testall(10:12,2:in(i)-1,i)))',(squeeze(testall(10:12,2:in(i)-1,i)-testall(7:9,2:in(i)-1,i)))')';
%                 
%                 % plot(  rssq(asd))
%                 % [juk,ind]=max(rssq(asd));
%                 %  ll=length(asd)-1;
%                 % plot( (-ind:ll-ind) ,rssq(asd)/max(rssq(asd)))
%                 
%                 %plot(  rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i))))%/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
%                 
%             end
%         end
%     end
%     ll=length(asd);
%     asda=reshape(asd,[3,ll*4]);
%     plot( rssq(asda)/max(rssq(asda)))
%     %/max(rssq(squeeze(testall(4:6,2:end-1,i))-(testall(4:6,1:end-2,i)))))
%     hold on
%     
%     
% end
% 
