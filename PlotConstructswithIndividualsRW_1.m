 %% plot construct mean fluorescence on top of individual traces of construct 
 %nc 14 AP bin 20
 EmbryoName=[];
 for qq=1:length(Data)
     nc14time=length(Data(qq).ElapsedTime)-Data(qq).nc14;
     nc14all(:,qq)=nc14time;
     %EmbryoName(:,qq)=Data(qq).SetName(9:end-1); %also for intention of
     %plotting
 end
 maxnc14=max(nc14all);
 
 nc14AP20=zeros(maxnc14,length(Data));
 nc14error=zeros(maxnc14,length(Data));
for jj=1:length(Data)   %collect info and make mean of just nc14
nc14=Data(jj).nc14;

    nc14AP20(1:nc14all(jj)+1,jj)=[Data(jj).MeanVectorAP(nc14:end,20)]; %each column is a different embryo where each row is the mean fluorescence at APbin 20 at that time point in nc14
nc14error(1:nc14all(jj)+1,jj)=Data(jj).SDVectorAP(nc14:end,20);
end
ConstructMean=[nanmean(nc14AP20,2)]; %mean of each row aka each time point across all embryos
ConstructError=[nanmean(nc14error,2)];
constructtype=input('Which construct?','s');

 figure
 h=[];
 legendname=[];
 for ii=1:length(Data)
     nc14=Data(ii).nc14;
     h(ii)=plot(Data(ii).MeanVectorAP(nc14:end,20))%h(ii)=errorbar(Data(ii).MeanVectorAP(nc14:end,20),Data(ii).SDVectorAP(nc14:end,20)) %included the h part to go back and figure out how to add in legend 
     e=plot(Data(ii).MeanVectorAP(nc14:end,20));
     e.Marker='.';
     %e.LineWidth=1;
     e.MarkerSize=10;
     %e.CapSize=5;
     %legendname(ii)=Data(ii).SetName(9:18);
     
     hold on 
    
     
 end
 
 %errorbar(ConstructMean,ConstructError);
 h1=plot(ConstructMean,'LineWidth',3,'Color','k')
 %e=errorbar(ConstructMean,ConstructError);
%      e.Marker='.';
%      e.MarkerSize=10;
%      e.CapSize=5;
%      e.Color='k';
 legend([h1],{'Construct mean'});
 xlabel('Frame'); ylabel('Mean Fluorescence');
 ylim([0 45000]);
 
 title(['Mean Fluorescence nc14',' ',constructtype,' ', 'n=',num2str(length(Data))]); 
 

 
     
 
    %for jj=length(SpotNuc)
        

% ConstructList={'KrBoth'; 'KrDist';'KrProx'}
% MinParticles=5;
% MinEmbryos=1;
% 
% for i=length(ConstructList)
%     construct=ConstructList{i}
%     Data=LoadMS2SetsCS(construct);
%     [TotalProd,TotalProdError,TotalProdN,MeanTotalProd,SDTotalProd,SETotalProd] =IntegratemRNA(Data,MinParticles,MinEmbryos);
%     AllConstructs(i).MeanTotalProd=MeanTotalProd;
%     AllConstructs(i).SETotalProd=SETotalProd;
%     AllConstructs(i).SDTotalProd=SDTotalProd; 
% end
% 
%   
%     
%     %% Total mRNA over space
%     %NC 14
%     
%     
%     
%     