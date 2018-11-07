
%% nc14 find spots associated with each nuclei and plot their individual fluorescences over time to look at degree of correlation 
% count number of nuclei nc14 with spots 
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
nuclei=unique([CompiledParticles_14.Nucleus]);
nc14time=CompiledParticles_14(end).NucEnd-nc14;
SpotCount=[];
x=input('create all spot graphs?','s')
if x=='y'
    j=1;
else 
    j=0;
end
q=0;
Graphs=input('Make Spot Graphs for certain AP Position?','s')
if Graphs=='y'
    q=1;
end
if q==1
    APtoTrack=input('Which AP position to track? 0.xxx');
end


for ii=1:length(nuclei)
    APEstm(ii)=round(CompiledParticles_14(ii).MeanAP,3);
    for jj=1:length(APbinID)
        if APEstm(ii) < APbinID(jj)
            CompiledParticles_14(ii).APBin=APbinID(jj);
            break;
        end
    end
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));
    CompiledParticles_14(ii).FullFrames=zeros(1,nc14time);  %make field as long as number of frames in nc14
    if length(temp)==1
        SpotCount(ii,1)=1;
    end
    
    
    %need to only look at frames both exist in 
    
    if length(temp)==2
        SpotCount(ii,2)=1;
%        
if j==1
        
         figure(ii) 
         err1=CompiledParticles_14(temp(1)).FluoError*ones(size(CompiledParticles_14(temp(1)).Fluo));
         err2=CompiledParticles_14(temp(2)).FluoError*ones(size(CompiledParticles_14(temp(2)).Fluo));
         errorbar(CompiledParticles_14(temp(1)).Frame,CompiledParticles_14(temp(1)).Fluo,err1)
       %errorbar([CompiledParticles_14(temp(1)).FullFrames],[CompiledParticles_14(temp(1)).FullFluo],err1)
         hold on 
         errorbar(CompiledParticles_14(temp(2)).Frame,CompiledParticles_14(temp(2)).Fluo,err2)
         %errorbar([CompiledParticles_14(temp(2)).FullFrames],[CompiledParticles_14(temp(2)).FullFluo],err2)
         xlabel('Frame of nc14');
         ylabel('Spot Fluorescence');
         
         title(['2-Spot Fluorescence ''nucleus' num2str(nuclei(ii)) '' 'APbin' num2str(CompiledParticles_14(ii).APBin)]); %brackets needed to get string to work 
         
end
if q==1
    if [CompiledParticles_14(temp(1)).APBin]==APtoTrack
        figure(ii) 
         err1=CompiledParticles_14(temp(1)).FluoError*ones(size(CompiledParticles_14(temp(1)).Fluo));
         err2=CompiledParticles_14(temp(2)).FluoError*ones(size(CompiledParticles_14(temp(2)).Fluo));
         errorbar(CompiledParticles_14(temp(1)).Frame,CompiledParticles_14(temp(1)).Fluo,err1)
         hold on 
         errorbar(CompiledParticles_14(temp(2)).Frame,CompiledParticles_14(temp(2)).Fluo,err2)
       xlabel('Frame of nc14');
       ylabel('Spot Fluorescence')
       title(['2-Spot Fluorescence ''nucleus' num2str(nuclei(ii)) '' 'AP' '' num2str(APtoTrack)]); %brackets needed to get string to work 
    end
end

   
    
    end
    
        
end
OneSpot=sum(SpotCount(:,1));
    TwoSpot=sum(SpotCount(:,2));
    AllSpots=[OneSpot,TwoSpot];
    figure
    bar(AllSpots)


%% nc14 calculate the difference in fluorescence between spots in the same nucleus 

nc14time=length(ElapsedTime)-nc14;  
for ii=1:length(nuclei)
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));
%     CompiledParticles_14(temp(1)).nc14frames=CompiledParticles_14(temp(1)).Frame-nc14;   %to avoid overwriting frames transform them to their frame # within nc14
%     CompiledParticles_14(temp(1)).nc14frames=CompiledParticles_14(temp(1)).nc14frames +1;
%     for z=1:length([CompiledParticles_14(temp(1)).nc14frames])
%         CompiledParticles_14(temp(1)).Fixnc14frames(CompiledParticles_14(temp(1)).nc14frames(z))=CompiledParticles_14(temp(1)).nc14frames(z);
%     end
%     for q=1:length(CompiledParticles_14(temp(1)).Fluo)
%         for qq=1:length(CompiledParticles_14(temp(1)).Fixnc14frames)
%             if CompiledParticles_14(temp(1)).Fixnc14frames(qq) ~=0
%                 CompiledParticles_14(temp(1)).nc14Fluo(qq)=CompiledParticles_14(temp(1)).Fluo(q)
%             end
%         end
%     end
    %need to only look at nuclei with 2 spots  
    if length(temp)==2
%         CompiledParticles_14(temp(2)).nc14frames=CompiledParticles_14(temp(2)).Frame -nc14;
%         CompiledParticles_14(temp(2)).nc14frames=CompiledParticles_14(temp(2)).nc14frames +1;
%         for z=1:length([CompiledParticles_14(temp(2)).nc14frames])
%         CompiledParticles_14(temp(2)).Fixnc14frames(CompiledParticles_14(temp(2)).nc14frames(z))=CompiledParticles_14(temp(2)).nc14frames(z);
%         end
%         for q=1:length(CompiledParticles_14(temp(2)).Fluo)
%         for qq=1:length(CompiledParticles_14(temp(2)).Fixnc14frames)
%             if CompiledParticles_14(temp(2)).Fixnc14frames(qq) ~=0
%                 CompiledParticles_14(temp(2)).nc14Fluo(qq)=CompiledParticles_14(temp(2)).Fluo(q)
%             end
%         end
%         end
%     CompiledParticles_14(temp(1)).nc14Fluo(CompiledParticles_14(temp(1)).nc14Fluo==0)=nan;
%     CompiledParticles_14(temp(2)).nc14Fluo(CompiledParticles_14(temp(2)).nc14Fluo==0)=nan;
%         tempN=[];
%         tempN2=[];
%         for jj=1:length(schnitzcells(nuclei(ii)).frames)%want to compare frames in schnitz to those of each spot in CompiledParticles  
%        interestframe=schnitzcells(nuclei(ii)).frames(jj);
%        interestframenc14=interestframe-nc14;           %want frames within nc14 for corresponding to particle frames of nc14
%        interestframenc14=interestframe +1;
%        if interestframenc14 >= 1
%            tempN=find([CompiledParticles_14(temp(1)).Fixnc14frames]==interestframenc14);
%            %tempN=find([CompiledParticles_14(temp(1)).Frame]==interestframe);  %make frame where nucleus exist but no spot equal to 1 so can find later on 
%        if isempty(tempN)
%            CompiledParticles_14(temp(1)).Fixnc14frames(interestframenc14)=1;
%            CompiledParticles_14(temp(1)).nc14Fluo(interestframenc14)=0;
%            %CompiledParticles_14(temp(1)).Fluo(CompiledParticles_14(temp(1)).Fluo==0)=nan;
%            %CompiledParticles_14(temp(1)).Fluo(interestframenc14)=0;
%            CompiledParticles_14(temp(1)).Fixnc14frames(CompiledParticles_14(temp(1)).Fixnc14frames==0)=nan;
%            %CompiledParticles_14(temp(1)).Frame(interestframenc14)=0;
%        
%        end
       
%        tempN2=find([CompiledParticles_14(temp(2)).Fixnc14frames]==interestframe);
%        if isempty(tempN2)
%            CompiledParticles_14(temp(2)).Fixnc14frames(interestframenc14)=1;
%            CompiledParticles_14(temp(2)).nc14Fluo(interestframenc14)=0;    %need to make fluo field as long as frame to match up correctly 
%            %CompiledParticles_14(temp(2)).Fluo(CompiledParticles_14(temp(2)).Fluo==0)=nan;    %make fields that created but not corresponding to frame looking at equal to nan instead of 0 so don't have all points as 0 fluorescence
%            %CompiledParticles_14(temp(2)).Fluo(interestframenc14)=0; 
%            CompiledParticles_14(temp(2)).Fixnc14frames(CompiledParticles_14(temp(2)).Frame==0)=nan;
%            %CompiledParticles_14(temp(2)).Frame(interestframenc14)=0;
% 
%        end
%        end
      
        DiffArray=zeros(2,nc14time);
        for jj=1:length(CompiledParticles_14(temp(1)).Frame)
%             
                 DiffArray(1,jj)=CompiledParticles_14(temp(1)).Fluo(jj); 
                 

            
            SpotDiff(ii).MeanAP=CompiledParticles_14(temp(1)).MeanAP;  %do I want mean here?
        SpotDiff(ii).Err1=CompiledParticles_14(temp(1)).FluoError*ones(size(CompiledParticles_14(temp(1)).Fluo)); %calculate error for each spot for later reference/plotting
        SpotDiff(ii).Frame1=CompiledParticles_14(temp(1)).Frame;
        %SpotDiff(ii).nc14Frame1=CompiledParticles_14(temp(1)).Fixnc14frames;
        SpotDiff(ii).APBin=CompiledParticles_14(temp(1)).APBin;
        end
        for qq=1:length(CompiledParticles_14(temp(2)).Frame)
%             
            DiffArray(2,qq)=CompiledParticles_14(temp(2)).Fluo(qq);
            %end
        SpotDiff(ii).Err2=CompiledParticles_14(temp(2)).FluoError*ones(size(CompiledParticles_14(temp(2)).Fluo));
        SpotDiff(ii).Frame2=CompiledParticles_14(temp(2)).Frame;
        %SpotDiff(ii).nc14Frame2=CompiledParticles_14(temp(2)).Fixnc14frames;
        end
%         for jj=1:length(schnitzcells(nuclei(ii)).frames)%want to compare frames in schnitz to those of each spot in CompiledParticles  
%         interestframe=schnitzcells(nuclei(ii)).frames(jj);
%         interestframenc14=interestframe-nc14time;           %want frames within nc14 for corresponding to particle frames of nc14
%         interestframenc14=interestframe +1;
%             tempN=find([CompiledParticles_14(temp(1)).Frame]==interestframe);
%         if isempty(tempN)
%             CompiledParticles_14(temp(1)).Frame(interestframenc14)=0;
%             CompiledParticles_14(temp(1)).Frame(CompiledParticles_14(temp(1)).Frame==0)=nan;
%             CompiledParticles_14(temp(1)).Frame(interestframenc14)=0;
%  
%         end
%         tempN2=find([CompiledParticles_14(temp(2)).Frame]==interestframe);
%         if isempty(tempN2);
%             CompiledParticles_14(temp(2)).Frame(interestframenc14)=0;
%             CompiledParticles_14(temp(2)).Frame(CompiledParticles_14(temp(2)).Frame==0)=nan;
%             CompiledParticles_14(temp(2)).Frame(interestframenc14)=0;
%  
%         end
%         end
       
       
           SpotDiff(ii).SpotOne=[DiffArray(1,:)];
         SpotDiff(ii).SpotTwo=[DiffArray(2,:)];
         SpotDiff(ii).AbsDiffBtwn=[abs(DiffArray(1,:)-DiffArray(2,:))];
         
         SpotDiff(ii).SpotCorr=[corrcoef(DiffArray')];
% DiffArrayAvg=mean(DiffArray,2);

    end
end

%% 
z=0;
DoIt=input('Create spot graphs for certain AP positions?', 's')
if DoIt=='y'
    z=1;
else 
    z=0;
end
if z==1
    Which=input('Which AP bin? 0.xx')
    temp=[find([SpotDiff.APBin]==Which)];
    for tt=1:length(temp)
      figure(tt) 
      err1=SpotDiff(temp(tt)).Err1;
         %err1=SpotDiff(temp(tt)).Err1*ones(size(SpotDiff(temp(tt)).SpotOne));
         err2=SpotDiff(temp(tt)).Err2;
         %err2=SpotDiff(temp(tt)).Err1*ones(size(SpotDiff(temp(tt)).SpotTwo));
         errorbar([1:length(SpotDiff(temp(tt)).SpotOne)],[SpotDiff(temp(tt)).SpotOne],err1)
       %errorbar([CompiledParticles_14(temp(1)).FullFrames],[CompiledParticles_14(temp(1)).FullFluo],err1)
         hold on 
         errorbar(1:length(SpotDiff(temp(tt)).SpotTwo),SpotDiff(temp(tt)).SpotTwo,err2)
         %errorbar([CompiledParticles_14(temp(2)).FullFrames],[CompiledParticles_14(temp(2)).FullFluo],err2)
         xlabel('Frame of nc14');
         ylabel('Spot Fluorescence');
         
         title(['2-Spot Fluorescence ' '']) %'APbin' num2str(SpotDiff(ii).APBin)]);
    end
end

         %round AP positions of nuclei to closest AP bin
APEstm=[];

for j=1:length(SpotDiff)
    if ~isempty(SpotDiff(j).MeanAP)
        
    APEstm(j)=round(SpotDiff(j).MeanAP,3);
    
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            SpotDiff(j).APBin=APbinID(jj);
            break;
        end
    end
    end
    
end

 for ii=1:length(SpotDiff)
      if isempty(SpotDiff(ii).APBin)
          SpotDiff(ii).APBin=0;
      end
  end

Bins=unique([SpotDiff.APBin]);%want to find mean correlation of spots at each APbin for fewer points/easier to see graph when combining
%APCorr=zeros(22,length(Bins));
%APBins=Data(1).APbinID;
        APCorr=zeros(1,length(APbinID));
for zz=1:length(Bins)
    temp3=find([APbinID==Bins(zz)]);
    temp2=find([SpotDiff.APBin]==Bins(zz));%issue with find I think has to do with empty values 
    for qq=1:length(temp2)
        if ~isempty(SpotDiff(temp2(qq)).SpotCorr)
    APCorr(qq,temp3)=[SpotDiff(temp2(qq)).SpotCorr(1,2)];
        end
    end 
end
APCorr(APCorr==0)=nan;
%APCorr=APCorr(:,2:end); %remove first column of 0's corresponding to bin of "0"
MeanAPCorr=nanmean(APCorr,1);

figure
h=plot(APbinID,MeanAPCorr,'LineWidth',2)
% for ii=1:length(MeanAPCorr)
% plot(Bins(2:end),MeanAPCorr(ii))
% hold on 
% end


hold on 
% figure
for qq=1:length(SpotDiff)
    if ~isempty(SpotDiff(qq).SpotCorr)
    plot(SpotDiff(qq).MeanAP,SpotDiff(qq).SpotCorr(2),'o')
    hold on 
    end
    
end
xlabel('Mean AP Position')
ylabel('Correlation of spot fluorescence')
ylim([-0.5 1])
xlim([0.3 0.75])
title('Correlation of loci transcription')
legend([h], 'Mean Correlation at each AP bin','Location','best');

save('SpotCorrelation','SpotDiff')
save('MeanAPCorrelation','MeanAPCorr')



% scatter(DiffArrayAvg(1,:),DiffArrayAvg(2,:))

