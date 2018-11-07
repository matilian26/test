%% nc14 find spots associated with each nuclei and plot their individual fluorescences over time to look at degree of correlation 

nc_number=[CompiledParticles.nc];
CompiledParticles_13=CompiledParticles(nc_number==13);
nuclei=unique([CompiledParticles_13.Nucleus]);
nc13time=nc14-nc13;%(max([CompiledParticles_13.NucEnd]))-nc13;
SpotCount=[];
RunAll=input('ask about graphs?','s');
if RunAll=='y'
x=input('create all spot graphs?','s')
else
    x='n'
end
if x=='y'
    j=1;
else 
    j=0;
end
q=0;
if RunAll=='y'
Graphs=input('Make Spot Graphs for certain AP Position?','s')
else
    Graphs='n'
end
if Graphs=='y'
    q=1;
end
if q==1
    APtoTrack=input('Which AP position to track? 0.xxx');
end


for tt=1:length(CompiledParticles_13)
    APEstm(tt)=round(CompiledParticles_13(tt).MeanAP,3);
    for jj=1:length(APbinID)
        if APEstm(tt) < APbinID(jj)
            CompiledParticles_13(tt).APBin=APbinID(jj);
            break;
        end
    end
end 
    for ii=1:length(nuclei)
    temp=find([CompiledParticles_13.Nucleus]==nuclei(ii));
    CompiledParticles_13(ii).FullFrames=zeros(1,nc13time);  %make field as long as number of frames in nc14
    if length(temp)==1
        SpotCount(ii,1)=1;
    end
    
    
    %need to only look at frames both exist in 
    
    if length(temp)==2
        SpotCount(ii,2)=1;
%        
if j==1
        
         figure(ii) 
         err1=CompiledParticles_13(temp(1)).FluoError*ones(size(CompiledParticles_13(temp(1)).Fluo));
         err2=CompiledParticles_13(temp(2)).FluoError*ones(size(CompiledParticles_13(temp(2)).Fluo));
         errorbar(CompiledParticles_13(temp(1)).Frame,CompiledParticles_13(temp(1)).Fluo,err1)
       %errorbar([CompiledParticles_14(temp(1)).FullFrames],[CompiledParticles_14(temp(1)).FullFluo],err1)
         hold on 
         errorbar(CompiledParticles_13(temp(2)).Frame,CompiledParticles_13(temp(2)).Fluo,err2)
         %errorbar([CompiledParticles_14(temp(2)).FullFrames],[CompiledParticles_14(temp(2)).FullFluo],err2)
         xlabel('Frame of nc14');
         ylabel('Spot Fluorescence');
         
         title(['2-Spot Fluorescence ''nucleus' num2str(nuclei(ii)) '' 'APbin' num2str(CompiledParticles_13(ii).APBin)]); %brackets needed to get string to work 
         
end
if q==1
    if [CompiledParticles_13(temp(1)).APBin]==APtoTrack
        figure(ii) 
         err1=CompiledParticles_13(temp(1)).FluoError*ones(size(CompiledParticles_13(temp(1)).Fluo));
         err2=CompiledParticles_13(temp(2)).FluoError*ones(size(CompiledParticles_13(temp(2)).Fluo));
         errorbar(CompiledParticles_13(temp(1)).Frame,CompiledParticles_13(temp(1)).Fluo,err1)
         hold on 
         errorbar(CompiledParticles_13(temp(2)).Frame,CompiledParticles_13(temp(2)).Fluo,err2)
       xlabel('Frame of nc14');
       ylabel('Spot Fluorescence')
       title(['2-Spot Fluorescence ''nucleus' num2str(nuclei(ii)) '' 'AP' '' num2str(APtoTrack)]); %brackets needed to get string to work 
    end
end

   
    
    end
    
        
    end
    
    if RunAll=='y'
BarAll=input('bar graph of all spots?','s');
    else
        BarAll='n'
    end
if BarAll=='y'
OneSpot=sum(SpotCount(:,1));
    TwoSpot=sum(SpotCount(:,2));
    AllSpots=[OneSpot,TwoSpot];
    figure
    bar(AllSpots)
end


%% nc14 calculate the difference in fluorescence between spots in the same nucleus 

 LastFrame=max([CompiledParticles_13.NucEnd]);
 nc13frames=[nc13:nc14];
for ii=1:length(nuclei)      % I have a feeling this is an issue
    temp=find([CompiledParticles_13.Nucleus]==nuclei(ii));  
    CompiledParticles_13(temp(1)).nc13Frame=zeros(1,length(ElapsedTime));       %make frame value sit at that index position in the .nc14Frames structure
    for qq=1:length(CompiledParticles_13(temp(1)).Frame)
        CompiledParticles_13(temp(1)).nc13Frame(CompiledParticles_13(temp(1)).Frame(qq))=CompiledParticles_13(temp(1)).Frame(qq);
    end
    counter=0;      %Also make nc14Fluo field same length as nc14Frames so can have 0 values for frames where schintz exists but not fluor 
    CompiledParticles_13(temp(1)).nc13Fluo=nan(1,length(ElapsedTime));
    
    for jj=1:length(CompiledParticles_13(temp(1)).nc13Frame)
        if any(CompiledParticles_13(temp(1)).nc13Frame(jj))%(~isempty(CompiledParticles_14(temp(1)).nc14Frame(jj))) & (CompiledParticles_14(temp(1)).nc14Frame(jj)~=0)
            counter=counter+1;
            CompiledParticles_13(temp(1)).nc13Fluo(jj)=CompiledParticles_13(temp(1)).Fluo(counter);
        end
    end
       CompiledParticles_13(temp(1)).nc13Fluo(CompiledParticles_13(temp(1)).nc13Fluo==0)=nan;

    
    if length(temp)==2
        CompiledParticles_13(temp(2)).nc13Frame=zeros(1,length(ElapsedTime));       %make frame value sit at that index position in the .nc14Frames structure
    for zz=1:length(CompiledParticles_13(temp(2)).Frame)
        CompiledParticles_13(temp(2)).nc13Frame(CompiledParticles_13(temp(2)).Frame(zz))=CompiledParticles_13(temp(2)).Frame(zz);
    end
        counter2=0;%[1:length(CompiledParticles_14(temp(2)).Fluo)];
        
        CompiledParticles_13(temp(2)).nc13Fluo=nan(1,length(ElapsedTime));
             
    for jj=1:length(CompiledParticles_13(temp(2)).nc13Frame)
        
        if any(CompiledParticles_13(temp(2)).nc13Frame(jj)) %(~isempty(CompiledParticles_14(temp(2)).nc14Frame(jj))) & (CompiledParticles_14(temp(2)).nc14Frame(jj)~=0)
            counter2=counter2+1;
            CompiledParticles_13(temp(2)).nc13Fluo(jj)=CompiledParticles_13(temp(2)).Fluo(counter2);
        end
    end
    
    CompiledParticles_13(temp(2)).nc13Fluo(CompiledParticles_13(temp(2)).nc13Fluo==0)=nan;
        
            
        %end
        
        
%        
    end
%    
    %need to only look at nuclei with 2 spots  
    
        for jj=1:length(nc13frames)
            tempa=[];tempb=[];tempa2=[];tempb2=[];
        tempa=find([schnitzcells(CompiledParticles_13(temp(1)).Nucleus).frames]==nc13frames(jj));
        tempb=find([CompiledParticles_13(temp(1)).nc13Frame]==nc13frames(jj));

        %schnitz exists but no spot at that frame
        if (~isempty(tempa)) & (isempty(tempb))%(tempa ~=0) & (isempty(tempb) | tempb==0) %~isempty(tempa) & isempty(tempb)
            CompiledParticles_13(temp(1)).nc13Frame(nc13frames(jj))=1;
        end
        end
        for jj=1:length(CompiledParticles_13(temp(1)).nc13Frame)   %For nuclei where found exist in schnitz but no particle, have fluo value of 0 (other indices where no frame in compiledpar and not id-ed in schnitz should be nan's rn
            if CompiledParticles_13(temp(1)).nc13Frame(jj)==1
                CompiledParticles_13(temp(1)).nc13Fluo(jj)=0;
            end
        end
        if length(temp)==2
            for jj=1:length(nc13frames)
        tempa2=find([schnitzcells(CompiledParticles_13(temp(2)).Nucleus).frames]==nc13frames(jj));
       
        tempb2=find([CompiledParticles_13(temp(2)).nc13Frame]==nc13frames(jj));
        if (~isempty(tempa2)) & (isempty(tempb2))%(tempa2 ~=0) & (isempty(tempb2) | tempb2==0)%any(tempa2) & ~any(tempb2)
            CompiledParticles_13(temp(2)).nc13Frame(nc13frames(jj))=1;
        end 
        end
       
%         for jj=1:length(CompiledParticles_14(temp(1)).nc14Frame)   %For nuclei where found exist in schnitz but no particle, have fluo value of 0 (other indices where no frame in compiledpar and not id-ed in schnitz should be nan's rn
%             if CompiledParticles_14(temp(1)).nc14Frame(jj)==1
%                 CompiledParticles_14(temp(1)).nc14Fluo(jj)=0;
%             end
%         end
        for jj=1:length(CompiledParticles_13(temp(2)).nc13Frame)
            if CompiledParticles_13(temp(2)).nc13Frame(jj)==1
                CompiledParticles_13(temp(2)).nc13Fluo(jj)=0;
            end
        end
        end 
        
       
        
        DiffArray=[];%zeros(2,length(ElapsedTime));
        for jj=1:length(CompiledParticles_13(temp(1)).nc13Frame)
%              
                 DiffArray(1,jj)=CompiledParticles_13(temp(1)).nc13Fluo(jj); 
                
            SpotDiff(ii).MeanAP=CompiledParticles_13(temp(1)).MeanAP;  %do I want mean here?
        SpotDiff(ii).OriginalParticle=CompiledParticles_13(temp(1)).OriginalParticle;
            SpotDiff(ii).Err1=CompiledParticles_13(temp(1)).FluoError*ones(size(CompiledParticles_13(temp(1)).nc13Frame)); %calculate error for each spot for later reference/plotting
        SpotDiff(ii).Frame1=[CompiledParticles_13(temp(1)).Frame];
        SpotDiff(ii).nc14Frame1=[CompiledParticles_13(temp(1)).nc13Frame];
        SpotDiff(ii).APBin=CompiledParticles_13(temp(1)).APBin;
        SpotDiff(ii).Nucleus=CompiledParticles_13(temp(1)).Nucleus;
        SpotDiff(ii).SpotOne=[DiffArray(1,:)];
        SpotDiff(ii).TotalmRNAOne=[CompiledParticles_13(temp(1)).TotalmRNA];
        end
        if length(temp)==2
        for qq=1:length(CompiledParticles_13(temp(2)).nc13Fluo)
%             
            DiffArray(2,qq)=CompiledParticles_13(temp(2)).nc13Fluo(qq);
            end
        
        SpotDiff(ii).Original2ndParticle=CompiledParticles_13(temp(2)).OriginalParticle;

        SpotDiff(ii).Err2=CompiledParticles_13(temp(2)).FluoError*ones(size(CompiledParticles_13(temp(2)).nc13Frame));
        SpotDiff(ii).Frame2=[CompiledParticles_13(temp(2)).Frame];
        SpotDiff(ii).nc14Frame2=[CompiledParticles_13(temp(2)).nc13Frame];
        SpotDiff(ii).TotalmRNATwo=[CompiledParticles_13(temp(2)).TotalmRNA];
        
%        
       
       
           
         SpotDiff(ii).SpotTwo=[DiffArray(2,:)];
         SpotDiff(ii).AbsDiffBtwn=[abs(DiffArray(1,:)-DiffArray(2,:))];
        
         SpotDiff(ii).SpotCorr=[corrcoef(DiffArray','Rows','complete')];
% DiffArrayAvg=mean(DiffArray,2);

    end
end
%Get rid of all empty field values 
for ss=1:length(SpotDiff)
    if isempty(SpotDiff(ss).SpotOne)
        SpotDiff(ss).SpotOne=nan;
    end
    if isfield(SpotDiff, 'SpotTwo')
    if isempty(SpotDiff(ss).SpotTwo)
        SpotDiff(ss).SpotTwo=nan;
    end
    if isempty(SpotDiff(ss).Original2ndParticle)
        SpotDiff(ss).Original2ndParticle=nan;
    end
    if isempty(SpotDiff(ss).Err2)
        SpotDiff(ss).Err2=nan;
    end
    if isempty(SpotDiff(ss).Frame2)
        SpotDiff(ss).Frame2=nan;
    end
    if isempty(SpotDiff(ss).nc13Frame2)
        SpotDiff(ss).nc14Frame2=nan;
    end
    if isempty(SpotDiff(ss).AbsDiffBtwn)
        SpotDiff(ss).AbsDiffBtwn=nan;
    end
    if isempty(SpotDiff(ss).SpotCorr)
        SpotDiff(ss).SpotCorr=nan;
    end
    end
end

%% 
z=0;
if RunAll=='y'
DoIt=input('Create spot graphs for certain AP positions?', 's')
else
    DoIt='n';
end
if DoIt=='y'
    z=1;
else 
    z=0;
end
if z==1
    Which=input('Which AP bin? 0.xxxx')
    temp=[find([SpotDiff.APBin]==Which)];
    for tt=1:length(temp)
        
      figure(tt) 
      err1=SpotDiff(temp(tt)).Err1;
     x=[1:length(SpotDiff(temp(tt)).SpotOne)];
         err2=SpotDiff(temp(tt)).Err2;
         y=[SpotDiff(temp(tt)).SpotOne];
         idx=~isnan(SpotDiff(temp(tt)).SpotOne);
         if ~isempty(x)
         y=interp1(x(idx),y(idx),x,'linear');
         end
         plot(1:length(SpotDiff(temp(tt)).SpotOne),y,'-*')
         %errorbar(1:length(SpotDiff(temp(tt)).SpotOne),SpotDiff(temp(tt)).SpotOne,err1)
         hold on 
         if ~isempty(SpotDiff(temp(tt)).SpotTwo);
         x2=[1:length(SpotDiff(temp(tt)).SpotTwo)];
         y2=[SpotDiff(temp(tt)).SpotTwo];
         idx2=~isnan(SpotDiff(temp(tt)).SpotTwo);
         y2=interp1(x2(idx2),y2(idx2),x2,'linear');
         plot(1:length(SpotDiff(temp(tt)).SpotTwo),y2,'-*')
         end
     %errorbar(x2(idxs2),SpotDiff(temp(tt)).SpotTwo(idxs2),err2)    
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
if isfield(SpotDiff, 'SpotTwo')
for zz=1:length(Bins)
    temp3=find([APbinID==Bins(zz)]);
    temp2=find([SpotDiff.APBin]==Bins(zz));%issue with find I think has to do with empty values 
    for qq=1:length(temp2)
        if ~isempty(SpotDiff(temp2(qq)).SpotCorr)
            if length(SpotDiff(temp2(qq)).SpotCorr)==1
                APCorr(qq,temp3)=[SpotDiff(temp2(qq)).SpotCorr(1)];
            else
    APCorr(qq,temp3)=[SpotDiff(temp2(qq)).SpotCorr(1,2)];
            end
        end
    end 
end
APCorr(APCorr==0)=nan;
%APCorr=APCorr(:,2:end); %remove first column of 0's corresponding to bin of "0"
MeanAPCorr=nanmean(APCorr,1);

%figure
if RunAll=='y'
h=plot(APbinID,MeanAPCorr,'LineWidth',2)
% for ii=1:length(MeanAPCorr)
% plot(Bins(2:end),MeanAPCorr(ii))
% hold on 
% end


hold on 
% figure

for qq=1:length(SpotDiff)
    if ~isempty(SpotDiff(qq).SpotCorr)
        if length(SpotDiff(qq).SpotCorr)==1
            plot(SpotDiff(qq).MeanAP,SpotDiff(qq).SpotCorr(1),'o')
        else
    plot(SpotDiff(qq).MeanAP,SpotDiff(qq).SpotCorr(2),'o')
        end
    hold on 
    end
    
end
xlabel('Mean AP Position')
ylabel('Correlation of spot fluorescence')
ylim([-0.5 1])
xlim([0.2 0.75])
title('Correlation of loci transcription')
legend([h], 'Mean Correlation at each AP bin','Location','best');
% saveas(gcf,'MeanCorrelationAPAdj.png')
end
end
% 
 %save('SpotCorrelationAdj','SpotDiff')
% save('MeanAPCorrelationAdj','MeanAPCorr')



% scatter(DiffArrayAvg(1,:),DiffArrayAvg(2,:))

