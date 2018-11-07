%% nc14
% count number of nuclei nc14 with spots 
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
nuclei=unique([CompiledParticles_14.Nucleus]);
NucleiLife=[];

for ii=1:length(nuclei)   %Find out how long each nucleus exists for 
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));
    temp2=length(schnitzcells(nuclei(ii)).frames);%temp2=length(schnitzcells(temp(1)).frames);
    NucleiLife(ii)=temp2;
%     
end
nucleithresh=[];   
%Find nuclei who exist for >=20 frames 
for jj=1:length(nuclei)
    if NucleiLife(jj) >=20
        nucleithresh(jj)=nuclei(jj);
    end
end
nucleithresh=unique([nucleithresh]); %Want to get rid of 0's
nucleithresh=[nucleithresh(2:end)]; %remove last "unique" 0 value

for ii=1:length(nucleithresh)    %only want to look at nuclei that exist >= 20 frames
    temp=find([CompiledParticles_14.Nucleus]==nucleithresh(ii));
    num_spots=length(temp);
    num_spots_ap(ii).spots=num_spots;
    num_spots_ap(ii).MeanAP=CompiledParticles_14(temp(1)).MeanAP;

% for ii=1:length(nuclei)
%     temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));  %could this miscount so that a nucleus that has 2 spots, then one for ex would be counted as having 3 spots?
%     num_spots=length(temp);
%     num_spots_ap(ii).spots=num_spots;
%     num_spots_ap(ii).MeanAP=CompiledParticles_14(temp(1)).MeanAP;
%     
end
%round AP positions of nuclei to closest AP bin
APEstm=[num_spots_ap.MeanAP];
for j=1:length(nucleithresh)
    APEstm(j)=round(num_spots_ap(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            num_spots_ap(j).APBin=APbinID(jj);
            break;
        end
    end
end

% Deal with nuclei with no spot 
%first translate x,y coordinates to AP position 
%Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
NucPosition=[];



for i=1:length(schnitzcells)     %find the AP position of each nucleus across time 
        for j=1:length(schnitzcells(i).frames)
            
            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((schnitzcells(i).ceny(j)-coordAZoom(2))./(schnitzcells(i).cenx(j)-coordAZoom(1)));
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-schnitzcells(i).ceny(j)).^2+(coordAZoom(1)-schnitzcells(i).cenx(j)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            NucTimeFrame=schnitzcells(i).frames(j);  %making the columns the time frame, each row is a nucleus 
            NucPosition(i,NucTimeFrame)=APPositions/APLength;
            
        end
end

 NoSpotNucLife=[];
  NoSpotPosition=[];
  for ii=1:length(NucPosition)    %only keep nuclei that exist >=20 frames
        temp=[length(schnitzcells(ii).frames)];
        NoSpotNucLife(ii)=temp;
    end
    for jj=1:length(NucPosition)
        for qq=1:size(NucPosition,2)
        if NoSpotNucLife(jj) >=20
            NoSpotPosition(jj,qq)=NucPosition(jj,qq);
        end
        end
    end
     

%remove Nuclei that identified as having spots in above 

    SpotNuc=unique([CompiledParticles_14.Nucleus]); %find the nuclei that have spots. Only want to remove those that have spots in nc14 (??) 
    for jj=1:length(SpotNuc)
        NoSpotPosition(SpotNuc(jj),:)=nan;      %make sure these aren't used later in plots
        
    end
   
  
    %Isolate nuclei in nc14
   NoSpotPosition_nc14=NoSpotPosition(:,nc14:end); 
    %change 0's indicating nuclei not present at that time point to nan's
    %so don't have all off nuclei mapping to first AP bin when take average
    %of AP position
    %for ii=1:length(NoSpotPosition_nc14)
        NoSpotPosition_nc14(NoSpotPosition_nc14==0)=nan; %for ii=1:length(NoSpotPosition)
        
        
   
    
            
    %incorporate no-spot nuclei into array of nuclei with spots and change
    %AP position to their mean positions as did for nuclei w spots
    for j=1:length(NoSpotPosition_nc14)
        MeanAPNoSpot=nanmean(NoSpotPosition_nc14(j,:));
    num_spots_ap(length(num_spots_ap)+1).MeanAP=MeanAPNoSpot;%[mean(NoSpotPosition(j,:),2,'omitnan')]
    end
    %fill in 0 values for no spot nuclei
    for jj=1:length(num_spots_ap)
         if isempty(num_spots_ap(jj).spots)
         num_spots_ap(jj).spots =0;
        end
    end
    
    %round AP positions of nuclei to closest AP bin
APEstm=[num_spots_ap.MeanAP];
for j=1:length(num_spots_ap)
    APEstm(j)=round(num_spots_ap(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            num_spots_ap(j).APBin=APbinID(jj);
            break;
        end
    end
end
    

    %put number of spots in a nuclei as column and AP position of
    %associated nuclei as row to make stacked bar graph 
    APbins = unique([num_spots_ap.APBin]);   %need to add part to add up all the zero spot nuclei, also make sure only the double counted nuclei were changed to nan??
SpotsAP=nan(length(APbinID),3);%SpotsAP=zeros(length(APbins),3);
DispNuc=input('Plot fraction or number of nuclei? f/n','s')

for ii=1:length(APbins)
    jj = find([num_spots_ap.APBin]==APbins(ii));
    APcolumn=find([APbinID]==APbins(ii));        %want row number to correspond to column number of AP bin in array of all AP bins across embryo
    
    Spotnumber=[];
   
    Spotnumber=[num_spots_ap(jj).spots];
    OneSpot=sum([Spotnumber]==1);
    TwoSpot=sum([Spotnumber]==2);
    NoSpot=sum([Spotnumber]==0);
    if DispNuc=='f'
    TotalNucs=OneSpot+TwoSpot+NoSpot;
    SpotsAP(APcolumn,1)=NoSpot/TotalNucs;
    SpotsAP(APcolumn,2)=OneSpot/TotalNucs;
    SpotsAP(APcolumn,3)=TwoSpot/TotalNucs;
    else
        SpotsAP(APcolumn,1)=NoSpot;
        SpotsAP(APcolumn,2)=OneSpot;
        SpotsAP(APcolumn,3)=TwoSpot;
    end
    
end
figure 
bar(APbinID,SpotsAP,'stacked')
xlabel('Egg Length'); 
if DispNuc=='f'
ylabel('fraction of nuclei')
ylim([0 1.5]);
else
    ylabel('# of nuclei')
    ylim([0 60]);
end
legend('no spot','one spot','two spots','Location','best');
title('ON nuclei nc14');
 
xlim([APbins(1)-0.25 APbins(end)+0.25]);
if DispNuc=='f'
save('SpotsAPRatioFract','SpotsAP');
saveas(gcf,'SpotOnRatioAPFract.jpg');
else
    save('SpotsAPRatio','SpotsAP');
    saveas(gcf,'SpotONRatioAP.jpg');
end
    
%% Limit of 10 frame existence 
% nc14
% count number of nuclei nc14 with spots 
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
nuclei=unique([CompiledParticles_14.Nucleus]);
NucleiLife=[];

for ii=1:length(nuclei)   %Find out how long each nucleus exists for 
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));
    temp2=length(schnitzcells(nuclei(ii)).frames);%temp2=length(schnitzcells(temp(1)).frames);
    NucleiLife(ii)=temp2;
%     
end
nucleithresh=[];   
%Find nuclei who exist for >=20 frames 
for jj=1:length(nuclei)
    if NucleiLife(jj) >=10
        nucleithresh(jj)=nuclei(jj);
    end
end
nucleithresh=unique([nucleithresh]); %Want to get rid of 0's
nucleithresh=[nucleithresh(2:end)]; %remove last "unique" 0 value

for ii=1:length(nucleithresh)    %only want to look at nuclei that exist >= 20 frames
    temp=find([CompiledParticles_14.Nucleus]==nucleithresh(ii));
    num_spots=length(temp);
    num_spots_ap(ii).spots=num_spots;
    num_spots_ap(ii).MeanAP=CompiledParticles_14(temp(1)).MeanAP;

% for ii=1:length(nuclei)
%     temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));  %could this miscount so that a nucleus that has 2 spots, then one for ex would be counted as having 3 spots?
%     num_spots=length(temp);
%     num_spots_ap(ii).spots=num_spots;
%     num_spots_ap(ii).MeanAP=CompiledParticles_14(temp(1)).MeanAP;
%     
end
%round AP positions of nuclei to closest AP bin
APEstm=[num_spots_ap.MeanAP];
for j=1:length(nucleithresh)
    APEstm(j)=round(num_spots_ap(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            num_spots_ap(j).APBin=APbinID(jj);
            break;
        end
    end
end

% Deal with nuclei with no spot 
%first translate x,y coordinates to AP position 
%Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
NucPosition=[];



for i=1:length(schnitzcells)     %find the AP position of each nucleus across time 
        for j=1:length(schnitzcells(i).frames)
            
            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((schnitzcells(i).ceny(j)-coordAZoom(2))./(schnitzcells(i).cenx(j)-coordAZoom(1)));
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-schnitzcells(i).ceny(j)).^2+(coordAZoom(1)-schnitzcells(i).cenx(j)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            NucTimeFrame=schnitzcells(i).frames(j);  %making the columns the time frame, each row is a nucleus 
            NucPosition(i,NucTimeFrame)=APPositions/APLength;
            
        end
end

 NoSpotNucLife=[];
  NoSpotPosition=[];
  for ii=1:length(NucPosition)    %only keep nuclei that exist >=20 frames
        temp=[length(schnitzcells(ii).frames)];
        NoSpotNucLife(ii)=temp;
    end
    for jj=1:length(NucPosition)
        for qq=1:size(NucPosition,2)
        if NoSpotNucLife(jj) >=20
            NoSpotPosition(jj,qq)=NucPosition(jj,qq);
        end
        end
    end
     

%remove Nuclei that identified as having spots in above 

    SpotNuc=unique([CompiledParticles_14.Nucleus]); %find the nuclei that have spots. Only want to remove those that have spots in nc14 (??) 
    for jj=1:length(SpotNuc)
        NoSpotPosition(SpotNuc(jj),:)=nan;      %make sure these aren't used later in plots
        
    end
   
  
    %Isolate nuclei in nc14
   NoSpotPosition_nc14=NoSpotPosition(:,nc14:end); 
    %change 0's indicating nuclei not present at that time point to nan's
    %so don't have all off nuclei mapping to first AP bin when take average
    %of AP position
    %for ii=1:length(NoSpotPosition_nc14)
        NoSpotPosition_nc14(NoSpotPosition_nc14==0)=nan; %for ii=1:length(NoSpotPosition)
        
        
   
    
            
    %incorporate no-spot nuclei into array of nuclei with spots and change
    %AP position to their mean positions as did for nuclei w spots
    for j=1:length(NoSpotPosition_nc14)
        MeanAPNoSpot=nanmean(NoSpotPosition_nc14(j,:));
    num_spots_ap(length(num_spots_ap)+1).MeanAP=MeanAPNoSpot;%[mean(NoSpotPosition(j,:),2,'omitnan')]
    end
    %fill in 0 values for no spot nuclei
    for jj=1:length(num_spots_ap)
         if isempty(num_spots_ap(jj).spots)
         num_spots_ap(jj).spots =0;
        end
    end
    
    %round AP positions of nuclei to closest AP bin
APEstm=[num_spots_ap.MeanAP];
for j=1:length(num_spots_ap)
    APEstm(j)=round(num_spots_ap(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            num_spots_ap(j).APBin=APbinID(jj);
            break;
        end
    end
end
    

    %put number of spots in a nuclei as column and AP position of
    %associated nuclei as row to make stacked bar graph 
    APbins = unique([num_spots_ap.APBin]);   %need to add part to add up all the zero spot nuclei, also make sure only the double counted nuclei were changed to nan??
SpotsAP=nan(length(APbinID),3);%SpotsAP=zeros(length(APbins),3);
for ii=1:length(APbins)
    jj = find([num_spots_ap.APBin]==APbins(ii));
    APcolumn=find([APbinID]==APbins(ii));        %want row number to correspond to column number of AP bin in array of all AP bins across embryo
    
    Spotnumber=[];
   
    Spotnumber=[num_spots_ap(jj).spots];
    OneSpot=sum([Spotnumber]==1);
    TwoSpot=sum([Spotnumber]==2);
    NoSpot=sum([Spotnumber]==0);
    SpotsAP(APcolumn,1)=NoSpot;
    SpotsAP(APcolumn,2)=OneSpot;
    SpotsAP(APcolumn,3)=TwoSpot;
    
end
figure 
bar(APbinID,SpotsAP,'stacked')
xlabel('Egg Length'); ylabel('# of nuclei')
legend('no spot','one spot','two spots','Location','best');
title('ON nuclei nc14');
ylim([0 60]); 
xlim([APbins(1)-0.25 APbins(end)+0.25]);

save('SpotsAPRatio10Frame','SpotsAP');
saveas(gcf,'SpotOnRatioAP10Frame.jpg');
%% nc13 
nc_number=[CompiledParticles.nc];
CompiledParticles_13=CompiledParticles(nc_number==13);
nuclei=unique([CompiledParticles_13.Nucleus]);

for ii=1:length(nuclei)   %Find out how long each nucleus exists for 
    temp=find([CompiledParticles_13.Nucleus]==nuclei(ii));
    temp2=length(schnitzcells(temp(1)).frames);
    NucleiLife(ii)=temp2;
%     
end

for jj=1:length(nuclei)
    if NucleiLife(jj) >=20
        nucleithresh(jj)=nuclei(jj);
    end
end
nucleithresh=unique([nucleithresh]); %Want to get rid of 0's
nucleithresh=[nucleithresh(2:end)];

for ii=1:length(nucleithresh)    %only want to look at nuclei that exist >= 20 frames
    temp=find([CompiledParticles_13.Nucleus]==nucleithresh(ii));
    num_spots=length(temp);
    num_spots_ap_13(ii).spots=num_spots;
    num_spots_ap_13(ii).MeanAP=CompiledParticles_13(temp(1)).MeanAP;
% for ii=1:length(nuclei)
%     temp=find([CompiledParticles_13.Nucleus]==nuclei(ii));  %could this miscount so that a nucleus that has 2 spots, then one for ex would be counted as having 3 spots?
%     num_spots=length(temp);
%     num_spots_ap_13(ii).spots=num_spots;
%     num_spots_ap_13(ii).MeanAP=CompiledParticles_13(temp(1)).MeanAP;
    
end

%round AP positions of nuclei to closest AP bin
APEstm=[num_spots_ap_13.MeanAP];
for j=1:length(nuclei)
    APEstm(j)=round(num_spots_ap_13(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            num_spots_ap_13(j).APBin=APbinID(jj);
            break;
        end
    end
end

% Deal with nuclei with no spot 
%first translate x,y coordinates to AP position 
%Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
NucPosition_13=[];



for i=1:length(schnitzcells)     %find the AP position of each nucleus across time 
        for j=1:length(schnitzcells(i).frames)
            
            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((schnitzcells(i).ceny(j)-coordAZoom(2))./(schnitzcells(i).cenx(j)-coordAZoom(1)));
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-schnitzcells(i).ceny(j)).^2+(coordAZoom(1)-schnitzcells(i).cenx(j)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            NucTimeFrame=schnitzcells(i).frames(j);  %making the columns the time frame, each row is a nucleus 
            NucPosition_13(i,NucTimeFrame)=APPositions/APLength;
            
        end
end

 NoSpotNucLife=[];
%  NoSpotPosition=NucPosition;
  for ii=1:length(NucPosition)    %only keep nuclei that exist >=20 frames
        temp=[length(schnitzcells(ii).frames)];
        NoSpotNucLife(ii)=temp;
    end
    for jj=1:length(NucPosition)
        if NoSpotNucLife(jj) >=20
            NoSpotPosition(jj)=NucPosition(jj);
        end
    end
    
%remove Nuclei that identified as having spots in above 
NoSpotPosition_13=NucPosition_13;
    SpotNuc=unique([CompiledParticles.Nucleus]); %find the nuclei that have spots 
    for jj=1:length(SpotNuc)
        NoSpotPosition_13(SpotNuc(jj),:)=nan;      %make sure these aren't used later in plots

        
    end
    
     %Isolate nuclei in nc14
     if nc13==0
         nc13time=1;
     else 
         nc13time=nc13;
     end
     
   NoSpotPosition_nc13=NoSpotPosition_13(:,nc13time:end); 
    %change 0's indicating nuclei not present at that time point to nan's
    %so don't have all off nuclei mapping to first AP bin when take average
    %of AP position
    %for ii=1:length(NoSpotPosition_nc14)
        NoSpotPosition_nc13(NoSpotPosition_nc13==0)=nan; 
    
 %incorporate no-spot nuclei into array of nuclei with spots and change
    %AP position to their mean positions as did for nuclei w spots
    for j=1:length(NoSpotPosition_nc13)
        MeanAPNoSpot=nanmean(NoSpotPosition_nc13(j,:));
    num_spots_ap_13(length(num_spots_ap_13)+1).MeanAP=MeanAPNoSpot;%[mean(NoSpotPosition(j,:),2,'omitnan')]
    end
    %fill in 0 values for no spot nuclei
    for jj=1:length(num_spots_ap_13)
         if isempty(num_spots_ap_13(jj).spots)
         num_spots_ap_13(jj).spots =0;
        end
    end
    
APEstm=[num_spots_ap_13.MeanAP];
for j=1:length(num_spots_ap_13)
    APEstm(j)=round(num_spots_ap_13(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            num_spots_ap_13(j).APBin=APbinID(jj);
            break;
        end
    end
end

   APbins = unique([num_spots_ap_13.APBin]);   %need to add part to add up all the zero spot nuclei, also make sure only the double counted nuclei were changed to nan??
SpotsAP=zeros(length(APbins),3);
for ii=1:length(APbins)
    jj = find([num_spots_ap_13.APBin]==APbins(ii));
    
    Spotnumber=[];
   
    Spotnumber=[num_spots_ap_13(jj).spots];
    OneSpot=sum([Spotnumber]==1);
    TwoSpot=sum([Spotnumber]==2);
    NoSpot=sum([Spotnumber]==0);
    SpotsAP(ii,1)=NoSpot;
    SpotsAP(ii,2)=OneSpot;
    SpotsAP(ii,3)=TwoSpot;
    
end
figure 
bar(APbins,SpotsAP,'stacked')
xlabel('Egg Length'); ylabel('# of nuclei')
legend('no spot','one spot','two spots','Location','best');
title('ON nuclei nc13');
ylim([0 80]);
    
    
    
