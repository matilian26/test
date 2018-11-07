%% Plot average fraction of ON nuclei in different nc 13 & nc 14
%% Fraction of ON nuclei (out of all nuclei) for given AP position vs Time 
DataSet=input('What data set to evaluate?', 's');
Data=LoadMS2Sets(DataSet);
NEllipses(:,:,1)=Data(1).NEllipsesAP;
for i=1:length(Data)
    NEllipses(1:length(Data(i).NEllipsesAP),:,i)=Data(i).NEllipsesAP;
    NParticles(1:length(Data(i).NEllipsesAP),:,i)=Data(i).NParticlesAP;
end
for i=1:length(Data)
    AllON(:,:,i)=NParticles(:,:,i) ./ NEllipses(:,:,i);
end
MeanAllON=nansum(AllON,3) ./ length(Data);
APbin=input('Which AP bin to evaluate?');
APbinName=num2str(APbin);
figure 
plot(MeanAllON(:,APbin))
xlabel('Time'); 
ylabel('Avg fraction of ON nuclei');
title(['ON nuclei AP bin',APbinName, DataSet]);
ylim([0 2.5]);

figure 
plot(MeanAllON)
xlabel('Time');
ylabel('Avg fraction of ON nuclei');
title(['ON nuclei all AP bins','',DataSet]);
ylim([0 2.5]);


%% Fraction of ON nuclei for a given time (hopefully I'll get this to be a certain nc) 
%% nc 13 
NEllipses=[];
NParticles=[];
for i=1:length(Data)
    NC13Start(i)=[Data(i).nc13].';
    NC13End(i)=[Data(i).nc14].'-1;
end
NC13length=NC13End-NC13Start;
NEllipses(:,:,1)=Data(1).NEllipsesAP;
AllParticles=[];
for i=1:length(Data)
    for n=1:length(Data(i).NParticlesAP)
        if isempty(Data(i).NParticlesAP(n))
            AllParticles(:,n,i)=0;
        end
    end
end

for i=1:length(Data)
    NEllipses(NC13Start(i):NC13length(i),:,i)=Data(i).NEllipsesAP(NC13Start(i):NC13length(i),:);
    NParticles(NC13Start(i):NC13length(i),:,i)=Data(i).NParticlesAP(NC13Start(i):NC13length(i),:);
end
for i=1:length(Data)
    AllON(:,:,i)=NParticles(:,:,i) ./ NEllipses(:,:,i);
end
MeanAllON=nansum(AllON,3) ./ length(Data);
% MeanONAllAP=MeanAllON';
% NC13length=max(NC13End-NC13Start);
% FirstNC13=min(NC13Start);
% LastNC13=max(NC13End);
% figure 
% plot(MeanONAllAP(:,FirstNC13:LastNC13))
% 

%% 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

for cc=1:length(ConstructList)
    ConRatioOn=[];
    Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        Filename=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat']
        Filename2=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat'];
        Filename3=[DropboxFolder filesep PrefixName filesep 'APDetection.mat'];
        load(Filename);
        load(Filename2);
        load(Filename3);
        RatioOn=nan(length(ElapsedTime),41);
        %specifically look at bursts in nc14
         nc_number=[CompiledParticles.nc];
         CompiledParticles_14=CompiledParticles(nc_number==14);
               % Put nucleus field in Schnitz 
        for ss=1:length(schnitzcells)
            schnitzcells(ss).Nucleus=ss;
        end
        %Limit to Nuclei that exist in nc14 
        nc14Frame=Data(ee).nc14;
        schnitzcells_n14=[schnitzcells];
        for ss=1:length(schnitzcells_n14)
            if (max(schnitzcells_n14(ss).frames)) < nc14Frame
                schnitzcells_n14(ss).Nucleus=nan;
               
            end
        end
        %Need to add AP bin info to schnitzcells
        
        %first translate x,y coordinates to AP position 
%Angle between the x-axis and the AP-axis
    if exist('coordPZoom', 'var')
        APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    else
        error('coordPZoom not defined. Was AddParticlePosition.m run?')
    end
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
NucPosition=[];



for i=1:length(schnitzcells_n14)     %find the AP position of each nucleus across time 
        for j=1:length(schnitzcells_n14(i).frames)
            
            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((schnitzcells_n14(i).ceny(j)-coordAZoom(2))./(schnitzcells_n14(i).cenx(j)-coordAZoom(1)));
            
            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-schnitzcells_n14(i).ceny(j)).^2+(coordAZoom(1)-schnitzcells_n14(i).cenx(j)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            NucTimeFrame=schnitzcells_n14(i).frames(j);  %making the columns the time frame, each row is a nucleus 
            NucPosition(i,NucTimeFrame)=APPositions/APLength;
            
        end
   %schnitzcells_nc14(i).APPos=[NucPosition(i,:)];     
end
for ii=1:length(NucPosition)
    schnitzcells_n14(ii).APPos=[NucPosition(ii,:)];
end
for ii=1:length(schnitzcells_n14)
    schnitzcells_n14(ii).MeanAP=nanmean(schnitzcells_n14(ii).APPos);
end
APbinID=[Data(ee).APbinID];                %estimate AP bin for later 
%APEstm=[schnitzcells_n14.MeanAP];
for j=1:length(schnitzcells_n14)
    APEstm(j)=round(schnitzcells_n14(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            schnitzcells_n14(j).APBin=APbinID(jj);
            break;
        end
    end
end

for j=1:length(CompiledParticles_14)
    APEstm(j)=round(CompiledParticles_14(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            CompiledParticles_14(j).APBin=APbinID(jj);
            break;
        end
    end
end
APParticleInfo=[CompiledParticles_14.APBin];
APNucInfo=[schnitzcells_n14.APBin];
for aa=1:length(APbinID)
    SpotAPsubset=CompiledParticles_14(APParticleInfo==APbinID(aa));
    NucAPsubset=schnitzcells_n14(APNucInfo==APbinID(aa));
    for tt=1:length(ElapsedTime)
        SpotAround=0;
        NucAround=0;
        for bb=1:length(NucAPsubset)
            if any(ismember(tt,NucAPsubset(bb).frames))
                NucAround=NucAround+1;
            end
        end
        for bb=1:length(SpotAPsubset)
            if any(ismember(tt,SpotAPsubset(bb).Frame))
                SpotAround=SpotAround+1;
            end
        end
    if SpotAround ~= 0 & (NucAround ~=0) 
    RatioOn(tt,aa)=SpotAround/NucAround;
    elseif (SpotAround==0) & (NucAround~=0)
        RatioOn(tt,aa)=0;
    elseif (SpotAround==0) & (NucAround==0)
        RatioOn(tt,aa)=nan;
    end
    end
end
ConRatioOn=[ConRatioOn; RatioOn];
WholeRatios(cc).Embryo(ee).RatioOn=RatioOn;
    end
    WholeRatios(cc).AllRatioOn=ConRatioOn;
end

            
% Plot stuff
for cc=1:length(ConstructList)
    plot(nanmean(WholeRatios(cc).AllRatioOn))
    hold on
end
    

