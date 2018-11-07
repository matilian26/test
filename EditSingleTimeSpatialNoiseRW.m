%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');
TimePoint=input('Which time frame to use?');

% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    ConTotalmRNA=[];
    EmbCV=[];
    ConCV=[];
    Numbernucs=nan(50,41);
    APTable=nan(50,41);
    for ee=1:NEmbryos
        nc14Frame=[];
        %EmbCV=[];
        EmbTotalmRNA=[];
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end
        Filename2=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat'];
        Filename3=[DropboxFolder filesep PrefixName filesep 'APDetection.mat'];
        Filename4=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        Filename5=[DropboxFolder filesep PrefixName filesep 'CompiledParticles.mat'];
        load(Filename);
        load(Filename2);
        load(Filename3);
        load(Filename4);
        load(Filename5);
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
    if isempty(schnitzcells_n14(j).frames)
        schnitzcells_n14(j).APBin=nan;
    end
    if isnan(schnitzcells_n14(j).Nucleus)
        schnitzcells_n14(j).APBin=nan;
    end
    %Set threshold of nucleus existing for 5 frames- if doesn't meet,
    %change APbin to nan so not called in next section 
%     if length(schnitzcells_n14(j).frames) < 10
%         schnitzcells_n14(j).APBin=nan;
%     end
    % Remove schnitz that exist and produce mRNA in BurstingProperties 
    NucBursts=find([BurstProperties.Nucleus]==schnitzcells_n14(j).Nucleus);
    if ~isempty(NucBursts)
        schnitzcells_n14(j).APBin=nan;
    end
    
end

        %if ncUse=='y'
%         nc_number=[CompiledParticles.nc];
%         CompiledParticles_14=CompiledParticles(nc_number==14);
        APstuff=[SpotDiff.APBin];
        APSchnitstuff=[schnitzcells_n14.APBin];
        WholeAPCV=[];
        for aa=1:length(APbinID)
            Counter=0;
            APsubset=[];
            APSchnitsubset=[];
            APmean=[];
            APStdDev=[];
            APCV=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            APSchnitsubset=schnitzcells_n14(APSchnitstuff==APbinID(aa));
            if isempty(APsubset) & isempty(APSchnitsubset)
                APCV=nan(1,length(ElapsedTime));
            elseif isempty(APsubset) & ~isempty(APSchnitsubset)
                for tt=1:length(ElapsedTime)
                    for bb=1:length(APSchnitsubset)
                        if ismember(tt,[APSchnitsubset(bb).frames])
                            FluoVal(bb,tt)=0;
                        else
                            FluoVal(bb,tt)=nan;
                        end
                    end
                end
                if size(FluoVal,2) > ElapsedTime
                    FluoVal=FluoVal(:,1:length(ElapsedTime))
                end
                APCV=(nanstd(FluoVal))./(nanmean(FluoVal));
            elseif ~isempty(APsubset)
                for tt=1:length(ElapsedTime)
                    for bb=1:length(APsubset)
                        if length(APsubset(bb).SpotOne)==length(ElapsedTime)
                        FluoVal(bb,tt)=APsubset(bb).SpotOne(tt);
                        else
                            for rr=1:length(APsubset(bb).SpotOne)
                                FluoVal(bb,rr)=APsubset(bb).SpotOne(rr);
                            end
                            FluoVal(bb,[rr+1:length(ElapsedTime)])=nan;
                        end
                    end
                end
                for tt=1:length(ElapsedTime)
                    if isfield(APsubset,'SpotTwo')
                        for bb=1:length(APsubset)
                            if length(APsubset(bb).SpotTwo)==length(ElapsedTime)
                                FluoVal(bb+length(APsubset),tt)=APsubset(bb).SpotTwo(tt);
                            else
                            for rr=1:length(APsubset(bb).SpotTwo)
                                FluoVal(bb+length(APsubset),rr)=APsubset(bb).SpotTwo(rr);
                            end
                            FluoVal(bb+length(APsubset),[rr+1:length(ElapsedTime)])=nan;
                            
                            end
                        end
                    end
                end
           
                APCV=((nanstd(FluoVal))./(nanmean(FluoVal)));
            end
            WholeAPCV=[WholeAPCV,[APCV]'];
        end
        WholeNoise(cc).Embryo(ee).CV=[WholeAPCV];
        EmbCV=[EmbCV;WholeAPCV];
    end
    ConCV=[ConCV;EmbCV];
    WholeNoise(cc).CV=[EmbCV];
end

%% Plotting noise
DistalColor=[8 180 238] ./ 255;
DistalEmptyColor=[8 210 238] ./ 255; 
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 230 60] ./ 255;
ProxEmptyColor=[251 250 50] ./255;
DoubProxColor=[251 190 80] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothColor=[12 195 82] ./ 255;
BothEmptyColor=[12 250 150] ./ 255;

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DistalEmptyColor;
Colors(5).Color=ProxEmptyColor;
Colors(6).Color=DoubDistColor;
Colors(7).Color=DoubProxColor;
Colors(8).Color=BothColor;
Colors(9).Color=BothEmptyColor;

fontsize=18;
fontname='Helvetica';
EggLength=(APbinID).*100;

figure
for cc=1:length(ConstructList)
    plot(EggLength,nanmean(WholeNoise(cc).CV),'Color',Colors(cc).Color);
    hold on
end