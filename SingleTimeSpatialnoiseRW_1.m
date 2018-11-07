%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
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
    Numbernucs=nan(50,41);
    APTable=nan(50,41);
    for ee=1:NEmbryos
        nc14Frame=[];
        EmbTotalmRNA=[];
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesSlope.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end
        Filename2=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat'];
        Filename3=[DropboxFolder filesep PrefixName filesep 'APDetection.mat'];
        Filename4=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        load(Filename);
        load(Filename2);
        load(Filename3);
        load(Filename4);
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
        for aa=1:length(APbinID)
            Counter=0;
            APsubset=[];
            mRNATime=nan(1,600);
            APSchnitsubset=[];
%             VarPart=[];
%             SumVarPart=[];
            APmean=[];
            APStdDev=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            APSchnitsubset=schnitzcells_n14(APSchnitstuff==APbinID(aa));
            if (isempty(APsubset)) & (isempty(APSchnitsubset))
                mRNATime(:)=nan;
            end
            if (isempty(APsubset)) & (~isempty(APSchnitsubset))
                for bb=1:length(APSchnitsubset)
                    NucOnlyFrames=[];
                    NucOnlyFrames=APSchnitsubset(bb).frames;
                    NucFrameInterest=ismember(TimePoint, NucOnlyFrames);
                    if any(NucFrameInterest)==1
                        mRNATime(bb)=0;
                    end
                end
                   
            end
            if ~isempty(APsubset)
                
            for bb=1:length(APsubset)
                Counter=Counter+1;
               if (length(APsubset(bb).SpotOne)) >= TimePoint 
                mRNATime(bb)=APsubset(bb).SpotOne(TimePoint);
                TotalmRNAerror(bb)=APsubset(bb).Err1(TimePoint);
               else
                   mRNATime(bb)=nan;
                   TotalmRNAerror(bb)=nan;
               end
               
            end
            
            if isfield(APsubset, 'SpotTwo')
            for bb=1:length(APsubset)
                if (length(APsubset(bb).SpotTwo)) >= TimePoint%1%~isempty(APsubset(bb).SpotTwo(TimePoint)) & (~isnan(APsubset(bb).SpotTwo(TimePoint)))
                    mRNATime(bb+Counter)= APsubset(bb).SpotTwo(TimePoint);
                    TotalmRNAerror(bb+Counter)=APsubset(bb).Err2(TimePoint);
                else
                    mRNATime(bb+Counter)=nan;
                    TotalmRNAerror(bb+Counter)=nan;
                end
            end
            end
            for bb=1:length(APSchnitsubset)
                mRNATime(bb+(2*Counter))=0;
            end
            end
            
            %if ~isnan(TotalmRNA)
            APmean=nanmean(mRNATime);
            APStdDev=nanstd(mRNATime);
            %APstdDev=((nanstd(TotalmRNA))/(sqrt(length(TotalmRNA))));   %Hur
%             else
%                 APmean=nan;
%             end
%             for bb=1:length(TotalmRNA)
%                 VarPart(bb)=((TotalmRNA(bb)-APmean)^2);
%             end
            %end
%             if ~isempty(VarPart)
%             SumVarPart=sum(VarPart);
            %APTable(ee,aa)=((sqrt((1/length(APsubset))*(SumVarPart)))/APmean);
            APTable(ee,aa)=((APStdDev)/(APmean));
            Numbernucs(ee,aa)=sum(~isnan(mRNATime));
            %end
            EmbTotalmRNA=[EmbTotalmRNA, (mRNATime)'];
        end
        for aa=1:length(APbinID)
            if sum(~isnan(EmbTotalmRNA(:,aa)))==1
                EmbTotalmRNA(:,aa)=nan;
            end
        end
        WholeAPTable(cc).Embryo(ee).EmbryoTable=EmbTotalmRNA;
        WholeAPTable(cc).Embryo(ee).EmbCV=((nanstd(WholeAPTable(cc).Embryo(ee).EmbryoTable))./(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable)));
        WholeAPTable(cc).Embryo(ee).EmbNoiseStrength=((var(WholeAPTable(cc).Embryo(ee).EmbryoTable,'omitnan'))./(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable)));
        WholeAPTable(cc).Embryo(ee).NumberNucs=[Numbernucs(ee,:)];
        ConTotalmRNA=[ConTotalmRNA; EmbTotalmRNA];
    end
    %remove points where only have one value for a whole AP bin
    for aa=1:length(APbinID)
        if sum(~isnan(ConTotalmRNA(:,aa)))==1
            ConTotalmRNA(:,aa)=nan;
        end
    end
    WholeAPTable(cc).APTable=APTable;
    WholeAPTable(cc).MeanCV=nanmean(WholeAPTable(cc).APTable);
    WholeAPTable(cc).AllNucNumbers=Numbernucs;
    WholeAPTable(cc).TotalmRNA=ConTotalmRNA;
    WholeAPTable(cc).AllCV=((nanstd(WholeAPTable(cc).TotalmRNA))./(nanmean(WholeAPTable(cc).TotalmRNA)));
    WholeAPTable(cc).AllNoiseStrength=((var(WholeAPTable(cc).TotalmRNA, 'omitnan'))./(nanmean(WholeAPTable(cc).TotalmRNA)));
%     WholeAPTable(cc).TestSD=nanstd(WholeAPTable(cc).TotalmRNA);
%     for aa=1:length(APbinID)
%         WholeAPTable(cc).TestSE(aa)=((WholeAPTable(cc).TestSD(aa))/(sqrt(sum(~isnan(WholeAPTable(cc).TotalmRNA(:,aa))))));
%     end
    %Calculate variance in CV for each CV value 
    for aa=1:length(APbinID)
        CVVariance=[];
        %NumberNucsAP=sum(WholeAPTable(cc).AllNucNumbers(:,aa));
        NumberNucsAP=sum(~isnan(WholeAPTable(cc).TotalmRNA(:,aa)));
        %SqCV=(WholeAPTable(cc).AllCV(aa))^2;
        %CVVariance(aa)=(SqCV)*((1/(2*(NumberNucsAP-1)))+((1/((NumberNucsAP)^2))*SqCV));
        CVVariance(aa)=(((WholeAPTable(cc).AllCV(aa))^2)*((1/(2*(NumberNucsAP-1)))+((1/(NumberNucsAP))*((WholeAPTable(cc).AllCV(aa))^2))));
        WholeAPTable(cc).CVSD(aa)=sqrt(CVVariance(aa));
        if NumberNucsAP ~=0
        WholeAPTable(cc).CVSE(aa)=(WholeAPTable(cc).CVSD(aa))/(sqrt(NumberNucsAP));
        else
            WholeAPTable(cc).CVSE(aa)=nan;
        end
    end
    
    
end

%% Look at individual embryo data 
for cc=1:length(ConstructList)
    EmbAPbox=[];
    for ee=1:length(WholeAPTable(cc).Embryo)
        EmbAPbox=[EmbAPbox; WholeAPTable(cc).Embryo(ee).EmbCV];
    end
    ConAPbox(cc).APbox=EmbAPbox;
end
% CV variance of each embryo 
for cc=1:length(ConstructList)
    for ee=1:length(WholeAPTable(cc).Embryo)
        for aa=1:length(APbinID) 
            CVEmbVariance=[];
            NumbEmbAPNucs=WholeAPTable(cc).Embryo(ee).NumberNucs(aa);
            CVEmbVariance=((WholeAPTable(cc).Embryo(ee).EmbCV(aa))^2)*((1/(2*(NumbEmbAPNucs-1)))+((1/(NumbEmbAPNucs)^2)*((WholeAPTable(cc).Embryo(ee).EmbCV(aa))^2)));
            WholeAPTable(cc).Embryo(ee).CVSD(aa)=(CVEmbVariance)^2;
            if NumbEmbAPNucs ~=0
            WholeAPTable(cc).Embryo(ee).CVSE(aa)=(WholeAPTable(cc).Embryo(ee).CVSD(aa))/(sqrt(NumbEmbAPNucs)); 
            else
                WholeAPTable(cc).Embryo(ee).CVSE(aa)=nan;
            end
        end
    end
end  

%% Plotting
DistalColor=[1 64 172]./255;
DistalEmptyColor=[8 210 238] ./ 255; 
Distal32CColor=[118 180 238] ./ 255;
DoubDistColor=[73 184 253] ./ 255;
ProxColor=[238 123 23]./255;
ProxEmptyColor=[251 250 50] ./255;
DoubProxColor=[215 183 58] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothSep32CColor=[150 250 81] ./255;
BothColor=[52 119 71]./255;
Both32CColor=[120 195 82] ./ 255;
BothEmptyColor=[12 250 150] ./ 255;
DoubProx32CColor=[200 150 100] ./ 255;

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DistalEmptyColor;
Colors(5).Color=ProxEmptyColor;
Colors(6).Color=DoubDistColor;
Colors(7).Color=DoubProxColor;
Colors(8).Color=BothColor;
Colors(9).Color=BothEmptyColor;
Colors(10).Color=DistalColor;
Colors(11).Color=ProxColor;
Colors(12).Color=BothSepColor;
Colors(13).Color=BothColor;
Colors(14).Color=DoubProxColor;

fontsize=18;
fontname='Helvetica';
EggLength=(APbinID).*100;
FigDirect=[DropboxFolder filesep 'Figures'];
%% Do by embryos

%%
%Plot noise vs mean
NoisevsMean=input('Plot noise vs mean?', 's');
if NoisevsMean=='y'
%     Culurz=rand(1,41);
%     Culurz2=linspace(1,5,41);
%     Markerz={'o','+','*','.','x','s','d','^','v','>','<','p','h'}
%     Widthz=linspace(1,5,15);
%     APColors=ones(1,41);
c=[1:41]';
c=num2str(c);
c=cellstr(c);
ColorMap=[ones(41,1), zeros(41,1), ones(41,1)];
for aa=1:length(APbinID)
    ColorMap(aa,:)=[ColorMap(aa,1).*rand(1),ColorMap(aa,2).*rand(1), ColorMap(aa,3)*rand(1)];
end
    for cc=1:length(ConstructList)
        figure
        for ee=1:length(WholeAPTable(cc).Embryo)
            scatter(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable), WholeAPTable(cc).Embryo(ee).EmbCV(:));
            %scatter(log2(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable)), log2(WholeAPTable(cc).Embryo(ee).EmbCV(:)));
            hold on 
            %s(ee).scat.Marker=Markerz{ee};
        end
        xlabel('Mean total mRNA production');
        ylabel('Spatial noise')
        title(ConstructList{cc});
        xlim([0 1200000]);
        ylim([0 8]);
    figure
    for ee=1:length(WholeAPTable(cc).Embryo)
            scatter(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable), log2(WholeAPTable(cc).Embryo(ee).EmbCV(:)));
            hold on
            labelpoints((nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable)),(log2(WholeAPTable(cc).Embryo(ee).EmbCV)),c)
    end
        xlabel(['log mean fluorescence time',num2str(TimePoint)]);
        ylabel('log CV');
        title(ConstructList{cc});
    figure
    for ee=1:length(WholeAPTable(cc).Embryo)
        scatter(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable), WholeAPTable(cc).Embryo(ee).EmbNoiseStrength(:));
        hold on 
    end
    xlabel(['Mean fluorescence time',num2str(TimePoint)])
    title(ConstructList{cc});
    ylabel('Noise strength')
    
    figure
    scatter(nanmean(WholeAPTable(cc).TotalmRNA),log2(WholeAPTable(cc).AllCV),[],[1:41],'LineWidth',1.5);
    labelpoints(nanmean(WholeAPTable(cc).TotalmRNA),log2(WholeAPTable(cc).AllCV),c);
    end
end
%%
%plot CVs vs AP 
for cc=1:length(ConstructList)
    plot(EggLength, WholeAPTable(cc).AllCV,'Color',Colors(cc).Color);
    hold on
end
xlabel('% Egg length')
ylabel('Spatial noise');
%%
%Compare singles
figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('coefficient of variation');
xlim([0 100]);
print('-painters',[FigDirect filesep 'ProxDistPtSpaceCV'],'-dsvg');


figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(EggLength, WholeAPTable(3).AllCV, WholeAPTable(3).CVSE, 'Color', Colors(3).Color, 'LineWidth', 2.5);
xlabel('% egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);

figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
print('-painters',[FigDirect filesep 'SinglesBothPtSpaceCV'],'-dsvg');

%% single allele 
figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(4).AllCV, WholeAPTable(4).CVSE, 'Color', Colors(4).Color, 'LineWidth', 2.5);
xlabel('% Egg length');
legend('Distal', '1x Distal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Coefficient of variation');
xlim([0 100]);

figure 
errorbar(EggLength, WholeAPTable(4).AllCV, WholeAPTable(4).CVSE, 'Color', Colors(4).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% Egg length');
legend('1x Distal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Coefficient of variation');
xlim([0 100]);

figure 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(5).AllCV, WholeAPTable(5).CVSE, 'Color', Colors(5).Color, 'LineWidth', 2.5);
xlabel('% Egg length');
legend('Proximal', '1x Proximal');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Coefficient of variation');
xlim([0 100]);

figure 
errorbar(EggLength, WholeAPTable(5).AllCV, WholeAPTable(5).CVSE, 'Color', Colors(5).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% Egg length');
legend('1x Proximal','Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Coefficient of variation');
xlim([0 100]);

figure 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(9).AllCV, WholeAPTable(9).CVSE, 'Color', Colors(9).Color, 'LineWidth', 2.5);
xlabel('% Egg length');
legend('Both', '1x Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Coefficient of variation');
xlim([0 100]);
%%
figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(7).AllCV, WholeAPTable(7).CVSE, 'Color', Colors(7).Color, 'LineWidth', 2.5);
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% egg length');
ylabel('coefficient of variation');
xlim([0 100]);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
print('-painters',[FigDirect filesep 'AllPtSpaceCV'],'-dsvg');


figure 
errorbar(EggLength, WholeAPTable(6).AllCV, WholeAPTable(6).CVSE, 'Color', Colors(6).Color, 'LineWidth', 2.5);
hold on
errorbar(EggLength, WholeAPTable(7).AllCV, WholeAPTable(7).CVSE, 'Color', Colors(7).Color, 'LineWidth', 2.5); 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% egg length');
ylabel('coefficient of variation');
xlim([0 100]);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
print('-painters',[FigDirect filesep '2xEnhancersvBothPtSpaceCV'],'-dsvg');

% Zoom in to center 20%
figure 
errorbar(EggLength, WholeAPTable(6).AllCV, WholeAPTable(6).CVSE, 'Color', Colors(6).Color, 'LineWidth', 2.5);
hold on
errorbar(EggLength, WholeAPTable(7).AllCV, WholeAPTable(7).CVSE, 'Color', Colors(7).Color, 'LineWidth', 2.5); 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% egg length');
ylabel('coefficient of variation');
xlim([40 60]);
ylim([0 3]);
set(gca, 'FontSize', fontsize, 'FontName', fontname);

%%
% Temperature alteration
figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(10).AllCV, WholeAPTable(10).CVSE, 'Color', Colors(10).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
%legend('Distal', 'Distal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 12]);
print('-painters',[FigDirect filesep 'DistTCompPtSpaceCV'],'-dsvg');

figure 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(11).AllCV, WholeAPTable(11).CVSE, 'Color', Colors(11).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 12]);
print('-painters',[FigDirect filesep 'ProxTCompPtSpaceCV'],'-dsvg');


figure 
errorbar(EggLength, WholeAPTable(3).AllCV, WholeAPTable(3).CVSE, 'Color', Colors(3).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(12).AllCV, WholeAPTable(12).CVSE, 'Color', Colors(12).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
%legend('Both Sep', 'Both Sep 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 12]);

figure 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(13).AllCV, WholeAPTable(13).CVSE, 'Color', Colors(13).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
%legend('Both', 'Both 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 12]);
print('-painters',[FigDirect filesep 'BothTCompPtSpaceCV'],'-dsvg');


figure 
errorbar(EggLength, WholeAPTable(7).AllCV, WholeAPTable(7).CVSE, 'Color', Colors(7).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(14).AllCV, WholeAPTable(14).CVSE, 'Color', Colors(14).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
%legend('2x Proximal', '2x Proximal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
print('-painters',[FigDirect filesep '2xProxTCompPtSpaceCV'],'-dsvg');
