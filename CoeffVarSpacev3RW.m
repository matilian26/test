%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty';'KrDist32C';'KrProx32C';'KrBothSep32C';'KrBoth32C';'Kr2xProx32C'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');

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
        Filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end
        Filename2=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat'];
        Filename3=[DropboxFolder filesep PrefixName filesep 'APDetection.mat'];
        load(Filename);
        load(Filename2);
        load(Filename3);
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
        APstuff=[BurstProperties.APBin];
        APSchnitstuff=[schnitzcells_n14.APBin];
        for aa=1:length(APbinID)
            APsubset=[];
            TotalmRNA=nan(1,600); %Go back and switch this to 1, 600
            APSchnitsubset=[];
%             VarPart=[];
%             SumVarPart=[];
            APmean=[];
            APStdDev=[];
            APsubset=BurstProperties(APstuff==APbinID(aa));
            APSchnitsubset=schnitzcells_n14(APSchnitstuff==APbinID(aa));
            if (isempty(APsubset)) & (isempty(APSchnitsubset))
                TotalmRNA(:)=nan;
            end
            if (isempty(APsubset)) & (~isempty(APSchnitsubset))
                for bb=1:length(APSchnitsubset)
                    TotalmRNA(bb)=0;
                end
            end
            if ~isempty(APsubset)
                
            for bb=1:length(APsubset)
                
                TotalmRNA(bb)=APsubset(bb).TotalmRNA;
                TotalmRNAerror(bb)=APsubset(bb).TotalmRNAError;
            end
                
            for bb=1:length(APSchnitsubset)
                TotalmRNA(bb+length(APsubset))=0;
            end
            end
            
            %if ~isnan(TotalmRNA)
            APmean=nanmean(TotalmRNA);
            APStdDev=nanstd(TotalmRNA);
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
            Numbernucs(ee,aa)=sum(~isnan(TotalmRNA));
            %end
            EmbTotalmRNA=[EmbTotalmRNA, (TotalmRNA)'];
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
            APTable(ee,aa)=nan;
        end
    end
    WholeAPTable(cc).APTable=APTable;
    WholeAPTable(cc).MeanCV=nanmean(WholeAPTable(cc).APTable);
    WholeAPTable(cc).AllNucNumbers=Numbernucs;
    WholeAPTable(cc).TotalmRNA=ConTotalmRNA;
    WholeAPTable(cc).AllCV=((nanstd(WholeAPTable(cc).TotalmRNA))./(nanmean(WholeAPTable(cc).TotalmRNA)));
    WholeAPTable(cc).AllNoiseStrength=((var(WholeAPTable(cc).TotalmRNA, 'omitnan'))./(nanmean(WholeAPTable(cc).TotalmRNA)));
    %Calculate variance in CV for each CV value 
    for aa=1:length(APbinID)
        CVVariance=[];
        %NumberNucsAP=sum(WholeAPTable(cc).AllNucNumbers(:,aa));
        NumberNucsAP=sum(~isnan(WholeAPTable(cc).TotalmRNA(:,aa)));
        %SqCV=(WholeAPTable(cc).AllCV(aa))^2;
        %CVVariance(aa)=(SqCV)*((1/(2*(NumberNucsAP-1)))+((1/((NumberNucsAP)^2))*SqCV));
        %CVVariance(aa)=(((WholeAPTable(cc).AllCV(aa))^2)*((1/(2*(NumberNucsAP-1)))+((1/((NumberNucsAP).^2))*((WholeAPTable(cc).AllCV(aa))^2))));
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
%%

%Plot noise vs mean
NoisevsMean=input('Plot noise vs mean?', 's');
if NoisevsMean=='y'
    for cc=1:length(ConstructList)
        figure
        for ee=1:length(WholeAPTable(cc).Embryo)
            scatter(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable), WholeAPTable(cc).Embryo(ee).EmbCV(:));
            %scatter(log2(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable)), log2(WholeAPTable(cc).Embryo(ee).EmbCV(:)));
            hold on 
        end
        xlabel('Mean total mRNA production');
        ylabel('Spatial noise')
        title(ConstructList{cc});
        xlim([0 1200000]);
        ylim([0 8]);
    figure
    for ee=1:length(WholeAPTable(cc).Embryo)
        scatter(nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable), WholeAPTable(cc).Embryo(ee).EmbNoiseStrength(:));
        hold on 
    end
    xlabel('Mean total mRNA production')
    title(ConstructList{cc});
    ylabel('Noise strength')
    end    
end

LSLine=input('Do least squares line?','s');
if LSLine=='y'
    for cc=1:length(ConstructList)
        EmbMeanProd=[];
        EmbCVVal=[];
        EmbNoiseStrengthVal=[];
        for ee=1:length(WholeAPTable(cc).Embryo)
            EmbMeanProd=[EmbMeanProd; nanmean(WholeAPTable(cc).Embryo(ee).EmbryoTable)];
            EmbCVVal=[EmbCVVal; WholeAPTable(cc).Embryo(ee).EmbCV];
            EmbNoiseStrengthVal=[EmbNoiseStrengthVal; WholeAPTable(cc).Embryo(ee).EmbNoiseStrength];
        end
        figure
        scatter(EmbMeanProd(:), EmbCVVal(:));
        xlabel('Mean total mRNA production')
        ylabel('Spatial noise');
        lsline
        title(ConstructList{cc})
        figure 
        scatter(EmbMeanProd(:), EmbNoiseStrengthVal(:));
        lsline
        title(ConstructList{cc});
        xlabel('Mean total mRNA production'); 
        ylabel('Noise strength');
    end
end
%Embryo by embryo data
Useembryos=input('Get embryo data?','s');
if Useembryos=='y'
    figure 
for cc=1:length(ConstructList)
    for ee=1:length(WholeAPTable(cc).Embryo)
        plot(EggLength,WholeAPTable(cc).Embryo(ee).EmbCV,'o', 'Color', Colors(cc).Color)
        hold on
    end
end

for cc=1:length(ConstructList)
    figure
    for ee=1:length(WholeAPTable(cc).Embryo)
        errorbar(EggLength,WholeAPTable(cc).Embryo(ee).EmbCV, WholeAPTable(cc).Embryo(ee).CVSE, 'o','Color',Colors(cc).Color);
        hold on 
    end
    errorbar(EggLength, WholeAPTable(cc).AllCV, WholeAPTable(cc).CVSE, 'Color', Colors(cc).Color, 'LineWidth', 2.5);
    xlabel('% Egg length')
    xlim([0 100]);
    ylabel('Spatial noise (CV)'); 
    ylim([0 8])
    title(ConstructList{cc});
end
end
%%
% Plot all CVs vs AP
figure
set(gca, 'FontSize', fontsize, 'FontName', fontname);
for cc=1:length(ConstructList)
    errorbar(EggLength, WholeAPTable(cc).AllCV, WholeAPTable(cc).CVSE, 'Color', Colors(cc).Color, 'LineWidth',1.5);
    hold on 
end
xlabel('% Egg length');
ylabel('Coefficient of variation');
xlim([0 100]);
%%
% Singles 
figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
xlabel('% egg length');
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 7]);
print('-painters',[FigDirect filesep 'ProxDistSpaceCV'],'-dsvg');


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
ylim([0 7]);
print('-painters',[FigDirect filesep 'SinglesBothSpaceCV'],'-dsvg');

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
errorbar(EggLength, WholeAPTable(3).AllCV, WholeAPTable(3).CVSE, 'Color', Colors(3).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% Egg length');
legend('Both Sep', 'Both');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('Coefficient of variation');
xlim([0 100]);
%%
%Doubles
figure 
errorbar(EggLength, WholeAPTable(1).AllCV, WholeAPTable(1).CVSE, 'Color', Colors(1).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
errorbar(EggLength, WholeAPTable(7).AllCV, WholeAPTable(7).CVSE, 'Color', Colors(7).Color, 'LineWidth', 2.5);
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% egg length');
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 7]);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
print('-painters',[FigDirect filesep 'AllSpaceCV'],'-dsvg');


figure 
errorbar(EggLength, WholeAPTable(6).AllCV, WholeAPTable(6).CVSE, 'Color', Colors(6).Color, 'LineWidth', 2.5);
hold on
errorbar(EggLength, WholeAPTable(7).AllCV, WholeAPTable(7).CVSE, 'Color', Colors(7).Color, 'LineWidth', 2.5);
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
xlabel('% egg length');
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 7]);
set(gca, 'FontSize', fontsize, 'FontName', fontname);
print('-painters',[FigDirect filesep '2xEnhnacersvBothSpaceCV'],'-dsvg');

%Zoom in to center 20%
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
ylim([0 7]);
print('-painters',[FigDirect filesep 'DistTCompSpaceCV'],'-dsvg');

figure 
errorbar(EggLength, WholeAPTable(2).AllCV, WholeAPTable(2).CVSE, 'Color', Colors(2).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(11).AllCV, WholeAPTable(11).CVSE, 'Color', Colors(11).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 7]);
print('-painters',[FigDirect filesep 'ProxTCompSpaceCV'],'-dsvg');

figure 
errorbar(EggLength, WholeAPTable(3).AllCV, WholeAPTable(3).CVSE, 'Color', Colors(3).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(12).AllCV, WholeAPTable(12).CVSE, 'Color', Colors(12).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
%legend('Both Sep', 'Both Sep 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 7]);

figure 
errorbar(EggLength, WholeAPTable(8).AllCV, WholeAPTable(8).CVSE, 'Color', Colors(8).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(13).AllCV, WholeAPTable(13).CVSE, 'Color', Colors(13).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
%legend('Both', 'Both 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 7]);
print('-painters',[FigDirect filesep 'BothTCompSpaceCV'],'-dsvg');


figure 
errorbar(EggLength, WholeAPTable(7).AllCV, WholeAPTable(7).CVSE, 'Color', Colors(7).Color, 'LineWidth', 2.5);
hold on 
errorbar(EggLength, WholeAPTable(14).AllCV, WholeAPTable(14).CVSE, 'Color', Colors(14).Color, 'LineWidth', 2.5,'LineStyle','-.');
xlabel('% egg length');
%legend('2x Proximal', '2x Proximal 32C');
set(gca, 'FontSize', fontsize, 'FontName', fontname);
ylabel('coefficient of variation');
xlim([0 100]);
ylim([0 7]);
print('-painters',[FigDirect filesep '2xProxTCompSpaceCV'],'-dsvg');
%%
% plot CV's with box plots of the CVs calculated per embryo
for cc=1:length(ConstructList)
    figure 
    errorbar(WholeAPTable(cc).AllCV, WholeAPTable(cc).CVSE, 'Color', Colors(cc).Color, 'LineWidth',1.5);
    hold on 
    boxplot(WholeAPTable(cc).APTable);
    xlabel('AP bin');
    ylabel('Spatial noise (CV)');
end
