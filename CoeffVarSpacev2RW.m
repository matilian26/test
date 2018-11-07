%% calculate coefficient of variance for each construct 
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
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
    Numbernucs=nan(50,41);
    APTable=nan(50,41);
    for ee=1:NEmbryos
        nc14Frame=[];
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
APEstm=[schnitzcells_n14.MeanAP];
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
            TotalmRNA=[];
            APSchnitsubset=[];
%             VarPart=[];
%             SumVarPart=[];
            APmean=[];
            APStdDev=[];
            APsubset=BurstProperties(APstuff==APbinID(aa));
            APSchnitsubset=schnitzcells_n14(APSchnitstuff==APbinID(aa));
            if (isempty(APsubset)) & (isempty(APSchnitsubset))
                TotalmRNA=nan;
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
            
            if ~isnan(TotalmRNA)
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
            end
        end
    end
    %remove points where only have one value for a whole AP bin
    for aa=1:length(APbinID)
        if sum(~isnan(APTable(:,aa)))==1
            APTable(:,aa)=nan;
        end
    end
    %Calculate variance in CV for each CV value 
    
    WholeAPTable(cc).APTable=APTable;
    WholeAPTable(cc).MeanCV=nanmean(WholeAPTable(cc).APTable);
    WholeAPTable(cc).AllNucNumbers=Numbernucs;
    
end

%% Plotting
DistalColor=[8 180 238] ./ 255;
DoubDistColor=[1 17 181] ./ 255;
ProxColor=[251 230 60] ./ 255;
DoubProxColor=[251 190 80] ./ 255;
BothSepColor=[94 250 81] ./ 255;
BothColor=[12 195 82] ./ 255;

Colors(1).Color=DistalColor;
Colors(2).Color=ProxColor;
Colors(3).Color=BothSepColor;
Colors(4).Color=DoubDistColor;
Colors(5).Color=DoubProxColor;
Colors(6).Color=BothColor;

WhichAP=input('Which AP bin to use?');
EggLength=APbinID .* 100;

% Plot all across AP
figure 
for cc=1:length(ConstructList)
    plot(EggLength,WholeAPTable(cc).MeanCV,'Color',Colors(cc).Color,'LineWidth',1.5)
    hold on
end
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('Distal', 'Proximal', 'Both Sep', '2x Distal', '2x Proximal', 'Both');
xlim([0 100])

% Singles vs Both Sep or SE
figure 
plot(EggLength, WholeAPTable(1).MeanCV, 'Color', Colors(1).Color, 'LineWidth',1.5);
hold on 
plot(EggLength, WholeAPTable(2).MeanCV, 'Color', Colors(2).Color, 'LineWidth',1.5);
plot(EggLength, WholeAPTable(3).MeanCV, 'Color', Colors(3).Color, 'LineWidth',1.5);
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('Distal', 'Proximal', 'Both Sep');
xlim([0 100])

figure 
plot(EggLength, WholeAPTable(1).MeanCV, 'Color', Colors(1).Color, 'LineWidth',1.5);
hold on 
plot(EggLength, WholeAPTable(2).MeanCV, 'Color', Colors(2).Color, 'LineWidth',1.5);
plot(EggLength, WholeAPTable(6).MeanCV, 'Color', Colors(6).Color, 'LineWidth',1.5);
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('Distal', 'Proximal', 'Both');
xlim([0 100])

figure 
plot(EggLength, WholeAPTable(1).MeanCV, 'Color', Colors(1).Color, 'LineWidth',1.5);
hold on 
plot(EggLength, WholeAPTable(2).MeanCV, 'Color', Colors(2).Color, 'LineWidth',1.5);
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('Distal', 'Proximal');
xlim([0 100])

%Doubles vs SE or each other
figure 
plot(EggLength, WholeAPTable(4).MeanCV, 'Color', Colors(4).Color, 'LineWidth',1.5);
hold on 
plot(EggLength, WholeAPTable(5).MeanCV, 'Color', Colors(5).Color, 'LineWidth',1.5);
plot(EggLength, WholeAPTable(6).MeanCV, 'Color', Colors(6).Color, 'LineWidth',1.5);
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('2x Distal', '2x Proximal', 'Both');
xlim([0 100])

figure 
plot(EggLength, WholeAPTable(4).MeanCV, 'Color', Colors(4).Color, 'LineWidth',1.5);
hold on 
plot(EggLength, WholeAPTable(5).MeanCV, 'Color', Colors(5).Color, 'LineWidth',1.5);
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('2x Distal', '2x Proximal');
xlim([0 100])

figure 
plot(EggLength, WholeAPTable(4).MeanCV, 'Color', Colors(4).Color, 'LineWidth',1.5);
hold on 
plot(EggLength, WholeAPTable(6).MeanCV, 'Color', Colors(6).Color, 'LineWidth',1.5);
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('2x Distal', 'Both');
xlim([0 100])

figure 
plot(EggLength, WholeAPTable(5).MeanCV, 'Color', Colors(5).Color, 'LineWidth',1.5);
hold on 
plot(EggLength, WholeAPTable(6).MeanCV, 'Color', Colors(6).Color, 'LineWidth',1.5);
xlabel('% Egg length')
ylabel('Coefficient of variation')
title('Relative noise across space');
legend('2x Proximal','Both');
xlim([0 100])

%Quick anova of CV's across constructs at specified AP bin
APbox=nan(20,[length(ConstructList)]);
for cc=1:length(ConstructList)
    for ee=1:length(WholeAPTable(cc).APTable)
        APbox(ee,cc)=WholeAPTable(cc).APTable(ee,WhichAP);
    end
end
[p,tbl,stats]=anova1(APbox);

%% 2way anova of Duration vs AP position vs genotype 
CVManova=[];
APManova=[];
ConManova=[];
for cc=1:length(ConstructList)
    for aa=1:length(APbinID)
        for bb=1:size(WholeAPTable(cc).APTable,1)
            CVManova=[CVManova, WholeAPTable(cc).APTable(bb,aa)];
            APManova=[APManova, aa];
            ConManova=[ConManova,cc];
        end
    end
end



[p,tbl,stats]=anovan(CVManova,{APManova, ConManova},'sstype',1,'varnames',{'AP bin','Construct'})
%multcompare(stats,'Dimension',2)