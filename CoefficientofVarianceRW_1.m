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
    APTable=zeros(50,41);
    for ee=1:NEmbryos
        
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'BurstPropertiesnc14.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'BurstProperties.mat'];
        end
        load(Filename);
        %if ncUse=='y'
%         nc_number=[CompiledParticles.nc];
%         CompiledParticles_14=CompiledParticles(nc_number==14);
        APstuff=[BurstProperties.APBin];
        for aa=1:length(APbinID)
            APsubset=[];
            VarPart=[];
            SumVarPart=[];
            APmean=[];
            APsubset=BurstProperties(APstuff==APbinID(aa));
            if ~isempty(APsubset)
            for bb=1:length(APsubset)
                
                TotalmRNA(bb)=APsubset(bb).TotalmRNA;
                TotalmRNAerror(bb)=APsubset(bb).TotalmRNAError;
            end
            else
                    TotalmRNA=nan;
                end
        
            if ~isnan(TotalmRNA)
            APmean=nanmean(TotalmRNA);
            %APstdDev=((nanstd(TotalmRNA))/(sqrt(length(TotalmRNA))));   %Hur
%             else
%                 APmean=nan;
%             end
            for bb=1:length(TotalmRNA)
                VarPart(bb)=((TotalmRNA(bb)-APmean)^2);
            end
            end
            if ~isempty(VarPart)
            SumVarPart=sum(VarPart);
            APTable(ee,aa)=((sqrt((1/length(APsubset))*(SumVarPart)))/APmean);
            
            end
        end
    end
    APTable(APTable==0)=nan;
    
    WholeAPTable(cc).APTable=APTable;
end

%% Visualizing 
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
EgglengthVector=APbinID .* 100;

for cc=1:length(ConstructList)
    for ee=1:size(WholeAPTable(cc).APTable,1)
        APBoxing(ee,cc)= WholeAPTable(cc).APTable(ee,WhichAP);
    end
end
APBoxing(APBoxing==0)=nan;
[p,tbl,stats]=anova1(APBoxing);

figure 
for cc=1:length(ConstructList)
    for ee=1:size(WholeAPTable(cc).APTable,1)
        plot(cc,[WholeAPTable(cc).APTable(ee,WhichAP)],'o','Color',Colors(cc).Color);
        hold on 
    end
end
boxplot(APBoxing,'Colors','k')
title(['Relative expression noise in space',' ', num2str(EgglengthVector(WhichAP)), '% egg length']);
xticks([1:6]);
xticklabels({'Distal','Proximal','Both Sep', '2x Distal', '2x Proximal', 'Both'});
xlabel('Construct')
ylabel('Coefficient of variation');
xlim([0 7]);

%plot for each construct across AP the mean spatial variance at each point 
figure
for cc=1:length(ConstructList)
   
    MeanSpatialVariance(cc).MeanVar =nanmean(WholeAPTable(cc).APTable);
%     figure
    plot(EgglengthVector,MeanSpatialVariance(cc).MeanVar,'Color',Colors(cc).Color,'LineWidth',1.5);
    hold on 
end
xlim([0 100])
xlabel('% Egg length');
title('Mean spatial variation');
legend('Distal','Proximal','Both Sep', '2x Distal','2x Proximal', 'Both');

%Compare SE to singles
figure
plot(EgglengthVector,MeanSpatialVariance(1).MeanVar,'Color',Colors(1).Color,'LineWidth',1.5)
hold on 
plot(EgglengthVector,MeanSpatialVariance(2).MeanVar,'Color',Colors(2).Color,'LineWidth',1.5)
plot(EgglengthVector,MeanSpatialVariance(6).MeanVar,'Color',Colors(6).Color,'LineWidth',1.5)
xlim([0 100])
xlabel('% Egg length');
title('Relative expression noise in space');
legend('Distal','Proximal','Both');
ylabel('Coefficient of variation');

% Compare SE to 2x constructs
figure
plot(EgglengthVector,MeanSpatialVariance(4).MeanVar,'Color',Colors(4).Color,'LineWidth',1.5)
hold on 
plot(EgglengthVector,MeanSpatialVariance(5).MeanVar,'Color',Colors(5).Color,'LineWidth',1.5)
plot(EgglengthVector,MeanSpatialVariance(6).MeanVar,'Color',Colors(6).Color,'LineWidth',1.5)
xlim([0 100])
xlabel('% Egg length');
title('Relative expression noise in space');
legend('2x Dist','2x Prox','Both');
ylabel('Coefficient of variation');

figure
plot(EgglengthVector,MeanSpatialVariance(4).MeanVar,'Color',Colors(4).Color,'LineWidth',1.5)
hold on 
plot(EgglengthVector,MeanSpatialVariance(6).MeanVar,'Color',Colors(6).Color,'LineWidth',1.5)
xlim([0 100])
xlabel('% Egg length');
title('Relative noise in space');
legend('2x Distal','Both');
ylabel('Coefficient of variation');

figure
plot(EgglengthVector,MeanSpatialVariance(5).MeanVar,'Color',Colors(5).Color,'LineWidth',1.5)
hold on
plot(EgglengthVector,MeanSpatialVariance(6).MeanVar,'Color',Colors(6).Color,'LineWidth',1.5)
xlim([0 100])
xlabel('% Egg length');
ylabel('Coefficient of variation');
title('Relative noise in space');
legend('2x Proximal', 'Both');

figure
plot(EgglengthVector,MeanSpatialVariance(3).MeanVar,'Color',Colors(3).Color,'LineWidth',1.5)
hold on 
plot(EgglengthVector,MeanSpatialVariance(6).MeanVar,'Color',Colors(6).Color,'LineWidth',1.5)
xlim([0 100])
xlabel('% Egg length');
ylabel('Coefficient of variation');
title('Relative noise in space');
legend('Both Sep', 'Both');

% Compare doubles to singles bc why not
figure 
plot(EgglengthVector,MeanSpatialVariance(1).MeanVar,'Color',Colors(1).Color,'LineWidth',1.5)
hold on 
plot(EgglengthVector,MeanSpatialVariance(4).MeanVar,'Color',Colors(4).Color,'LineWidth',1.5)
xlim([0 100])
xlabel('% Egg length');
ylabel('Coefficient of variation');
ylim([0 3])
title('Relative noise in space');
legend('Distal', '2x Distal');

figure 
plot(EgglengthVector,MeanSpatialVariance(2).MeanVar,'Color',Colors(2).Color,'LineWidth',1.5)
hold on 
plot(EgglengthVector,MeanSpatialVariance(5).MeanVar,'Color',Colors(5).Color,'LineWidth',1.5)
xlim([0 100])
xlabel('% Egg length');
ylim([0 3])
ylabel('Coefficient of variation');
title('Relative noise in space');
legend('Proximal', '2x Proximal');

%Save spatial noise info for the constructs 
save([DropboxFolder filesep 'Constructs' filesep 'AllSpatialNoise'],'WholeAPTable');
save([DropboxFolder filesep 'Constructs' filesep 'MeanSpatialNoise'],'MeanSpatialVariance');


% %ANOVA across space for each construct?
% for cc=1%:length(ConstructList)
%     [p,tbl,stats]=anova1(WholeAPTable(cc).APTable);
% end
