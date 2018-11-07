%% create scatter plot of spot 1 fluo vs spot 2  
%load constructs
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic';'KrBoth'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location
%Ask if only want nc14 info   for right now just going to use nc14
ncUse=input('Want to only use nc14?','s');
NumberNuclei=nan(20,41,[length(ConstructList)]);
% go through each embryo of each construct
for cc=1:length(ConstructList)
     Data= LoadMS2SetsCS(ConstructList{cc});
    Datalength(cc)=length(Data);
    NEmbryos = length(Data);
    APbinID=[Data(1).APbinID];
    Label = ConstructList(cc);
    for aa=1:length(APbinID)
        WholeSpots(cc).APbin(aa).Spots=[];
    end
    for ee=1:NEmbryos
        PrefixName=Data(ee).Prefix;
        if ncUse=='y'
        Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat'];
        Filename2=[DropboxFolder filesep PrefixName filesep 'compiledparticles.mat'];
        else
            Filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelation.mat'];
        end
        load(Filename);
        load(Filename2);
        APstuff=[SpotDiff.APBin];
        for aa=1:length(APbinID)
            APsubset=[];
            APSpots=[];
            APsubset=SpotDiff(APstuff==APbinID(aa));
            if isempty(APsubset)
                APSpots=[nan,nan,nan];
            else
            for bb=1:length(APsubset)
                NucSpots=[];
                if isfield(APsubset, 'SpotTwo') & (length(APsubset(bb).SpotOne)==length(APsubset(bb).SpotTwo))
                    NucSpots=[[APsubset(bb).SpotOne]',[APsubset(bb).SpotTwo]',[1:length(APsubset(bb).SpotOne)]'];
                    APSpots=[APSpots;NucSpots];
                end
            end
            end
                
                WholeSpots(cc).Embryo(ee).APbin(aa).Spots=APSpots;
                WholeSpots(cc).APbin(aa).Spots=[WholeSpots(cc).APbin(aa).Spots;WholeSpots(cc).Embryo(ee).APbin(aa).Spots];
        end
        WholeSpots(cc).Embryo(ee).ElapsedTime=ElapsedTime;
    end
 
        
end

%% Plot Stuff
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
APtoUse=input('Which AP bin to use?');

EmbryoUse=input('Do scatter by embryos?','s');
if EmbryoUse=='y'
    
for cc=1:length(ConstructList)
    figure
    for ee=1:length(WholeSpots(cc).Embryo)
        c=(WholeSpots(cc).Embryo(ee).APbin(APtoUse).Spots(:,3).*[1 2.5 1.5])./255;
        %c=linspace(1,10,length(WholeSpots(cc).Embryo(ee).APbin(APtoUse).Spots(:,1)));
        cmap=pmkmp(length(WholeSpots(cc).Embryo(ee).ElapsedTime),'LinLhot')
        colormap(cmap)
        c=colorbar
        c.Label.String='Frames into nc14'
        scatter(WholeSpots(cc).Embryo(ee).APbin(APtoUse).Spots(:,1), WholeSpots(cc).Embryo(ee).APbin(APtoUse).Spots(:,2),[],WholeSpots(cc).Embryo(ee).APbin(APtoUse).Spots(:,3),'LineWidth',1.5)
        hold on
    end
    
    title(ConstructList{cc});
    h=refline(1,0);
        h.Color='k';
    ylabel('Spot 2 fluorescence (A.U.)');
    xlabel('Spot 1 fluorescence (A.U.)');
end
else
    for cc=1:length(ConstructList)
        figure 
        scatter(WholeSpots(cc).APbin(APtoUse).Spots(:,1),WholeSpots(cc).APbin(APtoUse).Spots(:,2));
        h=refline(1,0);
        h.Color='k';
        h.LineWidth=2;
        title(ConstructList{cc})
        ylabel('Spot 2 fluorescence');
        xlabel('Spot 1 fluorescence');
    end
end

BestFits=input('Do least square lines?','s');
if BestFits=='y'
    for cc=1:length(ConstructList)
        figure 
        scatter(WholeSpots(cc).APbin(APtoUse).Spots(:,1),WholeSpots(cc).APbin(APtoUse).Spots(:,2));
        h=lsline;
        h.Color='k';
        h.LineWidth=2;
        title(ConstructList{cc})
        ylabel('Spot 2 fluorescence');
        xlabel('Spot 1 fluorescence');
    end
end