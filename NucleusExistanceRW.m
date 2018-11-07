function NucleiLife=NucleusExistanceRW

%% Load embryos of constructs 
ConstructList= {'KrDist';'KrProx';'KrBothSep';'KrDistEmpty';'KrProxEmpty';'KrDistDuplicN';'KrProxDuplic';'KrBoth';'KrBothEmpty'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1'); 
for cc=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{cc})
    APbinID=[Data(1).APbinID];
    NEmbryos = length(Data);
    ConNucLife=[];
    nuclei=[];
    for ee=1:NEmbryos
        CompiledParticles=Data(ee).CompiledParticles;
        nc14=Data(ee).nc14;
         PrefixName=Data(ee).Prefix;
        filename=[DropboxFolder filesep PrefixName filesep 'APDetection.mat'];
        filename2=[DropboxFolder filesep PrefixName filesep PrefixName '_lin.mat'];
        load(filename); 
        load(filename2);
% Count number of frames a nucleus exists for in nc14
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
for bb=1:length(CompiledParticles_14)
    APEstm(bb)=round(CompiledParticles_14(bb).MeanAP,2);
    for aa=1:length(APbinID)
        if APEstm(bb) < APbinID(aa)
            CompiledParticles_14(bb).APBin=APbinID(aa);
            break;
        end
    end
end
nuclei=unique([CompiledParticles_14.Nucleus]);
NucleiLife=[];
for ss=1:length(schnitzcells)
            schnitzcells(ss).Nucleus=ss;
        end
schnitzcells_14=schnitzcells;
        for ss=1:length(schnitzcells_14)
            if (max(schnitzcells_14(ss).frames)) < nc14
                schnitzcells_14(ss).Nucleus=nan;
            end
            FramesUse=schnitzcells_14(ss).frames >= nc14;
            schnitzcells_14(ss).frames=schnitzcells_14(ss).frames(FramesUse);            
        if isempty(schnitzcells_14(ss).frames)
            schnitzcells_14(ss).frames=nan;
        end
    end

NucleiLife=nan(length(nuclei),41);

for qq=1:length(nuclei)
    temp=find([CompiledParticles_14.Nucleus]==nuclei(qq));
    APofNuc=CompiledParticles_14(temp(1)).APBin;
    aa=(APbinID==APofNuc);
    aa=find(aa);
    temp2=length(schnitzcells_14(nuclei(qq)).frames);
    NucleiLife(qq,aa)=temp2;
     
end
   
WholeExistance(cc).Embryo(ee).NucleiLife=NucleiLife;
WholeExistance(cc).Embryo(ee).TotalNuclei=sum(~isnan(NucleiLife));
ConNucLife=[ConNucLife; NucleiLife];
    end
    WholeExistance(cc).AllNucLife=ConNucLife;
    WholeExistance(cc).TotalNuclei=sum(~isnan(ConNucLife));
    WholeExistance(cc).AvgNucs=WholeExistance(cc).TotalNuclei ./ (length(WholeExistance(cc).Embryo));
end

%% Plot it
for cc=1:length(ConstructList)
    figure
    histogram(WholeExistance(cc).AllNucLife)
    title(ConstructList{cc});
    ylabel('# of nuclei');
    xlabel('# frames nuclei present');
end

    
    
    
 
%% Limit to nuclei above threshold of 20 frames 
%Find nuclei who exist for >=20 frames 
for jj=1:length(nuclei)
    if NucleiLife(jj) >=20
        nucleithresh(jj)=nuclei(jj);
    end
end
nucleithresh=unique([nucleithresh]); %Want to get rid of 0's
nucleithresh=[nucleithresh(2:end)]; %remove last "unique" 0 value 
for ii=1:length(nucleithresh)
    temp=find([CompiledParticles_14.Nucleus]==nucleithresh(ii));
end
    
        
        
        
        