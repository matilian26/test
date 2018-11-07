%% 2 way ANOVA of spot correlation vs AP position & Construct 
ConstructList= {'KrBoth';'KrDist';'KrProx';'KrBothSep';'KrDistDuplicN';'KrProxDuplic'} %{'KrBoth';'KrDist';'KrProx';'KrProxAtDist';...
    %'KrProxDuplic';'KrDistAtProxN';'KrDistDuplicN'};
    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
 Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
 ] = readMovieDatabase('2017-08-03-mKr1_E1');    %just any random dataset to give us the dropbox folder location

Man=[];
Verz=input('Use adjusted correlation calculation?','s');


counter=0;
for q=1:length(ConstructList)
    Data= LoadMS2SetsCS(ConstructList{q});
    NEmbryos = length(Data);
    Label = q;   %in original matrix have constructs equal numbers so not trying to convert strings/doubles/cells

    for ii=1:NEmbryos
        PrefixName=Data(ii).Prefix;
        if Verz=='y'
        filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelationAdj.mat']; %adjusted corr values rn 
        else
           filename=[DropboxFolder filesep PrefixName filesep 'SpotCorrelation.mat']; %adjusted corr values rn 
        end 
        load(filename);
        for jj=1:length(SpotDiff)
            
            if ~isempty(SpotDiff(jj).SpotCorr)
                counter=counter+1;
                if length(SpotDiff(jj).SpotCorr)==1
                    Man(1,counter)=SpotDiff(jj).SpotCorr(1)
                else
            Man(1,counter)=SpotDiff(jj).SpotCorr(1,2);
                end
            Man(2,counter)=Label;
            Man(3,counter)=SpotDiff(jj).APBin;
            end
            
        end
    end
end    %created matrix where 1st row is correlation of each spot pair, 2nd row is corresponding construct (in number form), and 3rd row is corresponding AP bin. For 2way anova

SpotCorrValues=[Man(1,:)];
Constructs=[Man(2,:)];
APPosition=[Man(3,:)];
% ConstructsLetters=[];
% for qq=1:length(Constructs)
% if Constructs(qq)==1
% ConstructsLetters(qq)='Both';
% elseif Constructs(qq)==2
% ConstructsLetters(qq)='Dist';
% elseif Constructs(qq)==3
% ConstructsLetters(qq)='Prox';
% elseif Constructs(qq)==4
% ConstructsLetters(qq)='BothSep';
% end
% end

[p,tbl,stats]=anovan(SpotCorrValues,{Constructs,APPosition}, 'model',2,'varnames',{'Construct','AP Bin'});