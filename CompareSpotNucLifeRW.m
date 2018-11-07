%% compare the number of frames associated with each nucleus(spot) in compiledParticles to its corresponding nucleus in schnitzcells
    %want to know if lack of spot frame means nucleus was out of frame or
    %actually was there just no spot. consistent frames of schnitz > frames
    %of spot should indicate this is the (general) case that the nucleus is there
    %just not transcribing. 
% count number of nuclei nc14 with spots 
nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
nuclei=unique([CompiledParticles_14.Nucleus]);
NucleiLife=[];
NucSpot1Life=[];
NucSpot2Life=[];
TotalNucSpotLife=[];
for ii=1:length(nuclei)   %Find out how long each nucleus exists for 
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));
    NucSpot1Life(ii)=length(unique(CompiledParticles_14(temp(1)).Frame));
    if length(temp)==2
        NucSpot2Life(ii)=length(unique(CompiledParticles_14(temp(2)).Frame));
    end
end
if length(NucSpot1Life) > length(NucSpot2Life)
    NucSpot2Life(end:length(NucSpot1Life))=0;
else if length(NucSpot1Life) < length(NucSpot2Life)
        NucSpot1Life(end:length(NucSpot2Life))=0;
    end
end

TotalNucSpotLife=NucSpot1Life+NucSpot2Life;
%     if length(temp)==1
%         TotalNucSpotLife(ii)=NucSpot1Life(ii);
%     else if length(temp)==2
%             TotalNucSpotLife(ii)=NucSpot1Life(ii)+NucSpot2Life(ii);
%         end
%     end
    
    %temp2=length(schnitzcells(temp(1)).frames);
%     NucleiLife(ii)=temp2;
%     


for jj=nuclei(1):nuclei(end)
    nuctotallife(jj)=length(unique(schnitzcells(jj).frames));
end
AliveNoSpot=[];
for qq=1:length(nuclei)
    if nuctotallife(qq) > TotalNucSpotLife(qq)
        AliveNoSpot(qq,1)=1;
    else if nuctotallife(qq) <= TotalNucSpotLife(qq)
            AliveNoSpot(qq,2)=1;
        end
    end
end
OnFrameNoSpot=sum(AliveNoSpot(:,1));
OutFrame=sum(AliveNoSpot(:,2));
AliveorNot=[OnFrameNoSpot OutFrame];
Case={'In Frame','Exit Frame'};
Name=pwd;
Name=Name(58:end);
figure
bar(AliveorNot)
set(gca,'XTickLabel',Case)
title([Name]);