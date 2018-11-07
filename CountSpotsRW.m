nc_number=[CompiledParticles.nc];
CompiledParticles_14=CompiledParticles(nc_number==14);
nuclei=unique([CompiledParticles_14.Nucleus]);

for ii=1:length(nuclei)
    temp=find([CompiledParticles_14.Nucleus]==nuclei(ii));
    num_spots=length(temp);
    num_spots_ap(ii).spots=num_spots;
    num_spots_ap(ii).MeanAP=CompiledParticles_14(temp(1)).MeanAP;
end
APEstm=[num_spots_ap.MeanAP];
for j=1:length(nuclei)
    APEstm(j)=round(num_spots_ap(j).MeanAP,2);
    for jj=1:length(APbinID)
        if APEstm(j) < APbinID(jj)
            num_spots_ap(j).APBin=APbinID(jj);
            break;
        end
    end
end

APbins = unique([num_spots_ap.APBin]);   %need to add part to add up all the zero spot nuclei, also make sure only the double counted nuclei were changed to nan??
SpotsAP=zeros(length(APbins),3);
for ii=1:length(APbins)
    jj = find([num_spots_ap.APBin]==APbins(ii));
    
    Spotnumber=[];
   
    Spotnumber=[num_spots_ap(jj).spots];
    OneSpot=sum([Spotnumber]==1);
    TwoSpot=sum([Spotnumber]==2);
    NoSpot=sum([Spotnumber]==0);
    SpotsAP(ii,1)=OneSpot;
    SpotsAP(ii,2)=TwoSpot;
    SpotsAP(ii,3)=NoSpot;
    
end
figure 
bar(APbins,SpotsAP,'stacked')
xlabel('Egg Length'); ylabel('# of nuclei')
legend('one spot','two spots','Location','best');
title('ON nuclei nc14');


% 
%     
% %     SpotsAP(1,ii)=[sum([Spotnumber] ==1)];
% %     SpotsAP(2,ii)=[sum([Spotnumber] ==2)];
%     
% end
% 
%     Spotz=[num_spots_ap(jj).spots];
%     SpotsAP=[];
%     
%     SpotsAP(1,ii)=%# of nuclei with 1 spot at this APbin
%     SpotsAP(2,ii)= % # with 2 spots at this APbin
%     %count the number of 0/1/2 spots and store in data structure
%     %look up the right data structure for a bar graph
%     
% end
% bar()

