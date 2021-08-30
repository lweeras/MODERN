clc
clear
%%Load Input Data
clusterData=load('ClusterData.txt'); % Input cluster data
Species=load('SpeciesData.txt'); % input speciecis avaialbility data
NumClu=4; % prvide number of clusters
NumSpe=3; % provide number of species in the system


%%%%%%% Do not change below %%%%%%%%%%%%%%%
[m,n]=size(clusterData);
Label=1;
k=1;
Data=[]; ExtraSpecies=[];
for i=1:m
    for j=1:n
        if (clusterData(i,j)==1)
            clusterData(i,j)=0;
            
        end
        
    end
end
clusterData;
for p=1:NumClu
    CurrentCluster=zeros(m,n);
    for i=1:m
        for j=1:n
            if (clusterData(i,j)==p+1)
                CurrentCluster(i,j)=1;
                
            end
            
        end
    end
    Data=[Data;CurrentCluster];
end
Data;
for p=1:NumClu
    CurrentClu=Data((p-1)*n+1:p*n,:);
    MultiData(:,:,p)=CurrentClu;
end
% for Species
for p=1:NumSpe
    CurrentSpe=Species((p-1)*n+1:p*n,:);
    AllSpe(:,:,p)=CurrentSpe;
end
[speciescount,Occurrence]=get_SpeciesCount(Species,NumClu,MultiData,NumSpe,m,n);
dlmwrite('C:\CompactConnect2020\speciescount.txt',speciescount,'delimiter','\t');
dlmwrite('C:\CompactConnect2020\Occurrence.txt',Occurrence,'delimiter','\t');
AllSpe;
MultiData;
while (NumClu~=1)
    Data;
    NumClu;
    for p=1:NumClu
        CurrentClu=Data((p-1)*n+1:p*n,:);
        MultiData(:,:,p)=CurrentClu;
    end
    OldData=MultiData;
    MultiData;
    Label=1;
    [AlterPath,AllList,LevelSize,Label]=To_JoinCluctersAlternatives(MultiData, m,n,Label,NumClu);
    AlterPath;
    %CL1cell meeting place of the first cluster with the cluster number and
    %level
    %CL2cell meeting place of the second cluster with the cluster number and
    %level
    %AllList Need tp trace back to get the paths
    % dlmwrite('C:\CompactConnect2020\AllList.txt',AllList,'delimiter','\t');
    MeetingLevel=Label;
    UnmodifiedMD=MultiData;
    AlterPath=Clean_AlterPath(AlterPath,MeetingLevel); % need this to remove duplicates
    AlterPath ;% All meeting places with the same lengths
    [NumPath,~]=size(AlterPath);
    counter=0; AllConnectors=[];IndCov=[];
    NumPath;
    while (NumPath>0)
        NumPath;
        CL1cell=AlterPath(2*counter+1,:);% Intersection of First Clsuter
        CL2cel2=AlterPath(2*counter+2,:);% Intersection of Second Clsuter
        clus1= CL1cell(3);% First Clsuter
        clus2= CL2cel2(3);% First Clsuter
        MeetingLevel=[CL1cell(4)+1,CL2cel2(4)+1];
        [FinalPath,IndividualCoverage]=get_Coridors(clus1,clus2,AllList,CL1cell,CL2cel2,LevelSize,NumClu,MeetingLevel,UnmodifiedMD,AllSpe,NumSpe,m,n) ;
        FinalPath;
        [numfinalPath,~]=size(FinalPath);
        AddClus1=clus1*ones(numfinalPath,1);
        AddClus2=clus2*ones(numfinalPath,1);
        counter=counter+1;
        NumPath=NumPath-2;
        AllConnectors=[AllConnectors;AddClus1,AddClus2,FinalPath];
        IndCov=[IndCov;IndividualCoverage];
        
    end
    AllConnectors;
    [numPath,~]=size(AllConnectors);
    [BestPath,AdditionalC]=Get_BestPath(AllConnectors);
    clus1=BestPath(1);
    clus2=BestPath(2);
    BestPath=BestPath(3:end);
    MultiData(:,:,clus1);
    MultiData(:,:,clus2);
    [MergeCluster,clusterData]=Get_CombinedSystem(BestPath,clus1,clus2,MultiData,m,n,NumClu,clusterData);
    MergeCluster;
    clusterData;

    ExtraSpecies=[ExtraSpecies,AdditionalC];
    UpdateMultiData=[MergeCluster]; % UpdateMultiData contains the Clus 1 and Clus 2 data
    for i=1:NumClu
        if(i ~=clus1 && i~=clus2)
            UpdateMultiData=[UpdateMultiData;OldData(:,:,i)];
        end
    end
    UpdateMultiData;
    IndCov;
    [numpath,~]=size(AllConnectors);
    dlmwrite('C:\CompactConnect2020\AllConnectors.txt',AllConnectors,'delimiter','\t');
    dlmwrite('C:\CompactConnect2020\Allindividual coverage.txt',IndCov,'delimiter','\t');
    dlmwrite('C:\CompactConnect2020\ClusterStep.txt',clusterData,'delimiter','\t');
    Data=UpdateMultiData;
    MultiData=[];
    Label=1;
    NumClu=NumClu-1;
end
AllConnectors;
ExtraSpecies;
MergeCluster;
[Newspeciescount,Occurrence]=get_SpeciesCount(Species,NumClu,MergeCluster,NumSpe,m,n)
Coridorcoverage=Newspeciescount-speciescount
SP=[Newspeciescount;speciescount;Coridorcoverage]
dlmwrite('C:\CompactConnect2020\SpeciesAnalysis.txt',SP,'delimiter','\t');
dlmwrite('C:\CompactConnect2020\ClustermultiData.txt',clusterData,'delimiter','\t');
dlmwrite('C:\CompactConnect2020\multiData.txt',MergeCluster,'delimiter','\t');
dlmwrite('C:\CompactConnect2020\AllConnectors.txt',AllConnectors,'delimiter','\t');