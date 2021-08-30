clc
clear
%Data=load('data.txt')
FinalCluster=load('ClustermultiData.txt');
cost=load('cost.txt');
Species=load('SpeciesData.txt');
NumSpe=3; % Number of species
M=4; % Number of corridor sites
rhs=load('rhs.txt');
[m,n]=size(FinalCluster);
commonedges=0; BL=0;
FinalCluster=[zeros(1,n);FinalCluster;zeros(1,n)];
FinalCluster=[zeros(m+2,1),FinalCluster, zeros(m+2,1)];
[m,n]=size(FinalCluster)
Data=FinalCluster;
% end
BDCells=zeros(m,n);

for i=2:m-1
    for j=2:n-1
        if (Data(i,j)>=1)
            
            if (Data(i,j+1)==0)
                BL=BL+1;
                BDCells(i,j)=1;
           
            end
            
            if (Data(i,j-1)==0)
                BL=BL+1;
                BDCells(i,j)=1;
            end
           
        
            if (Data(i+1,j)==0)
                BL=BL+1;
                BDCells(i,j)=1;
       
            end
         
            if (Data(i-1,j)==0)
                BL=BL+1;
                BDCells(i,j)=1;
      
            end
        end
    end
    
  
end
BDCells; % potential BD

FeasibleBD=zeros(m,n);


for i=2:m-1
    for j=2:n-1
        if (BDCells(i,j)==1)
            BDCells(i,j)=0;
            FinalCluster(i,j)=0;
            %True= get_traverse(m,n, FinalCluster )
            True = get_AdjacencyMatrix (FinalCluster,m,n);
            if (True==0)
                FeasibleBD(i,j)=1;
            end
            FinalCluster(i,j)=1;
            BDCells(i,j)=1;
        end
    end
end
FeasibleBD; % Feasibel BD
FeasibleBD=FeasibleBD(2:m-1,2:n-1)
BDCells=BDCells(2:m-1,2:n-1)
m=m-2; n=n-2;
dlmwrite('BDCells.txt',BDCells,'delimiter','\t');
dlmwrite('FeasibleBD.txt',FeasibleBD,'delimiter','\t');
%[m,n]=size(BDCells)
for p=1:NumSpe
    CurrentSpe=Species((p-1)*n+1:p*n,:);
    AllSpe(:,:,p)=CurrentSpe;
end

for p=1:NumSpe
    Specificspecies=zeros(m,n);
    for i=1:m
        for j=1:n
            
            if (AllSpe(i,j,p)>0 && FeasibleBD(i,j)==1)
               % BDCellswithSpepcies(i,j)=BDCellswithSpepcies(i,j)+1;
                Specificspecies(i,j)=AllSpe(i,j,p);
            end

        end
    end
    IndividualBD(:,:,p)=Specificspecies;
end
IndividualBD;
%dlmwrite('C:\CompactConnect2020\BDCells.txt',BDCells,'delimiter','\t');

ConnectCost=51.36;
num_row=m;num_col=n;
numVar=num_row*num_col;% total number of variables
%rhs=[rhs,-M]; % upper bound for the number of cells 
rhs=[rhs,-ConnectCost]; % upper bound for the number of cells 
rhs;
UB=ones(1,numVar);
LB=zeros(1,numVar);
constraint=[];
for p=1:NumSpe
    species=IndividualBD(:,:,p);
    Spconst = get_const(num_row,num_col,species);
    constraint=[constraint;Spconst];

end
OldfeasibleBD=FeasibleBD;
for i=1:num_row
    for j=1:num_col
        FeasibleBD(i,j)=FeasibleBD(i,j)*cost(i,j);
    end
end
FeasibleBD;
FeasibleBDConst=get_const(num_row,num_col,FeasibleBD);
constraint=[constraint;-FeasibleBDConst];
FeasibleBD=OldfeasibleBD;
constraint;
[x, fval] = intlinprog(-FeasibleBDConst,[1:1:numVar], constraint, rhs,[],[], LB, UB,[]);% optimization
sum(x);

 CheckSol=[];
for l=1:num_row
    for m=1:num_col
        if (x((l-1)*num_col+m)==1)
            FeasibleBD(m,l)=2;
            coveredspe=[m,l];
         for p=1:NumSpe
             coveredspe=[coveredspe,IndividualBD(m,l,p)];
         end
         CheckSol=[CheckSol;coveredspe];
        end
    end
end
FeasibleBD;
CheckSol
dlmwrite('CheckSol.txt',CheckSol,'delimiter','\t');
dlmwrite('FeasibleBDSelected.txt',FeasibleBD,'delimiter','\t');


