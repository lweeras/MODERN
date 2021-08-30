clc
clear
FinalCluster=load('ClustermultiData.txt');
cost=load('cost.txt');
Species=load('SpeciesData.txt');
rhs=load('rhs.txt')
ConnectCost=51.36; % the additional cost (in excess of total budget available for the reserve system) 
NumSpe=3; % the number of species in the reserve system
[m,n]=size(FinalCluster);
commonedges=0; BL=0;
FinalCluster=[zeros(1,n);FinalCluster;zeros(1,n)];
FinalCluster=[zeros(m+2,1),FinalCluster, zeros(m+2,1)];
[m,n]=size(FinalCluster);
Data=FinalCluster;
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
FeasibleBD=FeasibleBD(2:m-1,2:n-1);
BDCells=BDCells(2:m-1,2:n-1);
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
M=4; % Number of corridor sites

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
        FeasibleBD(i,j)=FeasibleBD(i,j)*cost(i,j); % this is modfied const for model P5, if we use unform cost c(i,j)=1
    end
end
FeasibleBD;
FeasibleBDConst=get_const(num_row,num_col,FeasibleBD);
constraint=[constraint;-FeasibleBDConst];
ObjFun=get_const(num_row,num_col,OldfeasibleBD);% obj fun with nonunifromcost
FeasibleBD=OldfeasibleBD;
rhs;
constraint;
size(constraint);
%[x, fval] = intlinprog(-FeasibleBDConst,[1:1:numVar], constraint, rhs,[],[], LB, UB,[]);% optimization for the P_5 with uniform cost
[x, fval] = intlinprog(ObjFun,[1:1:numVar], constraint, rhs,[],[], LB, UB,[]);% optimization for non uniform cost
checkincost=0;
sloselect=[x';FeasibleBDConst];
for l=1:numVar
   
        checkincost=checkincost+x(l)*FeasibleBDConst(l);
   
end
checkincost;
selectedset=sum(x);
RemoveCost=0;
x';
 CheckSol=[];
for s=1:num_row
    for r=1:num_col
        if (x((s-1)*num_col+r)>0)
           [ s, r, (s-1)*num_col (s-1)*num_col+r cost(s,r)];
            FeasibleBD(r,s)=2;
            RemoveCost=RemoveCost+cost(r,s);
            coveredspe=[r,s];
         for p=1:NumSpe
             coveredspe=[coveredspe,IndividualBD(r,s,p)];
         end
         CheckSol=[CheckSol;coveredspe,cost(r,s)];
        end
    end
end
FeasibleBD;
CheckSol
RemoveCost ;% the cost of the removed sites
name = 'Alice';   
age = 12;
X = ['the cost of the removed sites is ',num2str(RemoveCost)];
disp(X)

dlmwrite('CheckSol.txt',CheckSol,'delimiter','\t');
dlmwrite('FeasibleBDSelected.txt',FeasibleBD,'delimiter','\t');
str = ["x-cordinate","x-cordinate","Spec 1",  "Spec 2","Spec 3","Cost"];
disp('The following table shows the removed sites along with the species coverage and the cost')
CheckSol=[str;CheckSol]

