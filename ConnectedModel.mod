int n=...; // Number of nodels
int c=...; // Number of Clusters
int SmallRes=...;//
float Lambda1=...;
float Lambda2=...;
range node=1..n;
range parclel=1..SmallRes;
int List[node,node]=...;
int Cost1[node,node]=...;
dvar boolean X[node][node]; //  DECISION VARIABLES X;
dvar boolean Y[node][node][node];//  DECISION VARIABLES Y;



 minimize        sum(j in node, k in node: (List[j,k]==1 && j<=SmallRes))Cost1[j,k]*X[j][k]+sum(j in node, k in node: (List[j,k]==1 && j>SmallRes))1000*X[j][k];  
                    subject to {

 
forall (r in node,j in node :(r>SmallRes && r<c+SmallRes) && (j==r)  )  
sum ( k in node) List[k,j]*Y[r,j,k] -sum ( i in node) List[j,i]*Y[r,i,j]  ==1;

forall (r in node,j in node :(r>SmallRes && r<c+SmallRes) && (j==c+SmallRes)  )  
 sum ( k in node) List[k,j]*Y[r,j,k] -sum ( i in node) List[j,i]*Y[r,i,j]  ==-1;

forall (r in node,j in node :(r>SmallRes && r<c+SmallRes) && (j!=r) && (j!=c+SmallRes) )  
sum ( i in node) List[j,i]*Y[r,i,j] - sum ( k in node) List[j,k]*Y[r,j,k] ==0;
   
forall (r in node,j in node :(r>SmallRes && r<c+SmallRes))
forall (k in node:List[j][k]==1)    Y[r,j,k]-X[j][k]<=0;  



}

 float temp; // Used to check the runtime necessary to get a solution.
 
execute{
 
  obj_final = cplex.getObjValue();
  
  
 }
 
  main {
 
  var before = new Date();
  temp = before.getTime(); // Start timer.
  
  thisOplModel.generate();
  cplex.exportModel("MathModel.lp")
  thisOplModel.unconvertAllIntVars();
  cplex.solve();
  
  var after = new Date();
  writeln("Solving Time ~= ", after.getTime()-temp); // End timer and print.
  
  writeln("Integer Answers");  
  thisOplModel.postProcess();
  
  var f=new IloOplOutputFile("result.txt");
  
  f.writeln(thisOplModel.printSolution());
  
  f.close();
  

  
}              