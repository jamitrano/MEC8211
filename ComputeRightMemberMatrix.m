function [rightMemberMatrix,rightMemberMatrixStationnary] = ComputeRightMemberMatrix(diameter,N,Deff,reactionConstant,dt,ordre)
dr = diameter/2/(N-1);
dr2 = dr*dr;
rInv = 1./(linspace(0,N,N)*diameter/(2*N));
rInv(1)=0;
if ordre ==1
    backwardTerm = - Deff/(dr2).*ones(1,N);
    centralTerm = 2*Deff/(dr2) + Deff/dr.*rInv + 1/dt + reactionConstant;
    forwardTerm = -Deff/(dr2) -Deff/dr.*rInv;

    
end

if ordre==2
     backwardTerm = - Deff/(dr2).*ones(1,N) + Deff/(2*dr).*rInv;
    centralTerm = (2*Deff/(dr2)  + 1/dt + reactionConstant).*ones(1,N);
    forwardTerm = -Deff/(dr2) -Deff/(2*dr).*rInv;
end
rightMemberMatrix = diag(backwardTerm(1:N-1),-1) + diag(centralTerm,0) + diag(forwardTerm(2:N),1);

 centralTermStationnary = centralTerm - 1/dt ;
rightMemberMatrixStationnary = diag(backwardTerm(2:N),-1) + diag(centralTermStationnary,0) + diag(forwardTerm(1:N-1),1);
end

