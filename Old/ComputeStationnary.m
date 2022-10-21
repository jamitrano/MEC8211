function stationnary = ComputeStationnary(rightMemberMatrix,sourceTerm,stateVector,dirichletCondition,newmannCondition,ordre)
 [rightMember,stateVector] = AddNewmannBorderCondition(rightMemberMatrix,stateVector,newmannCondition,ordre);
 [stateVector,rightMember] = AddDirichletBorderCondition(stateVector,rightMember,dirichletCondition);
 stationnary = rightMember\sourceTerm + stateVector;
end

