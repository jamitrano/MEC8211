function [rightMember,stateVector] = AddNewmannBorderCondition(rightMember,stateVector,Dx,newmannCondition)
%Prise en compte des conditions de Newmann 

    for i= length(newmannCondition(:,1))
        idx = newmannCondition(i,1);
        rightMember(idx,:) = Dx(i,:);
        stateVector(idx) = 0;
   
    end
end

