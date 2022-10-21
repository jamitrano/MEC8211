function [rightMember,stateVector] = AddNewmannBorderCondition(rightMember,stateVector,newmannCondition,ordre)
%Prise en compte des conditions de Newmann 

    for i= length(newmannCondition(:,1))
        idx = newmannCondition(i,1);
        if ordre==1
            rightMember(idx,1:2) = [1 -1];
        end
        if ordre==2
            rightMember(idx,1:3) = [-3 4 -1];
        end
        stateVector(idx) = 0;
   
    end
end

