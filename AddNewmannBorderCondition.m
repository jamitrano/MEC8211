function [rightMember,stateVector] = AddNewmannBorderCondition(rightMember,stateVector,Dx,dx,newmannCondition,ordre)
%Prise en compte des conditions de Newmann 

    for i= length(newmannCondition(:,1))
        idx = newmannCondition(i,1);
        if ordre==1
            rightMember(idx,:) = Dx(i,:);
        end
        if ordre==2
            rightMember(idx,1:3) = 1/(dx*dx).*[1 -2 1];
        end
        stateVector(idx) = 0;
   
    end
end

