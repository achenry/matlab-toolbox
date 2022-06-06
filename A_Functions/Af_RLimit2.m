function [ Xrl ] = Af_RLimit2( X,RLx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Xdiff=diff(X);
Xrl=X;
ding=0;
for t=1:length(Xdiff)
    if abs(Xdiff(t))>RLx
        ding=sign(Xdiff(t));
    end
    if ding==1
        Xrl(t+1)=Xrl(t)+RLx;
        if Xrl(t+1)>X(t+1)
            Xrl(t+1)=X(t+1);
            ding=0;
        end
    elseif ding==-1
        Xrl(t+1)=Xrl(t)-RLx;
        if Xrl(t+1)<X(t+1)
            Xrl(t+1)=X(t+1);
            ding=0;
        end
    else
      Xrl(t+1)=X(t+1);
    end
end
% for t=2:length(X)
%     if  abs((Xrl(t)-Xrl(t-1)))>(RLx)
%         Xrl(t)=Xrl(t-1)+RLx*sign((Xrl(t)-Xrl(t-1)));
%     end
% end
% 
% end

