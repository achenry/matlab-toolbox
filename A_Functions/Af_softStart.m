function [actual,generic] = Af_softStart(startTime,endTime,tt)

as = Af_sigma(startTime,endTime);

T = [tt.^3,tt.^2,tt,ones(length(tt),1)]';

actual  = zeros(size(tt));
generic = zeros(size(tt));

actual(tt<startTime)   = 0;
generic(tt<startTime)  = 1;

actual(tt>=startTime&tt<endTime)     = as * T(:,tt>=startTime&tt<endTime); 
generic(tt>=startTime&tt<endTime)    = 1-as * T(:,tt>=startTime&tt<endTime); 

actual(tt>=endTime)       = 1;
generic(tt>=endTime)      = 0;