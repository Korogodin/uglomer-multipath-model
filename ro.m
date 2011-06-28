function f=ro(x)
   TauChip = 300;
   f = (abs(x) < TauChip).*(1 - abs(x)/TauChip);
%    if (abs(x) < TauChip)
%        f=1-abs(x)/TauChip;
%    else
%        f=0;
%    end