% Usage is parsave(filename, number of variables, variable name with no quotations)

function parsave(fname,numvars,varargin)
for i = 1:numvars
   eval([inputname(i+2),'= varargin{i};']);  
end
save('-mat',fname,inputname(3),'-v7.3');
for i = 2:numvars    
    save('-mat',fname,inputname(i+2),'-append','-v7.3');
end