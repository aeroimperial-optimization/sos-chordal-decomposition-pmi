function combineData
% Combine data for individual runs in single file

fnames = dir('ex5_2_sparse_graph*.mat');

% Initalize containers
Gsize = length(fnames);
gSparse    = cell(5,length(Gsize));
exponents    = cell(5,length(Gsize));
objSparse  = zeros(5,length(Gsize));
timeSparse = zeros(5,length(Gsize));
% gStandard    = cell(3,length(Gsize));
% objStandard = zeros(3,length(Gsize));
% timeStandard = zeros(3,length(Gsize));

for i = 1:Gsize
   data = load([fnames(i).folder,filesep,fnames(i).name]) ;
   gSparse(:,i) = data.gSparse(:,i);
   objSparse(:,i) = data.objSparse(:,i);
   timeSparse(:,i) = data.timeSparse(:,i);
   exponents(:,i) = data.exponents(:,i); 
end

% Extract other stuff
Gsize = data.Gsize;
A = data.A;
B = data.B;

% Save
save('ex5_2_sparse.mat','gSparse','objSparse','timeSparse','Gsize','A','B','exponents');

end