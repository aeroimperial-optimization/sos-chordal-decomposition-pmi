function [R,MC2,SMC2] = chordalGen(n,thresh)
%
% From Richard Manson
%
%create a chordal graph with largest maximal clique size < thresh

%first it creates a connected tree and then merges adjacent maximal cliques
%together to produce denser chordal graphs.

A = diag(ones(n,1));
inTree = false(n,1);
inTree(1) = true;
v = 1:n;

for i=1:n-1
    u=v(~inTree);
    w = v(inTree);
    r = randi(n-i,1);
    s = randi(i,1);
    A(u(r),w(s))=1;
    A(w(s),u(r))=1;
    inTree(u(r))=true;
end

%L2 = A - diag(sum(A));

MC = maximalCliques(A);
SMC = sum(MC);
[p,q] = size(MC);

k=1;
while max(SMC)<thresh
    r = randi(q);
    w =MC(:,r)'*MC(:,1:q);
    w(r)=0;
   
    t = find(w);
    s = randi(length(t));
    v = find(MC(:,t(s)));
    z = find(MC(:,r));
 
    keep = min([r,t(s)]);
    del = max([r,t(s)]);
    MC([v;z],keep)=1;
    SMC(1,keep) = sum(MC(:,keep));
    
    if del<q
    MC(:,1:q-1) = MC(:,[1:del-1,del+1:q]);
    SMC(:,1:q-1) = SMC(:,[1:del-1,del+1:q]);
    q= q-1;
    else
    MC = MC(:,1:q-1);
    SMC = SMC(:,1:q-1);
    q = q-1;
    end
end

R = zeros(n,n);

for i=1:q
    indx = MC(:,i)>0;
    R(indx,indx) =1;
end

MC2 =MC(:,1:q);
SMC2= SMC(1:q);