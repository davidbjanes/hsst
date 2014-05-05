function drawPCs(s,z);
N=numel(z);
nz=sparse(z,ones(N,1),ones(N,1));
[nz2,rendx]=sort(nz,'descend');
%
thres=10;
numClus=sum(nz2>thres);
colors=hsv(numClus);
%
clf;
hold on
for clus=1:numClus
    ndx=z==rendx(clus);
    plot(s(ndx,1),s(ndx,2),'.','color',colors(clus,:),'markersize',15);
end
hold off
legend(num2str(rendx(1:numClus)));
xlabel('PC-1')
ylabel('PC-2')

