function drawWaveforms(x,z);

errorbarsint=4;
p=size(x,1);
time=(1:p);
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
    plot(time,mean(x(:,ndx),2),'-','color',colors(clus,:),'linewidth',2);
end
legend(num2str(rendx(1:numClus)));
for clus=1:numClus
    ndx=z==rendx(clus);
    rndx=clus:errorbarsint:p;
    t=time(rndx);
    m=mean(x(rndx,ndx),2);
    s=std(x(rndx,ndx)');
    errorbar(t,m,s,'.','color',colors(clus,:),'linewidth',2);
end
hold off
