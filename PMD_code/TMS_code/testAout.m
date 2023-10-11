clear all
rs(1:3,1)=[0 0 1] 
js(1:3,1)=[0 0 1] 

[X Y]=ndgrid(-1:.01:1,-1:.01:1);
robs(1,:)=X(:)';
robs(2,:)=Y(:)';
robs(3,:)=0;
[Aout,curlAout]=computeAprimary(rs,js,1,robs,numel(robs)/3);

A=10^-7./sqrt(X.^2+Y.^2+1);
curlAx=-10^-7*Y./sqrt(X.^2+Y.^2+1).^3;
curlAy=10^-7*X./sqrt(X.^2+Y.^2+1).^3;
subplot(2,3,1),imagesc(reshape(Aout(1,:),size(X))-curlAx);
subplot(2,3,2),imagesc(reshape(Aout(2,:),size(X))-curlAy);
subplot(2,3,3),imagesc(reshape(Aout(3,:),size(X)));

subplot(2,3,4),imagesc(reshape(curlAout(1,:),size(X)));
subplot(2,3,5),imagesc(reshape(curlAout(2,:),size(X)));
subplot(2,3,6),imagesc(reshape(curlAout(3,:),size(X)));


