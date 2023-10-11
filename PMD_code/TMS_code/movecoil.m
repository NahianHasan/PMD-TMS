function [rs,js]=movecoil(rs2,js2,theta,nhat,pt,Ncoils)
tic


%ASSUME COIL IS ON X,Y PLANE
vx=[nhat(3),0,-nhat(1)]/norm([nhat(3),0,-nhat(1)]);
vy=[0,-nhat(3),nhat(2)]/norm([0,-nhat(3),nhat(2)]);
vhat1=vx*cos(theta)+vy*sin(theta);
vhat2=vx*sin(theta)-vy*cos(theta);

rs=zeros(size(rs2));js=zeros(size(js2));
rs(1,1:end)=rs2(1,:)*vhat1(1)+rs2(2,:)*vhat2(1)+rs2(3,:)*nhat(1)+pt(1);
rs(2,1:end)=rs2(1,:)*vhat1(2)+rs2(2,:)*vhat2(2)+rs2(3,:)*nhat(2)+pt(2);
rs(3,1:end)=rs2(1,:)*vhat1(3)+rs2(2,:)*vhat2(3)+rs2(3,:)*nhat(3)+pt(3);
for i=1:Ncoils
js(1,1:end,i)=js2(1,:,i)*vhat1(1)+js2(2,:,i)*vhat2(1)+js2(3,:,i)*nhat(1);
js(2,1:end,i)=js2(1,:,i)*vhat1(2)+js2(2,:,i)*vhat2(2)+js2(3,:,i)*nhat(2);
js(3,1:end,i)=js2(1,:,i)*vhat1(3)+js2(2,:,i)*vhat2(3)+js2(3,:,i)*nhat(3);
end
end