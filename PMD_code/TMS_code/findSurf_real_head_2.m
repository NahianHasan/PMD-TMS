function [pp,Anor,tri]=findSurf_real_head_2(te2p,p,teid,scth,th_hair,EEG_data_points,EEG_channels)

    [pp,tri]=icosphere(4);
    pcen=(pp(tri(:,1),:)+pp(tri(:,2),:)+pp(tri(:,3),:))/3;
    th=acos(pcen(:,3));
    ph=atan2(pcen(:,2),pcen(:,1));
    tri=tri(th<pi/4,:);%what part of the spere to keep
    ntr=size(tri);
    [p1,~,tri]=unique(tri(:));
    pp=pp(p1,:);tri=reshape(tri,ntr);
    trisurf(tri,pp(:,1),pp(:,2),pp(:,3))
    th_pp=acos(pp(:,3));
    ph_pp=atan2(pp(:,2),pp(:,1));
    
    

    %%%%%%%%%%%%%%%%%% Construct EEG scatter interpolants%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F=fopen('EEG10-10_UI_Jurak_2007.csv','r');
    tline = fgetl(F);
    EEGPts1 = [];
    EEGLab1 = {};
    cnt=1;
    while ischar(tline)
        L=tline;
        t=strsplit(L,',');
        channel = t{end};
        EEGPts1 = [EEGPts1; str2double(t{2}), str2double(t{3}), str2double(t{4})];
        EEGLab1{cnt} = channel;
        cnt = cnt+1;
        tline = fgetl(F);
    end
    fclose(F);
    figure()
    trisphere=delaunay(EEGPts1(:,1),EEGPts1(:,2),EEGPts1(:,3));
    [trisphere]=surftri(EEGPts1,trisphere);
    trisurf(trisphere,EEGPts1(:,1),EEGPts1(:,2),EEGPts1(:,3))
    
    EEGPts1(:,1) = EEGPts1(:,1)./max(EEGPts1(:,1));
    EEGPts1(:,2) = EEGPts1(:,2)./max(EEGPts1(:,2));
    EEGPts1(:,3) = EEGPts1(:,3)./max(EEGPts1(:,3));
    
    
    th=acos(EEGPts1(:,3))/2;
    ph=atan2(EEGPts1(:,2),EEGPts1(:,1));
    trisphere=delaunay(EEGPts1(:,1),EEGPts1(:,2),EEGPts1(:,3));
    [trisphere]=surftri(EEGPts1,trisphere);

    eegcoo_x=scatteredInterpolant(th,ph,EEG_data_points(:,1),'natural','nearest');
    eegcoo_y=scatteredInterpolant(th,ph,EEG_data_points(:,2),'natural','nearest');
    eegcoo_z=scatteredInterpolant(th,ph,EEG_data_points(:,3),'natural','nearest');

    
    %EEG = [eegcoo_x(th(:),ph(:)),eegcoo_y(th(:),ph(:)),eegcoo_z(th(:),ph(:))];
    %trisurf(trisphere,EEG(:,1),EEG(:,2),EEG(:,3))
    
    
    v1=EEG_data_points(trisphere(:,2),:)-EEG_data_points(trisphere(:,1),:);
    v2=EEG_data_points(trisphere(:,3),:)-EEG_data_points(trisphere(:,1),:);
    nhat=cross(v1,v2,2);
    area=sqrt(nhat(:,1).^2+nhat(:,2).^2+nhat(:,3).^2);
    nhat(:,1)=nhat(:,1)./area;
    nhat(:,2)=nhat(:,2)./area;
    nhat(:,3)=nhat(:,3)./area;

    nhatp=zeros(size(EEG_data_points));
    for i=1:numel(nhat(:,1))
        nhatp(trisphere(i,1),:)=nhatp(trisphere(i,1),:)+nhat(i,:);
        nhatp(trisphere(i,2),:)=nhatp(trisphere(i,2),:)+nhat(i,:);
        nhatp(trisphere(i,3),:)=nhatp(trisphere(i,3),:)+nhat(i,:);
    end
    for i=1:numel(nhatp(:,1))
        nhatp(i,:)=nhatp(i,:)/norm(nhatp(i,:));
    end
    
    figure()
    trisurf(trisphere,EEG_data_points(:,1),EEG_data_points(:,2),EEG_data_points(:,3))
    hold on
    quiver3(EEG_data_points(:,1),EEG_data_points(:,2),EEG_data_points(:,3),nhatp(:,1),nhatp(:,2),nhatp(:,3))
    hold off
   
    for i=1:numel(nhat(:,1))
        nhat(i,:)=nhat(i,:)/norm(nhat(i,:));
    end
    eegcoo_nx=scatteredInterpolant(th,ph,nhatp(:,1),'linear','nearest');
    eegcoo_ny=scatteredInterpolant(th,ph,nhatp(:,2),'linear','nearest');
    eegcoo_nz=scatteredInterpolant(th,ph,nhatp(:,3),'linear','nearest');
    
    figure()
    trisurf(trisphere,EEG_data_points(:,1),EEG_data_points(:,2),EEG_data_points(:,3))
    hold on
    scatter3(EEG_data_points(:,1),EEG_data_points(:,2),EEG_data_points(:,3))
    quiver3(EEG_data_points(:,1),EEG_data_points(:,2),EEG_data_points(:,3),10*nhatp(:,1),10*nhatp(:,2),10*nhatp(:,3))
    axis equal
    hold off


    figure()
    pp = [eegcoo_x(th_pp(:),ph_pp(:)),eegcoo_y(th_pp(:),ph_pp(:)),eegcoo_z(th_pp(:),ph_pp(:))];
    trisphere_pp=delaunay(pp(:,1),pp(:,2),pp(:,3));
    [trisphere_pp]=surftri(pp,trisphere_pp);
    
    
    trisurf(trisphere_pp,pp(:,1),pp(:,2),pp(:,3))
    hold on
    nhatp = [eegcoo_nx(th_pp(:),ph_pp(:)),eegcoo_ny(th_pp(:),ph_pp(:)),eegcoo_nz(th_pp(:),ph_pp(:))];
    scatter3(pp(:,1),pp(:,2),pp(:,3))
    hold on
    pp=pp+nhatp*th_hair;
    quiver3(pp(:,1),pp(:,2),pp(:,3),10*nhatp(:,1),10*nhatp(:,2),10*nhatp(:,3))







    
    %extract point normals and local coordinate systems
    Anor=zeros([4 4 numel(pp)/3]);
    
    
    thatp1=zeros(size(nhatp));
    thatp1(:,1)=1-nhatp(:,1).*nhatp(:,1);
    thatp1(:,2)=-nhatp(:,1).*nhatp(:,2);
    thatp1(:,3)=-nhatp(:,1).*nhatp(:,3);
    for i=1:numel(thatp1(:,1))
        thatp1(i,:)=thatp1(i,:)/norm(thatp1(i,:));
    end
    thatp2=cross(nhatp,thatp1);
    pp=pp+nhatp*th_hair;
    np = numel(thatp1(:,1));
    Anor(1:3,1,:)=reshape(thatp1',[3 1 np]);
    Anor(1:3,2,:)=reshape(thatp2',[3 1 np]);
    Anor(1:3,3,:)=reshape(nhatp',[3 1 np]);
    Anor(1:3,4,:)=reshape(pp',[3 1 np]);
    Anor(4,4,:)=1;
end
