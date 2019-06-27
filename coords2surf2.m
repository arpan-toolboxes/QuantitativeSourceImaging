function coords2surf2(roi_coords1,roi_coords2,roi_coords3, roi_coords4)

% Plots the significant voxels on the surface of the MNI template (taken
% from SPM. Colours are fixed as (legend):
% roi_coords1 = red
% roi_coords2 = green
% roi_coords3 = blue
% roi_coords1 and roi_coords2 = yellow
% roi_coords1 and roi_coords3 = pink
% roi_coords2 and roi_coords3 = cyan
% roi_coords1 and roi_coords2 and roi_coords3 = black
% Required 'Template.gii'
% % % 
% roi_coords1=[];roi_coords2=[];roi_coords3=[];
% 
% if (isempty(source1)~=1)
%     sourceP=source1.pow(source1.inside);
%     sourceC=source1.pos(source1.inside,:);
%     stat_roi= quantile(sourceP,threshold);
%     indi_stat=find(sourceP>stat_roi);
%     source_c1=sourceC(indi_stat,:);
%     source_p1=sourceP(indi_stat,:);
%     roi_coords1=source_c1;
% end
% 
% if (isempty(source2)~=1)
%     sourceP=source2.pow(source2.inside);
%     sourceC=source2.pos(source2.inside,:);
%     stat_roi= quantile(sourceP,threshold);
%     indi_stat=find(sourceP>stat_roi);
%     source_c2=sourceC(indi_stat,:);
%     source_p2=sourceP(indi_stat,:);
%     roi_coords2=source_c2;
% end
% 
% if (isempty(source3)~=1)
%     sourceP=source3.pow(source3.inside);
%     sourceC=source3.pos(source3.inside,:);
%     stat_roi= quantile(sourceP,threshold);
%     indi_stat=find(sourceP>stat_roi);
%     source_c3=sourceC(indi_stat,:);
%     source_p3=sourceP(indi_stat,:);
%     roi_coords3=source_c3;
% end

%%

pos=[];
dd=gifti('Template.gii');

distanceThreshold = 5;


faces = dd.faces; vertices1 = repmat(dd.vertices,2,1);
vertices=vertices1(1:length(faces),:);
facecolor1 = repmat(dd.cdata, 2,1);
facecolor=facecolor1(1:length(faces),:);
facecolorx=facecolor;
facecolorxx=facecolor;
facecolorxxx=facecolor;
% colors1 =repmat(colormap(autumn),size(pqr,1),1); for colour change
colors1 =repmat([1 0 0],size(roi_coords1,1),1);


for ii=1:size(roi_coords1, 1)
    pos1 = find(abs(vertices(1:length(vertices),1) - roi_coords1(ii, 1)) <= distanceThreshold & abs(vertices(1:length(vertices),2) - roi_coords1(ii, 2)) <= distanceThreshold & abs(vertices(1:length(vertices),3) - roi_coords1(ii, 3)) <= distanceThreshold  );
    facecolor(pos1,:) = repmat(colors1(ii,:), length(pos1), 1);
end

if (isempty(roi_coords2)~=1)
    xyz= roi_coords2;
    colors2 =repmat([0 1 0],size(xyz,1),1);
    
    for ii=1:size(xyz, 1)
        pos2 = find(abs(vertices(1:length(vertices),1) - roi_coords2(ii, 1)) <= distanceThreshold & abs(vertices(1:length(vertices),2) - roi_coords2(ii, 2)) <= distanceThreshold & abs(vertices(1:length(vertices),3) - roi_coords2(ii, 3)) <= distanceThreshold  );
        facecolorx(pos2,:) = repmat(colors2(ii,:), length(pos2), 1);
    end
    
end

if (isempty(roi_coords3)~=1)
    abc= roi_coords3;
    colors3 =repmat([0 0 1],size(abc,1),1);
    
    for ii=1:size(abc, 1)
        pos3 = find(abs(vertices(1:length(vertices),1) - roi_coords3(ii, 1)) <= distanceThreshold & abs(vertices(1:length(vertices),2) - roi_coords3(ii, 2)) <= distanceThreshold & abs(vertices(1:length(vertices),3) - roi_coords3(ii, 3)) <= distanceThreshold  );
        facecolorxx(pos3,:) = repmat(colors3(ii,:), length(pos3), 1);
    end
    
end

if (isempty(roi_coords4)~=1)
    abc= roi_coords4;
    colors4 =repmat([1 0 1],size(abc,1),1);
    
    for ii=1:size(abc, 1)
        pos4 = find(abs(vertices(1:length(vertices),1) - roi_coords4(ii, 1)) <= distanceThreshold & abs(vertices(1:length(vertices),2) - roi_coords4(ii, 2)) <= distanceThreshold & abs(vertices(1:length(vertices),3) - roi_coords4(ii, 3)) <= distanceThreshold  );
        facecolorxxx(pos4,:) = repmat(colors4(ii,:), length(pos4), 1);
    end
    
end
% RED GREEN BLUE
%%
count=1;
f=repmat([0 0 0],size(facecolor,1),1);
for i=1:size(facecolor)
    if facecolor(i,:)==[1 0 0];
        f(i,:)=[1 0 0];
    end
    if facecolorx(i,:)==[0 1 0]
        f(i,:)=[0 1 0];
    end
    if facecolorxx(i,:)==[0 0 1]
        f(i,:)=[0 0 1];
    end
    if facecolorxxx(i,:)==[1 0 1]
        f(i,:)=[1 0 1];
    end
    if facecolor(i,:)==[1 0 0] & facecolorx(i,:)==[0 1 0]
        f(i,:)=[1 1 0];
        count=count+1;
    end 
end
%%
for i=1:size(f,1)
    if f(i,:)==[0 0 0]
        f(i,:)=facecolor(i,:);
    end
    if f(i,:)==[1 1 1]
        f(i,:)=[0 0 0];
        pos(count,:)=f(i,:);
    end
end


%Audio: 1 165/255 0
%AV1: 1 0 0
%Visual: 0 165/255 1
%AV2: 0 0 1
%CrossAudio P2: [0 240/255 0]
%CrossAudio P3: [0 120/255 0]
%CrossVisual: [238/255 155/255 238/255]

% count=1;
% f=repmat([0 0 0],size(facecolor,1),1);
% for i=1:size(facecolor)
%     if facecolor(i,1)==1;
%         f(i,:)=[0 165/255 1];
%     end
%     if facecolorx(i,2)==1
%         f(i,:)=[1 0 1];
%     end
%     if facecolorx(i,2)==1 && facecolor(i,1)==1;
%         f(i,:)=[0 0 50/255];
%     end
% %     if facecolorxx(i,3)==1
% %         f(i,3)=1;
% %     end
% end
% 
% for i=1:size(f,1)
%     if f(i,:)==[0 0 0]
%         f(i,:)=facecolor(i,:);
%     end
%     if f(i,:)==[1 1 1]
%         f(i,:)=[0 0 0];
%         pos(count,:)=f(i,:);
%     end
% end
%%
p = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',f,...
    'FaceColor','flat','CDataMapping','direct','EdgeColor','none','facealpha',1);
set(p,'AmbientStrength',1.0,'DiffuseStrength',.0001);

daspect([1 1 1])
axis tight
camlight
lighting phong
axis off;
shading interp;