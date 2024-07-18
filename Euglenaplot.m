function [S, E]= Euglenaplot(xc,yc,zc,length,thickness,bodyrot,bodyangle,eyerot,eyeangle,transparent,ratio)

[x, y, z] = ellipsoid(0,0,0,length/2,thickness/2,thickness/2,30);
t = hgtransform;
S = surfl(x, y, z);
set(S,'parent',t);
% rotate(S,bodyrot,bodyangle)
R=makehgtform('axisrotate',bodyrot,bodyangle);
Tx = makehgtform('translate',[xc yc zc]);
t.Matrix = Tx*R;

set(S,'FaceColor',[0 1 0],'FaceAlpha',transparent,'EdgeColor','none'); 


% rotate(S,[1 0 0],bodyrotx)
% rotate(S,[0 1 0],bodyroty)
% rotate(S,[0 0 1],bodyrotz)
% hold on;
% [x2, y2, z2] = ellipsoid(0,0,0,0.03,0.03,0.03,10);
% E = surfl(x2, y2, z2);
% t2 = hgtransform;
% set(E,'parent',t2);
% R2=makehgtform('axisrotate',eyerot,eyeangle);
% Tx2 = makehgtform('translate',[xc+0.3 yc zc+0.06]);
% t2.Matrix = Tx2*R2;
% 
% set(E,'FaceColor',[1 0 0],'FaceAlpha',transparent*ratio,'EdgeColor','none'); 
% rotate(E,eyerot,eyeangle)
% rotate(E,[1 0 0],eyerotx)
% rotate(E,[0 1 0],eyeroty)
% rotate(E,[0 0 1],eyerotz)