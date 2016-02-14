function show(pic_no,triangles,coordinates,u, viewpoint)
figure(pic_no)
trisurf(triangles',coordinates(1,:),coordinates(2,:),full(u)','facecolor','interp','FaceLighting','phong');
grid('on');
light
lighting phong
%  colorbar
view(viewpoint)
%   zlim([0,2.3]);
