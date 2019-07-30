% Write a .vtk file, readable by paraview, containing the
% mesh and the solution 
%
% u should be a vector of size ncell+nedge, with the cell values first
%
function out=write_solution_vtk(u,namefile,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex)

% We first compute the value of the solution at the vertices.
% The value at a vertex s is the average of the values in the cells K around s,
% ponderated by a weitgh equal to the angle of the cell K at s.

uvert=zeros(nvert,1);
sumcoefvert=zeros(nvert,1);

for i=1:ncell
	for j=1:size(cell_v{i},2)-1
		% The vertices around cell_v{i}(j)=j1 are j0 and j2 computed as follows
		j1=cell_v{i}(j);
		j2=cell_v{i}(j+1);
		if (j==1)
			j0=cell_v{i}(size(cell_v{i},2)-1);
		else
			j0=cell_v{i}(j-1);
		end
		% The angle of the cell i around j1 is arccos of the normalised scalar product
		% between j0j1 and j1j2
		vec2=vertex(j2,:)-vertex(j1,:);
		vec1=vertex(j0,:)-vertex(j1,:);
		scal=dot(vec1,vec2);
		% Cross-product between the vectors
		crossprod=vec2(1)*vec1(2)-vec2(2)*vec1(1);
		% If the cross product is negative, the angle is above pi and thus equal
		%	to 2pi-arcos(normalised scal). Otherwise, it's just arccos
		if (crossprod<0)
			coef=2*pi-acos(scal/(norm(vec1)*norm(vec2)));
		else
			coef=acos(scal/(norm(vec1)*norm(vec2)));
		end
%coef

		sumcoefvert(j1)=sumcoefvert(j1)+coef;
		uvert(j1)=uvert(j1)+coef*u(i);
	end
end

% Check if sumcoefvert is reasonable
if (min(sumcoefvert)<1e-3 | max(sumcoefvert)>2*pi+1e-4)
	disp('wrong sumcoefvert');
end
%min(sumcoefvert)
%max(sumcoefvert)
% Then we apply the total ponderation
uvert=uvert./sumcoefvert;

% We use edge boundary values to re-compute the values at boundary vertices,
% to get a better representation
% Start by re-putting all boundary values at 0
for i=1:ncell
	for j=1:size(cell_v{i},2)-1
		if (cell_n{i}(j)<=0)
			j1=cell_v{i}(j);
			j2=cell_v{i}(j+1);
			uvert(j1)=0;
			uvert(j2)=0;
		end
	end
end
% Here it's a brute average, we can do weighted later
for i=1:ncell
	for j=1:size(cell_e{i},2)
		% Boundary edge only
		if (cell_n{i}(j)<= 0)
			% vertices of the edge
			j1=cell_v{i}(j);
			j2=cell_v{i}(j+1);
			uvert(j1)=uvert(j1)+0.5*u(ncell+cell_e{i}(j));
			uvert(j2)=uvert(j2)+0.5*u(ncell+cell_e{i}(j));
		end
	end
end

%
% We write the vtk file
%

completename=strcat(namefile,'.vtk');
fid=fopen(completename,'w');


fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Grille\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET POLYDATA\n');
fprintf(fid,'POINTS %d float\n',nvert);
for i=1:nvert
	fprintf(fid,'%E %E %E\n',vertex(i,1),vertex(i,2),uvert(i));
end
sumvert=0;
for i=1:ncell
	sumvert=sumvert+size(cell_v{i},2)-1;
end
fprintf(fid,'POLYGONS %d %d\n',ncell,sumvert+ncell);
for i=1:ncell
% Should work for all type of mesh
		fprintf(fid,'%d\n',size(cell_v{i},2)-1,cell_v{i}(1:size(cell_v{i},2)-1)-1);
end
fprintf(fid,'CELL_DATA %d\n',ncell);
fprintf(fid,'SCALARS sol float 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
for i=1:ncell
	fprintf(fid,'%E\n',u(i));
end


fclose(fid);




