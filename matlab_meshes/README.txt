
% Mesh data:
%		ncell = number of cells
%		nedge = number of edges
%		nvert = number of vertices
%		cell_v	=	cell array of vertices
%				cell_v{i} = array [g,a,b,c,...,g] with the numbers g,a,b,c... of the cell i vertices, ordered counter-clockwise
%							Note that the last one (g) is repeated at the start, see below why
%		cell_e	=	cell array of edges
%				cell_e{i} = array [u,v,w,...] with the numbers u,v,w of the cell i edges, ordered counter-clockwise
%							The vertices of cell_e{i}(j) are given cell_v{i}(j) and cell_v{i}(j+1) (this is why we repeat
%							the vertex g above, so that this formula is always true. Note that the length of 
%							cell_v{i} is one more than the length of cell_e{i}).
%				Due to the counter-clockwise numbering of the vertices in cell_v{i}, rotating the vector
%					cell_v{i}(j),cell_v{i}(j+1) by -pi/2 gives the outer normal to the cell i on the
%					edge cell_e{i}(j)
%		cell_n =	cell array of neighbours
%				cell_n{i} = array [A,B,C...] with the numbers of the cells on the other side of the edges [u,v,...].
%											number=0 if boundary edge
%		center = (ncell x 2) array with the centers coordinates (usually gravity centers)
%				center(i,:) = coordinates of the center of cell i
%		area = array of areas of cells (area(i)=area of cell i)
%		diam = array of diameters of cells (diam(i)=diameter of cell i)
%				diam(i,:) = diameter of cell i
%		vertex = (nvert x 2) array of the vertices coordinates
%				vertex(i,:) = coordinates of vertex i

