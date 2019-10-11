#ifndef MYQUADMESH_MESH_H
#define MYQUADMESH_MESH_H

#include <string>
#include <dolfin.h>

using namespace dolfin;

void build(Mesh& mesh, int nx, int ny, double lx0, double lx1, double ly0, double ly1){
    MeshEditor editor;
    editor.open(mesh, CellType::quadrilateral, 2, 2);

    // Create vertices and cells:
    editor.init_vertices_global((nx + 1)*(ny + 1), (nx + 1)*(ny + 1));
    editor.init_cells_global(nx*ny, nx*ny);

    // Storage for vertices
    std::vector<double> x(2);

    double a = lx0, b = lx1;
    double c = ly0, d = ly1;

    // Create main vertices:
    std::size_t vertex = 0;
    for (std::size_t iy = 0; iy <= ny; iy++){
        x[1] = c + iy*(d - c)/ny;
        for (std::size_t ix = 0; ix <= nx; ix++){
            x[0] = a + ix*(b - a)/nx;
            editor.add_vertex(vertex, x);
            vertex++;
        }
    }

    // Create rectangles
    std::size_t cell = 0;
    std::vector<std::size_t> v(4);
    for (std::size_t iy = 0; iy < ny; iy++){
        for (std::size_t ix = 0; ix < nx; ix++){
            v[0] = iy*(nx + 1) + ix;
            v[1] = v[0] + 1;
            v[2] = v[0] + (nx + 1);
            v[3] = v[1] + (nx + 1);
            editor.add_cell(cell, v);
            ++cell;
        }
    }
    // Close mesh editor
    editor.close();
};
#endif
