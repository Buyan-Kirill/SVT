#include "inmost.h"
#include <stdio.h>
#include <string>


using namespace INMOST;
using namespace std;

void make_cells_count_tag(Mesh* m)
{
	Tag tagCells = m->CreateTag("Adjacent_cells", DATA_INTEGER, NODE, NONE, 1);
	for (Mesh::iteratorNode inode = m->BeginNode(); inode != m->EndNode(); inode++) {
		Node node = inode->getAsNode();
		node.Integer(tagCells) = node.getCells().size();
	}
}

void make_vector_to_center_tag(Mesh* m)
{
	Tag tagVecToCnt = m->CreateTag("Vector_to_center", DATA_REAL, CELL, NONE, 3);
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		Storage::real cnt[3];
		Storage::real mesh_cnt[3] = {37.64409676420321, 23.831310721372773, 0.};
		cell->Barycenter(cnt);
		cell.RealArray(tagVecToCnt)[0] = mesh_cnt[0] - cnt[0];
		cell.RealArray(tagVecToCnt)[1] = mesh_cnt[1] - cnt[1];
		cell.RealArray(tagVecToCnt)[2] = mesh_cnt[2] - cnt[2];
	}
}


Storage::real find_diameter(Mesh* m)
{
	Storage::real d = 0.;
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		auto nodes = cell.getNodes();
		for (int i = 0; i < 3; i++) {
			auto node1 = nodes[i].Coords();
			auto node2 = nodes[(i + 1) % 3].Coords();
			d = max(d, sqrt((node1[0] - node2[0]) * (node1[0] - node2[0]) + (node1[1] - node2[1]) * (node1[1] - node2[1])));
		}
	}
	return d;
}


int main(int argc, char *argv[])
{
	string mesh_file;
	if (argc != 2) {
		cout << "Usage: mesh <mesh.vtk>" << endl;
		mesh_file = "../mesh.vtk";
	} else {
		mesh_file = argv[1];
	}
    Mesh *m = new Mesh;
	m->Load(mesh_file);
	make_cells_count_tag(m);
	make_vector_to_center_tag(m);
	cout << "Diameter = " << find_diameter(m) << endl;
	m->Save("../res.vtk");
	delete m;
	cout << "Success!" << endl;
	return 0;
}
