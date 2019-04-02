#ifndef TRI_MESH_H
#define TRI_MESH_H


#include <surface_mesh/Surface_mesh.h>
#include <Eigen/Dense>
#include <fstream>
#include <stack>
// @brief given V and F, build the surface_mesh
// @param V: the vertex list of the mesh
// @param F: the face list of the mesh
// @param mesh: the constructed mesh
// @return void


static void build_mesh(const Eigen::Matrix3Xd &V,
                       const Eigen::Matrix3Xi &F,
                       surface_mesh::Surface_mesh &mesh){
    mesh.reserve(V.cols(), F.cols() * 2, F.cols());
    for(int i = 0;i < V.cols();++i){
        mesh.add_vertex(surface_mesh::Point(
                            V(0, i), V(1, i), V(2, i)));
    }
    for(int i = 0;i  < F.cols();++i){
		mesh.add_triangle(surface_mesh::Surface_mesh::Vertex(F(0, i)),
			surface_mesh::Surface_mesh::Vertex(F(1, i)),
			surface_mesh::Surface_mesh::Vertex(F(2, i)));
    }
}

static void find_halfedge(const surface_mesh::Surface_mesh &mesh,
                          const surface_mesh::Surface_mesh::Vertex &v1,
                          const surface_mesh::Surface_mesh::Vertex &v2,
                          surface_mesh::Surface_mesh::Halfedge &he)
{
    surface_mesh::Surface_mesh::Halfedge_around_vertex_circulator  vh_it, vh_end;
    vh_it = vh_end = mesh.halfedges(v1);
    do{
        const surface_mesh::Surface_mesh::Vertex vhj = mesh.from_vertex(*vh_it);
        if(vhj == v2){
            he = *vh_it;
            break;
        }
    }while (++vh_it != vh_end);
}

class tri_mesh_data : public surface_mesh::Surface_mesh
{
public:

    void build_topology(){
        build_mesh(nods_, tris_, *this);
    }

	int read_obj(const char *filename){
		if (this->read(filename)){
			nods_.resize(3, this->n_vertices());
			for (auto vit : this->vertices()){
				const surface_mesh::Point &p = this->position(vit);
				nods_.col(vit.idx()) << p[0], p[1], p[2];
			}
			tris_.resize(3, this->n_faces());
			for (auto fit : this->faces()){
				Surface_mesh::Vertex_around_face_circulator fvit = this->vertices(fit), fvend = fvit;
				size_t i = 0;
				do{
					tris_(i++, fit.idx()) = (*fvit).idx();
				} while (++fvit != fvend);
            }

			return 1;
		}
    }

    double calc_avg_edge_length(const Eigen::Matrix3Xd &V) const{
        double avg_edge_length_ = 0;
        for (auto eit : this->edges()){
            const Surface_mesh::Halfedge he = this->halfedge(eit, 0);
            avg_edge_length_ += (V.col(this->from_vertex(he).idx()) - V.col(this->to_vertex(he).idx())).norm();
        }
        avg_edge_length_ /= this->n_edges();
    }
	Eigen::Matrix3Xd nods_;
	Eigen::Matrix3Xi tris_;
    Eigen::Matrix2Xd tnods_;
    Eigen::Matrix3Xi ttris_;
};

static Eigen::Vector3d cal_vertex_normal(const tri_mesh_data &m,
                                         const Eigen::Matrix3Xd &V,size_t vid)
{
    Eigen::Vector3d n(0, 0, 0);
    surface_mesh::Surface_mesh::Face_around_vertex_circulator vfit, vfend;
    vfit = vfend = m.faces(surface_mesh::Surface_mesh::Vertex(vid));
    do{
        const int idx = (*vfit).idx();
        const Eigen::Vector3d v1 = V.col(m.tris_(1, idx)) - V.col(m.tris_(0, idx));
        const Eigen::Vector3d v2 = V.col(m.tris_(2, idx)) - V.col(m.tris_(0, idx));
        n += v1.cross(v2);
    } while (++vfit != vfend);
    return n.normalized();
}


static double calc_avg_edge_length(const Eigen::Matrix3Xd &V,
                                   const Eigen::Matrix3Xi &F)
{
    double length = 0;
    for(int i = 0;i  < F.cols();++i){
        for(size_t j = 0;j  < 3;++j){
            length += (V.col(F(j, i)) - V.col(F((j + 1)%3, i))).norm();
        }
    }
    return length / F.size();
}

static int label_vertices_by_connect_component(const tri_mesh_data &m,
                                               const Eigen::VectorXi &vflags,
                                               size_t seed,
                                               Eigen::VectorXi &flags)
{
    std::stack<int>  S;
    S.push(seed);
    flags.setOnes(vflags.size());
    flags[seed] = 0;
    const int flag = vflags[seed];
    while(!S.empty()){
        const int current = S.top();
        S.pop();
        const surface_mesh::Surface_mesh::Vertex v(current);
        for (auto vh : m.halfedges(v)){
            const surface_mesh::Surface_mesh::Vertex &vj = m.to_vertex(vh);
            const int idx = vj.idx();
            if(vflags[idx] == flag && flags[idx]){
                flags[idx] = 0;
                S.push(idx);
            }
        }
    }
    return 1;
}


#endif
