#ifndef MESH_IO_H_H_
#define MESH_IO_H_H_

#include <string>
#include <Eigen/Dense>

namespace  common {
// @brief read and write obj mesh
// @param filename: the file name of the mesh
// @param V: the vertex list of the mesh
// @param F: the face list of the mesh
// @param tV: the texture vertex list of the mesh
// @param tF: the texture face list of the mesh
// @return 0: read or write failed; 1: read or write successufully

int read_obj(const std::string &filename, 
            Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F) 
{
    std::ifstream is(filename.c_str());
    if(!is)
        return 0;
    std::vector<double>  vs;
    std::vector<int>  fs;
    std::string line, pair[3];
    double  node[3];
    int  tri;
    while (!is.eof()) {
        std::getline(is, line);
        if (line.empty() || 13 == line[0])
            continue;
        std::istringstream instream(line);
        std::string word;
        instream >> word;
        if ("v" == word || "V" == word) {
            instream >> node[0] >> node[1] >> node[2];
            for (size_t j = 0; j < 3; ++j){
                vs.push_back(node[j]);
            }
        }
        else if ('f' == word[0] || 'F' == word[0]) {
            instream >> pair[0] >> pair[1] >> pair[2];
            for (size_t j = 0; j < 3; ++j) {
                tri = strtoul(pair[j].c_str(),NULL,10) - 1;
                fs.push_back(tri);
            }
        }
    }
    is.close();
    V.resize(3, vs.size() / 3);
    for (size_t i = 0, k = 0; i < V.cols(); ++i){
        for (size_t j = 0; j < 3; ++j){
            V(j, i) = vs[k++];
        }
    }
    F.resize(3, fs.size() / 3);
    for (size_t i = 0, k = 0; i < F.cols(); ++i){
        for (size_t j = 0; j < 3; ++j){
            F(j, i) = fs[k++];
        }
    }
    return 1;
}
int save_obj(const std::string &filename, 
            const Eigen::Matrix3Xd &nods,
            const Eigen::Matrix3Xi &tris)
{
    std::ofstream os(filename.c_str());
    if(!os)
        return 0;
    for (size_t i = 0; i < nods.cols(); ++i){
        os << "v " <<  nods.col(i).transpose() << "\n";
    }

    for (size_t i = 0; i < tris.cols(); ++i){
        const Eigen::Vector3i f = tris.col(i) + Eigen::Vector3i::Ones();
        os << "f " << f.transpose() << "\n";
    }
    os.close();
    return 1;
}

int read_obj(const std::string &in_file,
             Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F,
             Eigen::Matrix2Xd &tV, Eigen::Matrix3Xi &tF)
{
    std::ifstream is(in_file);
    if(!is)
        return 0;
    std::vector<double>  vs, vt;
    std::vector<int>  fs, ft;
    std::string line, pair[3];
    double  node[3];
    while (!is.eof()) {
        std::getline(is, line);
        if (line.empty() || 13 == line[0])
            continue;
        std::istringstream instream(line);
        std::string word;
        instream >> word;
        if ("v" == word || "V" == word) {
            instream >> node[0] >> node[1] >> node[2];
            for (size_t j = 0; j < 3; ++j){
                vs.push_back(node[j]);
            }
        }else if ("vt" == word || "VT" == word) {
            instream >> node[0] >> node[1];
            for (size_t j = 0; j < 2; ++j){
                vt.push_back(node[j]);
            }
        }
        else if ('f' == word[0] || 'F' == word[0]) {
            instream >> pair[0] >> pair[1] >> pair[2];
            //get vertex id in the this triangle face
            size_t tri;
            for (size_t j = 0; j < 3; ++j) {
                std::string vertex_str = pair[j].substr(0, pair[j].find('/'));
                sscanf_s(vertex_str.c_str(), "%lu", &tri);
                fs.push_back(--tri);

                vertex_str = pair[j].substr(pair[j].find('/')+1);
                sscanf_s(vertex_str.c_str(), "%lu", &tri);
                ft.push_back(--tri);
            }
        }
    }
    V.resize(3, vs.size() / 3);
    for (size_t i = 0, k = 0; i < V.cols(); ++i){
        for (size_t j = 0; j < 3; ++j){
            V(j, i) = vs[k++];
        }
    }
    F.resize(3, fs.size() / 3);
    for (size_t i = 0, k = 0; i < F.cols(); ++i){
        for (size_t j = 0; j < 3; ++j){
            F(j, i) = fs[k++];
        }
    }
    tV.resize(2, vt.size() / 2);
    for (size_t i = 0, k = 0; i < tV.cols(); ++i){
        for (size_t j = 0; j < 2; ++j){
            tV(j, i) = vt[k++];
        }
    }
    tF.resize(3, ft.size() / 3);
    for (size_t i = 0, k = 0; i < tF.cols(); ++i){
        for (size_t j = 0; j < 3; ++j){
            tF(j, i) = ft[k++];
        }
    }
    return 1;
}
int save_obj(const std::string &out_file,
             const Eigen::Matrix3Xd &V,const Eigen::Matrix3Xi &F,
             const Eigen::Matrix2Xd &tv,const Eigen::Matrix3Xi &tf)
{
    std::ofstream os(out_file);
    if(!os)
        return 0;
    os <<"mtllib ./make_human.mtl"<<std::endl;
    for (int i = 0; i < V.cols(); ++i){
        os << "v " << V(0,i) << " " << V(1,i) << " " << V(2,i) << "\n";
    }
    for (size_t i = 0; i < tv.cols(); ++i){
        os << "vt " << tv(0, i) << " " << tv(1, i)<< "\n";
    }
    for (size_t i = 0; i < tf.cols(); ++i){
        os << "f " << F(0, i) + 1 << "/" << tf(0, i) + 1 << " "
           << F(1, i)  + 1 << "/" << tf(1, i) + 1 << " "
           << F(2, i)  + 1 << "/" << tf(2, i) + 1 << "\n";
    }
    os.close();
    return 1;
}

}


#endif
