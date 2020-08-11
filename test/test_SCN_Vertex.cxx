#include "WCPPyUtil/SCN_Vertex.h"
#include "WCPPyUtil/read_numpy.h"

#include <iostream>

int main()
{
    typedef float VType;

    // std::vector<VType> x = {00.0, 01.0, 02.0};
    // std::vector<VType> y = {10.0, 11.0, 12.0};
    // std::vector<VType> z = {20.0, 21.0, 22.0};
    // std::vector<VType> q = {1.0, 2.0, 3.0};

    const std::string data_path = "../wire-cell-pydata/scn_vtx/nuecc-sample.npz";
    auto x = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "coords", 0, true);
    auto y = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "coords", 1);
    auto z = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "coords", 2);
    auto q = WCPPyUtil::read_npz("read_numpy", "read_npz", data_path, "ft", 0);

    const std::string module = "SCN_Vertex";
    const std::string function = "SCN_Vertex";
    //const std::string function = "Test";
    const std::string weights =
        "./input_data_files/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth";

    
    auto vtx = WCPPyUtil::SCN_Vertex(module, function, weights, {x, y, z, q}, "float32", true);

    return 0;
}
