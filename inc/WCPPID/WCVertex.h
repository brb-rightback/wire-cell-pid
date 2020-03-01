#ifndef WCPPID_WCVERTEX_H
#define WCPPID_WCVERTEX_H

#include <vector>

namespace WCPPID{
  class WCVertex{
  public:
    WCVertex(int vertex_id=0);
    ~WCVertex();

    void set_vertex_id(int id){vertex_id=id;};
    int get_vertex_id(){return vertex_id;};
    
  protected:
    int vertex_id;
  };

  typedef std::vector<WCVertex*> WCVertexSelection;
}

#endif
