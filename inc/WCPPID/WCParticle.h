#ifndef WCPPID_WCPARTICLE_H
#define WCPPID_WCPARTICLE_H

#include <vector>

namespace WCPPID{
  class WCParticle{
  public:
    WCParticle(int particle_id=0);
    ~WCParticle();

    void set_particle_id(int id){particle_id = id;};
    int get_particle_id(){return particle_id;};
    
  protected:
    int particle_id;
  };
  
  typedef std::vector<WCParticle*> WCParticleSelection;
}

#endif
