#ifndef WCPPID_WCPARTICLE_H
#define WCPPID_WCPARTICLE_H

#include <vector>

namespace WCPPID{
  class WCParticle{
  public:
    WCParticle();
    ~WCParticle();
    
  protected:
    int particle_id;
  };
  
  typedef std::vector<WCParticle*> WCParticleSelection;
}

#endif
