#ifndef OTOO_EOS_H
#define OTOO_EOS_H
#include "OTOO.hpp"
namespace OTOO {
  class EOS {
  public:
    EOS() {};
    ~EOS() {};

    virtual float GetS(float d, float e) = 0;
    virtual float GetP(float d, float e) = 0;
    virtual float GetE(float d, float e) = 0;
    virtual float GetT(float d, float e) = 0;
    virtual void TestEOS(uint64 n, float *d, float *e= NULL) = 0;
    virtual float GetEmin(float d) = 0; 
    virtual float GetEmin2(float d) = 0; 
  };
}
#endif
