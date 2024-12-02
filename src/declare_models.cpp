// ======================================================================
// Model declaration (gets compiled in models.cpp)
// model headers and declare_model must be consistent!


// model headers

#include "models/drylyotropicF.hpp"

//#include "models/activemodelhplus.hpp"

void DeclareModels()
{
        declare_model<DryLyotropicF>(
      "drylyotropicF",
      "Biphasic, lyotropic, nematic model as presented as described in "
      "10.1103/PhysRevLett.113.248303. We refer the user to this reference "
      "for further information. Extended to dry limit and with polarity "
      "alignment to the velocity. "
      );

}
