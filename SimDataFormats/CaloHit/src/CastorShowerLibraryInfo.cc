#include "SimDataFormats/CaloHit/interface/CastorShowerLibraryInfo.h"
#include <iostream>

ClassImp(CastorShowerLibraryInfo)

CastorShowerLibraryInfo::CastorShowerLibraryInfo() {
   // Clear();
   // std::cout << "\n    *** CastorShowerLibraryInfo object created ***    " << std::endl;
}

CastorShowerLibraryInfo::~CastorShowerLibraryInfo() {}

    
void CastorShowerLibraryInfo::Clear() {
   Energy.Clear();
   Eta.Clear();
   Phi.Clear();
}
