
#include "TMath.h"

Float_t deltaPhi(Float_t phi1, Float_t phi2) {
    Float_t substr = abs(phi1 - phi2);
    if (substr <= TMath::Pi()) {
        return substr;
    } else {
        //(substr > TMath::Pi() && substr <= 2 * TMath::Pi()) {
        return 2 * TMath::Pi() - substr;
    }
}

// Float_t deltaPhi(Float_t phi1, Float_t phi2) {
// Float_t substr = abs(phi1 - phi2);
// while (true) {
// if (substr < TMath::Pi()) {
// return substr;
//} else {
// substr -= TMath::Pi();
//}
//}
//}
