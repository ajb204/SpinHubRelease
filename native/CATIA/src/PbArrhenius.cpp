/*
 *  PbStandard.cpp
 *
 *
 *  Created by Guillaume Bouvignies on 25/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <Catia.h>
#include <Dataset.h>

#define PLANCKS_H (6.62606896e-34) /* kg m^2 / s */
#define BOLTZMANN_KB (1.3806504e-23) /* kg m^2 / K s^2 */
#define MOLAR_GAS_R (8.314472e0) /* kg m^2 / K mol s^2 */

double Catia::PbArrhenius(const Dataset& dset) const {
    const double deltaSb = ResolveParam(GlobalParam,"deltaSb");
    const double deltaHb = ResolveParam(GlobalParam,"deltaHb");
    const double deltaSab = ResolveParam(GlobalParam,"deltaSab");
    const double deltaHab = ResolveParam(GlobalParam,"deltaHab");


    const double T = dset._temperature + 273.15;

    const double A = BOLTZMANN_KB * T / PLANCKS_H;
    const double RT = MOLAR_GAS_R * T;

    const double pA = 1.0;
    const double pB = exp( -(deltaHb - T * deltaSb) / RT );

    const double pb = pB / (pA + pB);

    /*
    const double deltaGb=(deltaHb-T*deltaSb)/RT;
    const double deltaGab=(deltaHab-T*deltaSab)/RT;
    const double deltaGba=(deltaGb-deltaGab)/RT;



    const double tits1=exp((-(deltaHb-T*deltaSb)+(deltaHab-T*deltaSab))/RT);
    const double tits2=exp((-(deltaHab-T*deltaSab))/RT);
    const double tits3=tits1/(tits1+tits2);
    const double tits4=tits2/(tits1+tits2);


    std::cout<< "     deltaH: " << deltaHb << std::endl;
    std::cout<< "     deltaS: " << deltaSb << std::endl;
    std::cout<< "     deltaG: " << deltaHb-T*deltaSb << std::endl;
    std::cout<< "    deltaGb : " << deltaGb << std::endl;
    std::cout<< "    deltaGab: " << deltaGab << std::endl;
    std::cout<< "    deltaGba: " << deltaGba << std::endl;


    std::cout<< "  deltaG/RT: " << (deltaHb-T*deltaSb)/RT << std::endl;
    std::cout<< "         pA: " << pA << std::endl;
    std::cout<< "         pB: " << pB << std::endl;
    std::cout<< "         pb: " << pb << std::endl;
    std::cout<< "          T: " << T << std::endl;

    std::cout<< "       tits1: " << tits1 << std::endl;
    std::cout<< "       tits2: " << tits2 << std::endl;
    std::cout<< "       tits3: " << tits3 << std::endl;
    std::cout<< "       tits4: " << tits4 << std::endl;

*/

    return pb;
}

double Catia::PbArrhenius_3st(const Dataset& dset) const {

    const double deltaSb = ResolveParam(GlobalParam,"deltaSb");
    const double deltaSc = ResolveParam(GlobalParam,"deltaSc");

    const double deltaHb = ResolveParam(GlobalParam,"deltaHb");
    const double deltaHc = ResolveParam(GlobalParam,"deltaHc");

    const double T = dset._temperature + 273.15;
	const double RT = MOLAR_GAS_R * T;

    const double pA = 1.0;
    const double pB = exp( -(deltaHb - T * deltaSb) / RT );
    const double pC = exp( -(deltaHc - T * deltaSc) / RT );

    const double pb = pB / (pA + pB + pC);

    return pb;

}

double Catia::PcArrhenius_3st(const Dataset& dset) const {

    const double deltaSb = ResolveParam(GlobalParam,"deltaSb");
    const double deltaSc = ResolveParam(GlobalParam,"deltaSc");

    const double deltaHb = ResolveParam(GlobalParam,"deltaHb");
    const double deltaHc = ResolveParam(GlobalParam,"deltaHc");

    const double T = dset._temperature + 273.15;
	const double RT = MOLAR_GAS_R * T;

    const double pA = 1.0;
    const double pB = exp( -(deltaHb - T * deltaSb) / RT );
    const double pC = exp( -(deltaHc - T * deltaSc) / RT );

    const double pc = pC / (pA + pB + pC);

    return pc;

}
