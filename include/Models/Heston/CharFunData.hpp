/**
* CharFunData.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Utils/Types.hpp"

namespace uv::models::heston
{
	struct CharFunData
	{
		// Original CF outputs
		Complex<Real> psi;                 // psi(u) := exp( A + v0 * B )
		Complex<Real> A;                   // A := (kappa*theta/sigma^2)*( r*T − 2*log(1 − r*y) )
		Complex<Real> B;                   // B := u(u+i)*y / (1 − r*y)
		Complex<Real> beta;                // beta := kappa − sigma*rho*(i*u)
		Complex<Real> D;                   // D := sqrt( beta^2 + sigma^2*u(u+i) )
		Complex<Real> DT;                  // DT := D*T
		Complex<Real> betaPlusD;           // beta + D
		Complex<Real> betaMinusD;          // beta − D  (stable)
		Complex<Real> ui;                  // ui := u*i
		Complex<Real> kFac;                // kFac := (kappa*theta)/sigma^2
		Real invSigma2;					   // 1 / sigma^2
		Real kappaTheta;				   // kappa * theta
		Real sigma2;					   // sigma^2

		// Rescued intermediates (for gradient path)
		Complex<Real> uu;                  // uu := u(u+i)
		Complex<Real> eDT;                 // eDT := exp( −DT )
		Complex<Real> g;                   // g := (beta−D)/(beta+D)
		Complex<Real> Q;                   // Q := 1 − g*eDT
		Complex<Real> invQ;                // 1 / Q
		Complex<Real> invQ2;               // 1 / Q^2
		Complex<Real> R;                   // R := 1 − g
		Complex<Real> S;                   // S := (beta−D)*T − 2*log(Q/R)
		Complex<Real> fracB;               // fracB := (1 − eDT) / Q
		Complex<Real> denomG;              // denomG := (beta + D)^2
		Complex<Real> betaMinusDinvSigma2; // (betaMinusD)/sigma^2
	};
}
