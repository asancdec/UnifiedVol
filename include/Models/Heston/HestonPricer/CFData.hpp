/**
* CFData.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <complex>

namespace uv
{
	using cplx = ::std::complex<long double>;

	struct CFData
	{
		// Original CF outputs
		cplx psi;                 // psi(u) := exp( A + v0 * B )
		cplx A;                   // A := (kappa*theta/sigma^2)*( r*T − 2*log(1 − r*y) )
		cplx B;                   // B := u(u+i)*y / (1 − r*y)
		cplx beta;                // beta := kappa − sigma*rho*(i*u)
		cplx D;                   // D := sqrt( beta^2 + sigma^2*u(u+i) )
		cplx DT;                  // DT := D*T
		cplx betaPlusD;           // beta + D
		cplx betaMinusD;          // beta − D  (stable)
		cplx ui;                  // ui := u*i
		cplx kFac;                // kFac := (kappa*theta)/sigma^2
		long double invSigma2;    // 1 / sigma^2
		long double kappaTheta;   // kappa * theta
		long double sigma2;       // sigma^2

		// Rescued intermediates (for gradient path)
		cplx uu;                  // uu := u(u+i)
		cplx eDT;                 // eDT := exp( −DT )
		cplx g;                   // g := (beta−D)/(beta+D)
		cplx Q;                   // Q := 1 − g*eDT
		cplx invQ;                // 1 / Q
		cplx invQ2;               // 1 / Q^2
		cplx R;                   // R := 1 − g
		cplx S;                   // S := (beta−D)*T − 2*log(Q/R)
		cplx fracB;               // fracB := (1 − eDT) / Q
		cplx denomG;              // denomG := (beta + D)^2
		cplx betaMinusDinvSigma2; // (betaMinusD)/sigma^2
	};
}
