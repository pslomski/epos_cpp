#ifndef LIB_EPOS_FIT_HPP_INCLUDED
#define LIB_EPOS_FIT_HPP_INCLUDED

#include <nr.hpp>
#include <epos.hpp>

namespace epos {

	bool fit_fun(VectorFloat& x, VectorFloat& y, VectorFloat& a, VectorFloat& sig_a, fitfun_type function);

	class FitParams {
		epos::Float x;
		epos::Float y;
		epos::Float a0;
	};

} // namespace epos

#endif // LIB_EPOS_FIT_HPP_INCLUDED