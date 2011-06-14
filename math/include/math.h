#ifndef __math_h
#define __math_h

#include "math/include/wavelets.h"
#include "math/include/cal_algorithms.h"
#include "math/include/dataobject.h"
#include "math/include/Horn.h"
#include "math/include/matrix.h"
#include "math/include/vectorops.h"
#include "math/include/ecodes.h"
#include "math/include/cal_reply_codes.h"
#include "math/include/trc_reply_codes.h"
#include "math/include/rob_reply_codes.h"
#include "math/include/distortion.h"

namespace Math {
	inline std::string IntToString(int x) {
		char buf[128];
		sprintf(buf, "%d", x);
		return buf;
	}
}

#endif
