#ifndef __ecodes_h
#define __ecodes_h

#include "math/include/trc_reply_codes.h"
#include "math/include/cal_reply_codes.h"
#include "math/include/rob_reply_codes.h"

typedef int ECODE;

std::string ToString(ECODE code);

#endif

