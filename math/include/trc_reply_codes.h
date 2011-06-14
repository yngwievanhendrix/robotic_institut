/**

 \file trc_reply_codes.h
 (c) 2007 Institut f�r Robotik, Universit�t zu L�beck, +49-(0)451-5005201
 
 PROJECT:        Tracking system communication<br>
 COMPONENT:      trc_reply_codes.h
 \brief Error and exit codes for tracking system communication
 
 **************************************************************************/

#ifndef __trc_replycodes
#define __trc_replycodes

#include <string>

namespace Tracking {

	//! Tracking exit code data type
	typedef int TRC_ECODE;

	//! Tracking exit codes
	typedef enum {
		//! Success
		TRC_SUCCESS = 1,

		//! The tracking class is already connected
		TRC_CON_ALREADY_CONNECTED 	=	210001,
		//! The connection has failed (on the socket level)
		TRC_CON_CONNECTION_FAILED	=	210002,
		//! Communication with the tracking server has failed
		TRC_CON_COMMUNICATION_FAILED=	210003,
		//! Unknown tracking system type
		TRC_CON_UNKNOWN_TRACKER_TYPE=	210004,
		//! Sending of command failed
		TRC_CON_COMMAND_FAILED		=	210005,
		//! Not connected
		TRC_CON_NOT_CONNECTED		=	210006,
		//! Not initialised
		TRC_CON_NOT_INITIALISED		=	210007,

		//! Invalid vector dimension
		TRC_VEC_INVALID_DIM			=	220001,
		//! Invalid matrix dimension
		TRC_MAT_INVALID_DIM			=	220002,

		//! Thread not running
		TRC_THR_NOT_RUNNING			=	230001,
		//! Thread already running
		TRC_THR_ALREADY_RUNNING		=	230002,
		//! Couldn't start or resume thread
		TRC_THR_START_FAILED		=	230003,
		//! Tracking thread is NULL
		TRC_THR_NULL				=	230004,
		//! Tracking thread has not handshaked yet
		TRC_THR_NOT_HANDSHAKED		=	230005,
		//! No tracking values available yet
		TRC_THR_NO_VALUES_YET		=	230006,

		//! Marker is not visible
		TRC_MRK_INVISIBLE			=	240001,

		//! Invalid averaging mode selected
		TRC_CMD_INVALID_MODE		=	250001,
		//! Invalid logging level
		TRC_CMD_INVALID_LEVEL		=	250002,
		//! Command not supported
		TRC_CMD_NOT_SUPPORTED       =   250003,

		//! Invalid data type queried
		TRC_DAT_INVALID_TYPE		=   260001
	} TRC_ECODES;
}

#endif
