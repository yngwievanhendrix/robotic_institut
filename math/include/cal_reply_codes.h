/**

 \file cal_reply_codes.h
 (c) 2007 Institut für Robotik, Universität zu Lübeck, +49-(0)451-5005201
 
 PROJECT:        Tracking system calibration<br>
 COMPONENT:      rob_reply_codes.h
 \brief Error and exit codes for calibration
 
 **************************************************************************/

#ifndef __cal_reply_codes_h
#define __cal_reply_codes_h

namespace Calibration {
	//! Calibration class exit code data type
	typedef int CAL_ECODE;

	//! Calibration class exit codes
	typedef enum {
		//! Success
		CAL_SUCCESS					=	1,

		//! Communication problems with the robot
		CAL_COM_ROBOT_NOT_WORKING	=	310001,
		//! Communication problems with the tracking system
		CAL_COM_TRACKING_NOT_WORKING=	310002,

		//! Calibration: collecting points failed
		CAL_COLLECT_FAILED			=	320001,
		//! Calibration: opening output file failed
		CAL_FILE_OPEN_FAILED		=	320002,
		//! Calibration: writing to the output file failed
		CAL_FILE_WRITE_FAILED		=	320003,
		//! Calibration: reading from the input file failed
		CAL_FILE_READ_FAILED		=	320004,
		//! Calibration: not calibrated
		CAL_NOT_CALIBRATED			=	320005,
		//! No files specified
		CAL_FILE_NO_FILES			=	320006,
		//! No points specified
		CAL_NOT_INITIALISED			=	320007,
		//! Algorithm does not exist
		CAL_ALG_DOESNT_EXIST		=	320008,
		//! Points do not exist
		CAL_ALG_POINTS_DONT_EXIST	=   320009,
		//! XML read failed
		CAL_XML_PARSING_FAILED		=	320010,
		//! Incorrect number of trackingsystems in XML file (need 2)
		CAL_XML_INCORRECT_NR_OF_SYSTEMS = 320011,
		//! Incorrect number of trackers in XML file (need 2)
		CAL_XML_INCORRECT_NR_OF_TRACKERS = 320012,
		//! No robot found in calibration file
		CAL_XML_NO_ROBOT_PRESENT	=	320013,
		//! Number of calibration matrices is wrong (should be 2)
		CAL_XML_INCORRECT_NR_OF_MATRICES	=	320014,
		//! Too many robots in calibration file (should be 1)
		CAL_XML_TOO_MANY_ROBOTS		=	320015,
		//! Currently not supported rotation type detected
		CAL_UNSUPPORTED_ROTTYPE		=	320016,
		//! Currently not supported rotation conversion
		CAL_CONVERSION_NOT_SUPPORTED=	320017,
		//! Invalid index
		CAL_VEC_INVALID_INDEX       =   320018,
		//! No matrices available
		CAL_MAT_NO_MATRICES         =   320019,
		//! Algorithm does not support the data type
		CAL_ALG_DOESNT_SUPPORT_DTYPE =	320020,
		//! Invalid system index
		CAL_SYS_INVALID             =   320021
	} CAL_ECODES;
}

#endif
