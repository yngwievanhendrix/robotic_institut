/**

 \file rob_reply_codes.h
 (c) 2007 Institut für Robotik, Universität zu Lübeck, +49-(0)451-5005201
 
 PROJECT:        Robot communication<br>
 COMPONENT:      rob_reply_codes.h
 \brief Error and exit codes for robot communication
 
 **************************************************************************/

#ifndef __rob_replycodes
#define __rob_replycodes

namespace Robot {
	//! Robot exit code data type
	typedef int ROB_ECODE;

	//! Robot exit codes
	typedef enum {
		//! Success
		ROB_SUCCESS					=	1,

		//! The robot class is already connected
		ROB_CON_ALREADY_CONNECTED	=	110001,
		//! The connection has failed (on the socket level)
		ROB_CON_CONNECTION_FAILED	=	110002,
		//! Communication with the robot server has failed
		ROB_CON_COMMUNICATION_FAILED=	110003,
		//! Unknown robot type
		ROB_CON_UNKNOWN_ROBOT_TYPE	=	110004,
		//! Sending of command failed
		ROB_CON_COMMAND_FAILED		=	110005,
		//! Not connected
		ROB_CON_NOT_CONNECTED		=	110006,
		//! Wrong robot type
		ROB_CON_WRONG_ROBOT_TYPE	=	110007,

		//! Movement failed
		ROB_MOV_FAILED				=	120001,
		//! Not in PTP mode
		ROB_MOV_NOT_IN_PTP			=	120002,
		//! Not in real time mode
		ROB_MOV_NOT_IN_RT			=	120003,

		//! Invalid vector dimension
		ROB_VEC_INVALID_DIM			=	130001,
		//! Invalid matrix dimension
		ROB_MAT_INVALID_DIM			=	130002,

		//! Warning, if position of Kuka robots is queried in RT mode
		ROB_POS_KUKA_RT_WARN		=	140001,
		//! Querying the position has failed
		ROB_POS_FAILED				=	140002,
		//! Computing the quaternion has failed
		ROB_POS_QUAT_FAILED			=	140003,

		//! Setting a parameter failed
		ROB_SET_FAILED				=	150000,
		//! Invalid joint selected
		ROB_SET_INVALID_JOINT		=	150001,
		//! Invalid speed
		ROB_SET_INVALID_SPEED		=	150002,
		//! Invalid acceleration
		ROB_SET_INVALID_ACCELERATION=	150003,
		//! Not in correct operations mode
		ROB_SET_WRONG_MODE			=	150004,
		//! Duration of movement can only be set for adept
		ROB_SET_WRONG_ROBOT			=	150005,

		//! Invalid command
		ROB_INVALID_COMMAND			=	160000,

		//! Getting a parameter failed
		ROB_GET_FAILED				=	170000,

		//! Gripper not present
		ROB_GRIPPER_NOT_PRESENT		=   180000,

		//! Gripper command failed
		ROB_GRIPPER_FAILED          =   180001

	} ROB_ECODES;
}

#endif
