#include "math/include/ecodes.h"

std::string ToString(ECODE code) {

	using namespace Tracking;
	using namespace Robot;
	using namespace Calibration;

	switch(code) {
		case TRC_SUCCESS:
		// case CAL_SUCCESS: -- all success cases have the same ENUM value - 1
		// case ROB_SUCCESS: -- all success cases have the same ENUM value - 1
			return "No error";

			/* TRACKING CODES */
		case TRC_CON_ALREADY_CONNECTED:
			return "The tracking class is already connected";
		case TRC_CON_CONNECTION_FAILED:
			return "The connection has failed (on the socket level)";
		case TRC_CON_COMMUNICATION_FAILED:
			return "Communication with the tracking server has failed";
		case TRC_CON_UNKNOWN_TRACKER_TYPE:
			return "Unknown tracking system type";
		case TRC_CON_COMMAND_FAILED:
			return "Sending of command failed";
		case TRC_CON_NOT_CONNECTED:
			return "Not connected";
		case TRC_CON_NOT_INITIALISED:
			return "Not initialised";
		case TRC_VEC_INVALID_DIM:
			return "Invalid vector dimension";
		case TRC_MAT_INVALID_DIM:
			return "Invalid matrix dimension";
		case TRC_THR_NOT_RUNNING:
			return "Thread not running";
		case TRC_THR_ALREADY_RUNNING:
			return "Thread already running";
		case TRC_THR_START_FAILED:
			return "Couldn't start or resume thread";
		case TRC_THR_NULL:
			return "Tracking thread is NULL";
		case TRC_MRK_INVISIBLE:
			return "Marker is not visible";
		case TRC_CMD_INVALID_MODE:
			return "Invalid averaging mode selected";
		case TRC_CMD_INVALID_LEVEL:
			return "Invalid logging level";
		case TRC_CMD_NOT_SUPPORTED:
			return "Command not supported";
		case TRC_THR_NOT_HANDSHAKED:
			return "Tracking thread has not handshaked yet";
		case TRC_THR_NO_VALUES_YET:
			return "No tracking values available yet";
		case TRC_DAT_INVALID_TYPE:
			return "Invalid data type queried";

			/* CALIBRATION CODES */
		case CAL_COM_ROBOT_NOT_WORKING:
			return "Communication problems with the robot";
		case CAL_COM_TRACKING_NOT_WORKING:
			return "Communication problems with the tracking system";

		case CAL_COLLECT_FAILED:
			return "Calibration: collecting points failed";
		case CAL_FILE_OPEN_FAILED:
			return "Calibration: opening output file failed";
		case CAL_FILE_WRITE_FAILED:
			return "Calibration: writing to the output file failed";
		case CAL_FILE_READ_FAILED:
			return "Calibration: reading from the input file failed";
		case CAL_NOT_CALIBRATED:
			return "Calibration: not calibrated";
		case CAL_FILE_NO_FILES:
			return "No files specified";
		case CAL_NOT_INITIALISED:
			return "No points specified";
		case CAL_ALG_DOESNT_EXIST:
			return "Algorithm does not exist";
		case CAL_ALG_POINTS_DONT_EXIST:
			return "Points do not exist";
		case CAL_XML_PARSING_FAILED:
			return "XML read failed";
		case CAL_XML_INCORRECT_NR_OF_SYSTEMS:
			return "Incorrect number of trackingsystems in XML file (need 2)";
		case CAL_XML_INCORRECT_NR_OF_TRACKERS:
			return "Incorrect number of trackers in XML file (need 2)";
		case CAL_XML_NO_ROBOT_PRESENT:
			return "No robot found in calibration file";
		case CAL_XML_INCORRECT_NR_OF_MATRICES:
			return "Number of calibration matrices is wrong (should be 2)";
		case CAL_XML_TOO_MANY_ROBOTS:
			return "Too many robots in calibration file (should be 1)";
		case CAL_UNSUPPORTED_ROTTYPE:
			return "Currently not supported rotation type detected";
		case CAL_CONVERSION_NOT_SUPPORTED:
			return "Currently not supported rotation conversion";
		case CAL_VEC_INVALID_INDEX:
			return "Invalid index";
		case CAL_MAT_NO_MATRICES:
			return "No matrices available";
		case CAL_ALG_DOESNT_SUPPORT_DTYPE:
			return "Algorithm does not support the data type";
		case CAL_SYS_INVALID:
			return "Invalid system index";

			/* ROBOT CODES */
		case ROB_CON_ALREADY_CONNECTED:
			return "The robot class is already connected";
		case ROB_CON_CONNECTION_FAILED:
			return "The connection has failed (on the socket level)";
		case ROB_CON_COMMUNICATION_FAILED:
			return "Communication with the robot server has failed";
		case ROB_CON_UNKNOWN_ROBOT_TYPE:
			return "Unknown robot type";
		case ROB_CON_WRONG_ROBOT_TYPE:
			return "Wrong robot type";
		case ROB_CON_COMMAND_FAILED:
			return "Sending of command failed";
		case ROB_CON_NOT_CONNECTED:
			return "Not connected";
		case ROB_MOV_FAILED:
			return "Movement failed";
		case ROB_MOV_NOT_IN_PTP:
			return "Not in PTP mode";
		case ROB_MOV_NOT_IN_RT:
			return "Not in real time mode";
		case ROB_VEC_INVALID_DIM:
			return "Invalid vector dimension";
		case ROB_MAT_INVALID_DIM:
			return "Invalid matrix dimension";
		case ROB_POS_KUKA_RT_WARN:
			return "Warning, if position of Kuka robots is queried in RT mode";
		case ROB_POS_FAILED:
			return "Querying the position has failed";
		case ROB_POS_QUAT_FAILED:
			return "Computing the quaternion has failed";
		case ROB_SET_FAILED:
			return "Setting a parameter failed";
		case ROB_SET_INVALID_JOINT:
			return "Invalid joint selected";
		case ROB_SET_INVALID_SPEED:
			return "Invalid speed";
		case ROB_SET_INVALID_ACCELERATION:
			return "Invalid acceleration";
		case ROB_SET_WRONG_MODE:
			return "Not in correct operations mode";
		case ROB_SET_WRONG_ROBOT:
			return "Duration of movement can only be set for adept";
		case ROB_INVALID_COMMAND:
			return "Invalid command";
		case ROB_GET_FAILED:
			return "Getting a parameter failed";

		default:
			return "Unknown exit code";
	}
}

