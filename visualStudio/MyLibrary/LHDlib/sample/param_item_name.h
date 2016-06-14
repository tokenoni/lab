#ifndef PARAM_ITEM_NAME_H_
#define PARAM_ITEM_NAME_H_

/**  Channel Common Parameter's name
*/
#define ARC_PARAM_DIAG_NAME		"DiagnosticsName"
#define ARC_PARAM_SHOT			"Shot"
#define ARC_PARAM_SUB_SHOT		"SubShot"
#define ARC_PARAM_DATATYPE		"DataType"
#define ARC_PARAM_MOD_GROUP		"ModuleGroup"
#define ARC_PARAM_MOD_TYPE		"ModuleType"
#define ARC_PARAM_MANAGE		"ManagementVersion"

#define ARC_PARAM_CH			"ChannelNumber"
#define ARC_PARAM_IMAGE_TYPE	"ImageType"
#define ARC_PARAM_DATA_LEN		"DataLength(byte)"
#define ARC_PARAM_COMP_LEN		"CompLength(byte)"
#define ARC_PARAM_COMP_METHOD	"CompressionMethod"
#define ARC_PARAM_FRAME			"SampledFrame"
#define ARC_PARAM_FRAME_SIZE	"FrameByteSize"
#define ARC_PARAM_FRAME_X		"AcquisitionWindowWidth"
#define ARC_PARAM_FRAME_Y		"AcquisitionWindowHeight"
#define ARC_PARAM_CRC32			"CRC32"
#define ARC_PARAM_COMMENT		"Comment"

#define ARC_PARAM_DTS_SOURCE	"DTSsource"
#define ARC_PARAM_DTS_HOST		"DTShostID"
#define ARC_PARAM_DTS_MODULE	"DTSmoduleID"
#define ARC_PARAM_DTS_TRIGGER	"DTStriggerChannel"
#define ARC_PARAM_DTS_CLOCK		"DTSClockChannel"
#define ARC_PARAM_DTS_USER		"DTSuserDefinedClock"
#define ARC_PARAM_SIGNAL_NAME	"SignalName"

#define ARC_PARAM_CLK_SOURCE	"CLKsource"
#define ARC_PARAM_CLK_HOST		"CLKhostID"
#define ARC_PARAM_CLK_MODULE	"CLKmoduleID"

//#define ARC_PARAM_ADLER32		"ADLER32"

/**  Channel Parameter's name
*/
#define OTH_PARAM_VALUE_LEN		"Resolution(bit)"

/**  Shot Parameter's name
*/
#define SHOT_PARAM_DIAG_NAME 	"DiagnosticsName"
#define SHOT_PARAM_SHOT			"Shot"
#define SHOT_PARAM_SUB_SHOT		"SubShot"
#define SHOT_PARAM_MANAGE		"ManagementVersion"
#define SHOT_PARAM_MOD_GROUP	"ModuleGroup"
#define SHOT_PARAM_MOD_NUM		"ModuleNumber"
#define SHOT_PARAM_CH_NUM		"ChannelNumber"
#define SHOT_PARAM_COL_CH_NUM	"CollectedChannel"
#define SHOT_PARAM_DATA_LEN		"TotalDataLength"
#define SHOT_PARAM_COMP_LEN		"TotalCompLength"
#define SHOT_PARAM_ACQ_DATE		"AcquisitionDate"
#define SHOT_PARAM_ACQ_DATE_X	"AcquistionDate"
#define SHOT_PARAM_ARC_DATE		"ArchiveDate"
#define SHOT_PARAM_COMMENT		"Comment"
#define SHOT_PARAM_SITE_NAME 	"SiteName"

/**  Parameter's value
*/
#define DBS_MOD_GROUP		"ANY"
#define DBS_MOD_TYPE		"General"
#define DBS_MANAGE_VER		"10.0.0"

/**  CAMAC Parameter's name
*/
#define CAMAC_PARAM_LAM					"LAM"

/**  CPHA1219 Parameter's name
*/
#define CPHA1219_PARAM_CH_P_LIST		"Channel/Page"
#define CPHA1219_PARAM_TOTAL_LISTS		"TotalPages"
#define CPHA1219_PARAM_INTERVAL_SRC		"PagingClock"
#define CPHA1219_PARAM_INTERVAL_TIME	"ClockInterval(mSec)"
#define CPHA1219_PARAM_M_STATUS			"ModuleStatus"

/**  WE Parameter's name
*/
#define WE_PARAM_ST_NO				"StationNumber"
#define WE_PARAM_SLOT_NO			"SlotNumber"
#define WE_PARAM_CH_NO_MD			"ChannelNumberInModule"
#define WE_PARAM_TRIG_SRC			"TriggerSource"
#define WE_PARAM_CH_STATUS			"Status"
#define WE_PARAM_PHA_OPE_MODE		"OperationMode"
#define WE_PARAM_PHA_MES_MODE		"MeasureMode"
#define WE_PARAM_PHA_LINKED_IN		"LinkedInput"
#define WE_PARAM_PHA_CH_PAGE		"Channel/Page"
#define WE_PARAM_PHA_WIN_TIME		"WindowTime(uSec)"
#define WE_PARAM_PHA_ERR_CHK		"ErrorCheck"
#define WE_PARAM_PHA_ROI			"ROI_Enabled"
#define WE_PARAM_PHA_ROI_LLD		"ROI_LowerLimit"
#define WE_PARAM_PHA_ROI_ULD		"ROI_UpperLimit"
#define WE_PARAM_PHA_PAGING_SRC		"PagingClock"
#define WE_PARAM_PHA_PAGING_CLK		"ClockInterval(mSec)"
#define WE_PARAM_PHA_TRIG_SEL		"TriggerSelect"
#define WE_PARAM_PHA_GATING			"Gating"
#define WE_PARAM_PHA_BIT_LENGTH		"CounterLength"
#define WE_PARAM_PHA_TOTAL_PAGES	"TotalPages"
#define WE_PARAM_PHA_COUNT_SAMP		"CountOfSampledData"
#define WE_PARAM_OVER_RUN			"OverRun"

#define WE_PARAM_PHA_FR_TRIG_EVENT		"TriggerEvents"
#define WE_PARAM_PHA_FR_MEASURE_TM		"MeasuringTime(mSec)"
#define WE_PARAM_PHA_FR_ERROR_EVENT		"ErrorEvents"

/**  VME FPGA DEMOD Parameter's name
*/
#define DTS_PARAM_DTS_STATUS	"Status"				// Up/Down
#define DTS_PARAM_CARD_TYPE		"CardType"
#define DTS_PARAM_COL_TIME		"TimeOfCollected"		//,"Wed Jan 16 18:26:45 2008",1
#define DTS_PARAM_REG_EVENT 	"EventOutMode"			//,On,1
#define DTS_PARAM_BASE_RATE 	"BaseRate"				//,1000000,4
#define DTS_PARAM_EXT_TRIG		"ExternalTrigger"		//,Valid,1
#define DTS_PARAM_CLOCKSRC		"ClockSource"			//,External,1
#define DTS_PARAM_INHIBIT		"Inhibit"				// DMODSS
#define DTS_PARAM_BUSIF			"BusIF"					// DMODSS
#define DTS_PARAM_VME_INT		"VMEbusInterrupt"		// DMOD ,Disable
#define DTS_PARAM_MODE 			"MoveMode"				//,MODE0,1

#define DTS_PARAM_INT_TRIG		"InterruptTrigger"		//,Disable,1
#define DTS_PARAM_INT_EVENT		"InterruptEvent"		//,Disable,1
#define DTS_PARAM_INT_UNH		"InterruptUninhibit"	//,Disable,1
#define DTS_PARAM_INT_INH		"InterruptInhibit"		//,Disable,1
#define DTS_PARAM_INT_ERROR		"InterruptError"		//,Disable,1
#define DTS_PARAM_INT_CLK_ERR	"InterruptClockError"	//,Disable,1
#define DTS_PARAM_INT_SETUP		"InterruptSetup"		//,Disable,1
#define DTS_PARAM_INT_STOP		"InterruptStop"			//,Disable,1

#define DTS_PARAM_TIMER_TRIG	"TimerTriggerChannel"	//,1,4
#define DTS_PARAM_INTER_DELAY	"DelayLinePreset"		//,0,4
#define DTS_PARAM_MECH_DELAY	"InternalDelay"			//,15,4

#define DTS_PARAM_CH_NO_IN		"ChannelNumberInModule"	//,1,4 

#define DTS_PARAM_DVD_RANGE		"DvdRange"				//,1usec,1
#define DTS_PARAM_DVD_MAG		"DvdMagnification"		//,2,4
#define DTS_PARAM_DVD_SIG		"DvdSignal"				// for VME

#define DTS_PARAM_DVD_CON_DIAG	"DvdConnectDiag."		//,SX8O,1

#define DTS_PARAM_DELAY_R		"DelayTime_s"			//,-1,6
#define DTS_PARAM_PULSE_R		"PulseWidth_s"			//,10,6
#define DTS_PARAM_REP_TIME_R	"RepetitionTime_s"		//,11.1,6
#define DTS_PARAM_DELAY			"DelayTime"				//,118999985,4
#define DTS_PARAM_PULSE			"PulseWidth"			//,10000000,4
#define DTS_PARAM_REP_TIME		"RepetitionTime"		//,11100000,4
#define DTS_PARAM_REP_CNT		"RepetitionCount"		//,1,4
#define DTS_PARAM_TRIG_CH		"TriggerSelect"			//,1,4

#define DTS_PARAM_CON_DIAG	"ConnectDiag."				//,Langmuir2,1
#define DTS_PARAM_HOST 		"HostName"					//,133.75.175.100,1

#endif /*PARAM_ITEM_NAME_H_*/
