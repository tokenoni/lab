/*
*
*/

#if defined(WIN32) || defined(WIN64)
 #ifndef _WIN32_WINNT
 #define _WIN32_WINNT	(0x0400)
 #endif
 #pragma warning(disable:4786)
#else
 #include <string>
 #include <errno.h>
#endif // WIN32

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "retrieve.h"

#define RAW_DATA	1
#define ANA_DATA	2

#define	startIndexNoClosing()	IndexNoClosing(1)
#define	endIndexNoClosing()		IndexNoClosing(0)

void Usage( const char* av0)
{
	printf(" Retrieve version %s, Library version %s, zlib %s\n\n",RETRIEVE_VERSION,retrieveVersion(),retrieveZversion());
	printf("Usage: %s diag shot subshot chno [filename] [options...]\n",av0);

	printf("\n== Primary arguements ==\n");
	printf("	(1) diag name ( Magnetics, Thomson, ...)\n");
	printf("	(2) shot number\n");
	printf("	(3) subshot number\n");
	printf("	(4) channel number(s) (e.g. 1,3,5-8  1,3,5:8  0=all )\n");
	printf("	         or a signal tag name (e.g. sig-tag1 )\n");
	printf("	(5) [ file or path(.../) name ]\n");
	printf("	    ( files will be \"'filename'-'shot'-'subshot'-'chno'.dat\" )\n");
	printf("	    ( or \"'pathname(.../)'+'diag'-'shot'-'subshot'-'chno'.dat\" )\n");

	printf("\n== Options ==\n");
	printf("	[-f frame number ( 0=all, 'start:end'=ranged )\n"
	       "		Framed 2-D or other data retrieval.\n");
	printf("	[-h 'Transd server name' -p 'base path' [-n port_no] ]\n"
	       "		Direct connection mode without referring the Index database.\n");
	printf("	[-w | --wait [-t timeout (sec)] ]\n"
	       "		Data waiting, i.e. blocking, mode. Default timeout: 180.\n");
	printf("	[-R | --real ]\n"
	       "		Files are restored with real shot+subshot numbers, not alias ones.\n");
	printf("	[-V | -VF | -V4 | --volts ] | [-VV | -VD | -V8 | --double-volts ]\n"
	       "		Data are restored into voltage values in FLOAT or DOUBLE precision.)\n");
	printf("	[ --silent ]\n"
	       "		No working messages output.\n");
	printf("	[-T | --timearray ]\n"
	       "		Make Files of Time Array and Time Parameter.\n");
	printf("	[-D | --double ]\n"
	       "		Make Time Array in double type (ignored if -T/--timearray option does not exist).\n");
	printf("	[ --original-name ]\n"
	       "		Files are restored with the original names.\n");
	printf("	[ --signal-name ]\n"
	       "		Restored file names are converted to 'signal-name' from 'chno'.\n");
}

int main(int argc, char* argv[])
{
	if(argc < 5 || argc > 21) {
		Usage( argv[0] );
		return -1;
	}
	int exit_status = 0;
	char*			Diag_Name		= argv[1];
	unsigned int	AShot_Number	= atoi(argv[2]);
	unsigned short	ASub_Shot_Number= atoi(argv[3]);
	char*			Channel			= argv[4];
	unsigned int	arc_shot;
	unsigned short	arc_sub;
	const char*		func_name;
	const char*		error_mes;

	char*	IndexServer	= getenv("INDEXSERVERNAME");
	char	Filename[256];
	char*	Frame		= NULL;
	char	Path[256]	={'\0'};
	char	Server[32]	={'\0'};
	int		PortNumber	=0;
	int		TimeoutSec	=180;
	bool	SwVerbose		= true;
	bool	bSilent			= false;
	bool	SwDirectAccess	= false;
	bool	SwWait			= false;
	bool	SwArcRealShot	= false;
	bool	SwVoltsFlt		= false;
	bool	SwVoltsDbl		= false;
	bool	SwOriginalName	= false;		// add 13.2.1
	bool	SwMakeTimeArray = false;		// add 14.1.1
	bool	SwDoubleTimeArray = false;		// add 14.1.1
	bool	SwSignalName	= false;		// add 16.0.0
	char*	path_for_orgname= NULL;			// add 13.2.1
	char	orgname[256];					// add 13.2.1
	char*	parg;
	int		sframe,eframe;
	int		ret;
	int		chno;
	char* 		signal_name= NULL;

	strcpy(Filename, Diag_Name);
	for(int i=5; i < argc; i++)
	{
		parg = argv[i];
		if(*parg=='-')
		{
			bool	swChkOpt1Len = true;
			switch(*(parg+1))
			{
				case 'f':
					i++;
					if( i < argc ) Frame=argv[i];
					break;
				case 'h':
					i++;
					SwDirectAccess = true;
					if( i < argc ) strcpy(Server, argv[i]);
					break;
				case 'p':
					i++;
					if( i < argc ) strcpy(Path, argv[i]);
					break;
				case 'n':
					i++;
					if( i < argc ) PortNumber=atoi(argv[i]);
					break;
				case 't':
					i++;
					if( i < argc ) TimeoutSec=atoi(argv[i]);
					break;
				case 'T':
					SwMakeTimeArray = true;
					break;
				case 'D':
					SwDoubleTimeArray = true;
					break;
				case 'w': // 2009-07-29
					SwWait = true;
					break;
				case 'R': // 2009-09-14
					SwArcRealShot = true;
					break;
				case '-': // 2009-09-14
					swChkOpt1Len = false;
					if( 0 == strcmp(parg+1,"-wait") ) {
						SwWait = true;
					}
					else if( 0 == strcmp(parg+1,"-real") ) {
						SwArcRealShot = true;
					}
					else if( 0 == strcmp(parg+1,"-volts") ) {
						SwVoltsFlt = true;
						break;
					}
					else if( 0 == strcmp(parg+1,"-double-volts") ) {
						SwVoltsDbl = true;
						break;
					}
					else if( 0 == strcmp(parg+1,"-silent") ) {
						SwVerbose = false;
						bSilent	= true;
						break;
					}
					else if( 0 == strcmp(parg+1,"-original-name") ) {
						SwOriginalName = true;
						break;
					} else if( 0 == strcmp(parg+1,"-timearray") ) {
						SwMakeTimeArray = true;
						break;
					} else if( 0 == strcmp(parg+1,"-double") ) {
						SwDoubleTimeArray = true;
						break;
					} else if( 0 == strcmp(parg+1,"-signal-name") ) {
						SwSignalName = true;
						break;
					}
					else {
						swChkOpt1Len = true;
					}
					break;
				case 'V': // 2010-03-10
					if( 3 > strlen(parg+1) ){
						swChkOpt1Len = false;
						switch(*(parg+2))
						{
							case '\0':
							case 'F':
							case '4':
								SwVoltsFlt = true;
								break;
							case 'V':
							case 'D':
							case '8':
								SwVoltsDbl = true;
								break;
							default:
								swChkOpt1Len = true;
								break;
						}
					}
					break;
				default:
					Usage( argv[0] );
					return -1;
					break;
			}
			if( swChkOpt1Len && 1 != strlen(parg+1) || (SwVoltsDbl && SwVoltsFlt)) {
				Usage( argv[0] );
				return -1;
			}
		}
		else {
			size_t len = strlen(parg);
#if defined(WIN32) || defined(WIN64)
			if( parg[len-1] == '\\' ) {
#else
			if( parg[len-1] == '/' ) {
#endif
				path_for_orgname = new char[len+1];
				strcpy(path_for_orgname, parg);
				sprintf(Filename ,"%s%s", parg, Diag_Name);

			}
			else {
				strcpy(Filename, parg);
			}
		}
	}
	ret = 0;
	int chs_cnt = 0;
	short* st_chs = NULL;
	short* en_chs = NULL;
	for(;;) {
		if( !SwDirectAccess ) retrieve_SetTestPort(PortNumber);

		if( !SwDirectAccess && NULL == IndexServer ) {
			error_mes = "Environment variable INDEXSERVERNAME unknown!";
			ret = -1;
			break;
		}
		chs_cnt = ChannelDecode(Channel, NULL, NULL);
		if( 0 > chs_cnt ) {
			error_mes = "<Error: Illegal channel number!>";
			ret = -2;
			break;
		}
		else if( 0 == chs_cnt ) {
			chs_cnt = 1;
			st_chs = new short[chs_cnt];
 			en_chs = new short[chs_cnt];
			signal_name = Channel;
			st_chs[0] = 1;
			en_chs[0] = 1;
		}
		else {
			st_chs = new short[chs_cnt];
 			en_chs = new short[chs_cnt];
      		chs_cnt = ChannelDecode(Channel, st_chs, en_chs);
		}
		
		if( NULL == Frame ) {
			sframe = eframe = 0;
		}
		else if( NULL != strchr(Frame,':') )	{
			if(sscanf(Frame,"%d:%d",&sframe,&eframe)!=2 || 0 > sframe || 0 > eframe || eframe < sframe )	{
				error_mes = "<Error: Illegal frame number!>";
				ret = -5;
				break;
			}
		}
		else	{
			if(sscanf(Frame,"%d",&sframe)!=1 || 0 > sframe)	{
				error_mes = "<Error: Illegal frame number!>";
				ret = -6;
				break;
			}
			eframe = sframe;
		}

		if( SwOriginalName && ( SwVoltsDbl || SwVoltsFlt ) ) {
			error_mes = "<Error: Cannot use optional -V group and --original_name at the same time.>";
			ret = -7;
		}
		break;
	}
	
	if( 0 != ret ) {
		if( SwVerbose )	printf("%s\n" ,error_mes);
		return ret;
	}

	retrieveResetIndexServerName();

	if( SwOriginalName && NULL == path_for_orgname ) {
		path_for_orgname = new char[1];
		path_for_orgname[0] = '\0';
	}

	startIndexNoClosing();
	int rd = -1;
	if( !SwWait ) TimeoutSec = 0;
	if ( SwDirectAccess) {
		func_name = "retrieveOpenDirectWait";
		rd=retrieveOpenDirectWait(Diag_Name, Server, Path, AShot_Number, ASub_Shot_Number, PortNumber,TimeoutSec);
	}
	else {
		func_name = "retrieveOpenWait";
		rd=retrieveOpenWait(Diag_Name,IndexServer,AShot_Number,ASub_Shot_Number,TimeoutSec);
	}
	if(0 > rd ) {
		if( SwVerbose )	printf("Error %d [ %s ] from DLL(%s)!\n" ,
						rd ,retrieveErrorMessage(rd),func_name);
		endIndexNoClosing();
		return rd;
	}

	unsigned int n_channel;
	short	year ,month ,day ,hour ,min ,sec;
	char	management[32];
	char	comment[128];
	int		comment_size=128;
	char	cserver[32];

	for(;;) {
		func_name = "retrieveShotInfo";
		ret = retrieveShotInfo( rd ,&n_channel ,&year ,&month ,&day ,&hour ,&min ,&sec ,
								   management ,comment ,comment_size ,cserver ); 
		if(0 > ret ) break;
		int real_channel = n_channel;
		if( SwDirectAccess && 0 == n_channel )	{
			n_channel = 9999;
		}
		if( SwArcRealShot ) {
			func_name = "retrieveRealShotNumber";
			ret = retrieveRealShotNumber( rd, &arc_shot, &arc_sub );
			if( 0 != ret ) break;
		}
		else {
			arc_shot= AShot_Number;
			arc_sub	= ASub_Shot_Number;
		}
		if( SwVerbose ) {
			printf("---Arc Shot Infomation---\n");
			printf("Index Server Name  = %s\n",IndexServer);
			printf("Trans Server Name  = %s\n",cserver);
			printf("Shot Number        = %10u\n",arc_shot);
			printf("Sub Shot Number    = %10hu\n",arc_sub);
			printf("Number Of Channel  = %10u\n",real_channel);
			printf("Channel            = %s\n", Channel);
	//		printf("Channel            = %3d:%3d\n",sch,ech);
			if( 0 != sframe )
				printf("Frame              = %3d:%3d\n",sframe,eframe);
			printf("DATE               = %4d-%2d-%2d\n",year,month,day);
			printf("TIME               = %4d:%2d:%2d\n",hour,min,sec);
		}
		break;
	}
	if(0 > ret ) {
		if( SwVerbose )	printf("Error %d [ %s ] from DLL(%s)!\n" ,
						ret, retrieveErrorMessage(ret), func_name);
		retrieveClose( rd );
		endIndexNoClosing();
		return -1;
	}

	for(int i=0;i<chs_cnt;i++) {
		if( 0 == en_chs[i] )	en_chs[i] = n_channel;
		if( 0 == st_chs[i] )	st_chs[i] = 1;
		if( n_channel < (unsigned int)st_chs[i] || n_channel < (unsigned int)en_chs[i]) {
			if( SwVerbose ) printf("<Error: Illegal start or end channel number!>\n");
			retrieveClose( rd );
			endIndexNoClosing();
			return -1;
		}
	}
	{
		char tmpwork[256];
		strcpy(tmpwork, Filename);
		sprintf(Filename, "%s-%u-%hu", tmpwork, arc_shot, arc_sub);
	}
	unsigned int data_length;
	unsigned int comp_length;
	unsigned short param_count;
	short data_type;
	char image_type[32];
	unsigned short value_len;
	int is_nframe;
	const char *data_type_name;

	exit_status = 0;
	int	needDTSlink = 1;
	if(arc_shot < 56221) {
		needDTSlink = 0;
	}
	for(int loop=0; loop<chs_cnt; loop++ ) {
		int st_ch = st_chs[loop];
		int en_ch = en_chs[loop];
		int inc =1;
		chno=st_ch;
		if( st_ch > en_ch ) {
			inc = -1;
			st_ch = en_ch;
			en_ch = chno;
		}
		char *filepath	= new char[256];
	for(int cloop=st_ch ;cloop<=en_ch ;cloop++,chno += inc)	{
		func_name = "retrieveChInfo";
		if( NULL == signal_name ) {
			ret=retrieveChInfo( rd ,chno ,
					&data_length ,&comp_length ,&param_count ,&data_type ,image_type ,&value_len ,&is_nframe ,
					management ,comment ,comment_size );
			if( -20003 == ret && 9999 == n_channel ) break;
			if(0 > ret) {
				if( SwVerbose ) printf("Error %d [ %s ] from DLL(%s CH=%d)!\n" ,
								ret, retrieveErrorMessage(ret), func_name, chno);
				exit_status = ret;
				continue;
			}
		}
		else {
			ret=retrieveChInfo2( rd ,signal_name ,
					&data_length ,&comp_length ,&param_count ,&data_type ,image_type ,&value_len ,&is_nframe ,
					management ,comment ,comment_size, &chno );
			if(0 > ret) {
				if( SwVerbose ) printf("Error %d [ %s ] from DLL(%s SignalName=%s)!\n" ,
								ret, retrieveErrorMessage(ret), func_name, signal_name);
				exit_status = ret;
				continue;
			}
		}

		
		switch(data_type)
		{
			case RAW_DATA:
				data_type_name="RAW";
				break;
			case ANA_DATA:
				data_type_name="ANA";
				break;
			default:
				data_type_name="???";
				break;
		}
		unsigned int volts_buf_size	= 0;
		if( SwVoltsDbl || SwVoltsFlt ) {
			if(       9 > value_len ) {	volts_buf_size = data_length;	}
			else if( 17 > value_len ) {	volts_buf_size = data_length/2;	}
			else if( 33 > value_len ) {	volts_buf_size = data_length/4;	}
		}
		if( SwVerbose ) {
			printf("---Channel Information---\n");
			printf("Channel Number     = %d\n"	,chno);
			printf("Data Type          = %s\n"	,data_type_name);
			printf("Data Size          = %u\n"	,data_length);
			printf("Comp Size          = %u\n"	,comp_length);
			if( 0 != is_nframe )
				printf("Number of Frame    = %d\n"	,is_nframe);
			printf("Comment            = %s\n"	,comment);
		}
		
		orgname[0]='\0';
		ret = 0;
		int	j = 0;
		char**	Key = NULL;
		char**	Val = NULL;
		int*	Typ = NULL;
		char*	module_type = "unkown";
		if( 0 < param_count ) {
			Key=new char *[param_count];
			Val=new char *[param_count];
			Typ=new int[param_count];
			for(j=0;j<param_count;j++)	{
				Key[j]=new char[128];
				Val[j]=new char[128];
			}
			for(;;) {
				func_name = "retrieveChParams";
				ret =  retrieveChParams(rd, chno, Key ,Val , Typ );
				if( 0 > ret ) {
					if( SwVerbose ) printf("Error %d [ %s ] from DLL(%s CH=%d)!\n" ,
									ret, retrieveErrorMessage(ret), func_name, chno);
					break;
				}
				ret = 0;
				for(j=0;j<param_count;j++) {
					if( 0 == strcmp("ModuleType",Key[j]) ) {
						module_type = Val[j];
						break;
					}
				}
				if( SwOriginalName ) {
					for(j=0;j<param_count;j++) {
						if( 0 == strcmp("FileName",Key[j]) ) {
							strcpy(orgname,Val[j]);
							break;
						}
					}
				}
				break;
			}
		}
		FILE *fp		= NULL;

		filepath[0] = '\0';
		if( SwOriginalName ) {
			if( '\0' != orgname[0] ) {
				sprintf(filepath,"%s%s.prm",path_for_orgname,orgname);
			}
			else {
				if( SwVerbose ) printf("Warning undefined original file name.\n");
			}
		}
		if( '\0' == filepath[0] ) {
			if( NULL != signal_name && SwSignalName ) {
				sprintf(filepath,"%s-%s.prm", Filename, signal_name);
			}
			else {
				sprintf(filepath,"%s-%d.prm", Filename, chno);
			}
		}
		if( SwVerbose ) printf("-----Save to parameter file      [%s ]...",filepath);
		if((fp=fopen(filepath,"w"))==NULL)
		{
			if( SwVerbose ) printf("Error %d Can't open the file[%s]",errno,filepath);
			retrieveClose( rd );
			exit_status = ERROR_FILE_OPEN;
			break;		// -- exit for loop ch.
		}

		for(j=0;j<param_count;j++) {
			if( 1 == (int)Typ[j] && retrieve_check_string(Val[j]) ) {
				fprintf(fp,"%s,%s,\"%s\",%d\n" ,module_type ,Key[j] ,Val[j] ,(int)Typ[j]);
			}
			else {
				fprintf(fp,"%s,%s,%s,%d\n" ,module_type ,Key[j] ,Val[j] ,(int)Typ[j]);
			}
		}
		fclose(fp);
		if( -1 < ret ) {
			if( SwVerbose ) printf("OK\n");
		}
		else {
			exit_status = ret;
		}
		if( SwMakeTimeArray ) {
			unsigned int	retShot;
			unsigned short	retSubShot;
			char		DTSsource[64], DTShostID[64], DTSmoduleID[64], CLKsource[64], CLKhostID[64], CLKmoduleID[64];
			short		DTStriggerChannel, DTSclockChannel;
			char		strDTStriggerChannel[64], strDTSclockChannel[64], ExtOrInt[64], InternalClock[64], SamplingInterval[64], PreSampling[64];
			int			DTSuserDefine, DTStimeArraySize;
			int 		iPreSamplings, iSamplingInterval, ClkDTSid;
			char		ClkDTSmodule[64];
			int64_t		LastCount;

			void* timeArray=NULL;
			void* tmp_clock=NULL;
			void* tmp_start=NULL;

			int bytes_per_data;
			if(SwDoubleTimeArray) {
				tmp_clock = new double;
				tmp_start = new double;
				bytes_per_data = sizeof(double);
			}
			else {
				tmp_clock = new float;
				tmp_start = new float;
				bytes_per_data = sizeof(float);
			}
			ret = 0;
			if(arc_shot < 56221) {
				ret=retrieveGetDTSinfox(Diag_Name, arc_shot, arc_sub, chno, chno, IndexServer,
									DTSsource, DTShostID, DTSmoduleID,
									&DTStriggerChannel, &DTSclockChannel, &DTSuserDefine, &DTStimeArraySize,
									ExtOrInt, &iPreSamplings, &iSamplingInterval, &ClkDTSid, ClkDTSmodule);

				if(ret == 0) {
					strcpy(CLKsource, DTSsource);
					strcpy(CLKhostID, DTShostID);
					strcpy(CLKmoduleID, DTSmoduleID);
					sprintf(strDTStriggerChannel,"%d",DTStriggerChannel);
					sprintf(strDTSclockChannel,"%d",DTSclockChannel);
				}

			}
			if( 0 == ret ) {
				//	パラメータからDTSリンク情報を取り出す
				ret = GetDTSinfoFromChParameters(param_count, Key, Val, needDTSlink, &retShot, &retSubShot,
									DTSsource, DTShostID, DTSmoduleID, strDTStriggerChannel,  
									CLKsource, CLKhostID, CLKmoduleID, strDTSclockChannel,
									&DTSuserDefine, &DTStimeArraySize, ExtOrInt,
									InternalClock, SamplingInterval, PreSampling, &LastCount);

				//	DTSリンク情報から時間軸配列を作成
				if(0 == ret && DTStimeArraySize > 0) {
					if(SwDoubleTimeArray)	timeArray = new double[DTStimeArraySize];
					else					timeArray = new float[DTStimeArraySize];

					ret=retrieveGetDTSdatax2(IndexServer, 
									DTSsource, DTShostID, DTSmoduleID, strDTStriggerChannel, 
									CLKsource, CLKhostID, CLKmoduleID, strDTSclockChannel,
									&DTSuserDefine, DTStimeArraySize, retShot, retSubShot,
									ExtOrInt, InternalClock, SamplingInterval, PreSampling, LastCount, 1,
									SwDoubleTimeArray, timeArray, tmp_clock, tmp_start);
					if(ret) {
						if(bSilent == 0) printf("Error from DLL( retrieveGetDTSdatax2 )! %d\n",ret);
					}
				}
			}

			for( ;; ) {
				if( 0 != ret ) break;	// --exit for loop

				if( NULL != signal_name && SwSignalName ) {
					sprintf(filepath,"%s-%s.time", Filename, signal_name);
				}
				else {
					sprintf(filepath,"%s-%d.time", Filename, chno);
				}

				if((fp=fopen(filepath,"wb")) == NULL) {
					if(bSilent == 0)	printf("Error!! : Can't open the file[%s]\n", filepath);
					ret = ERROR_FILE_OPEN;
					break;	// --exit for loop ch
				}
				if(bSilent == 0)		printf("-----Save to time data file      [%s]...", filepath);

				fwrite(timeArray, bytes_per_data, DTStimeArraySize, fp);
				fclose(fp);

				if(bSilent == 0)		printf("OK\n");

				if( NULL != signal_name && SwSignalName ) {
					sprintf(filepath,"%s-%s.tprm", Filename, signal_name);
				}
				else {
					sprintf(filepath,"%s-%d.tprm", Filename, chno);
				}
			
				if((fp=fopen(filepath,"w"))==NULL)	{
					if(bSilent == 0)	printf("Error!! : Can't open the file[%s]\n", filepath);
					ret = ERROR_FILE_OPEN;
					break;	// --exit for loop ch
				}
				if(bSilent == 0)		printf("-----Save to time parameter file [%s]...", filepath);

				fprintf(fp,"DTSsource,%s\n", DTSsource);
				fprintf(fp,"DTShostID,%s\n", DTShostID);
				fprintf(fp,"DTSmoduleID,%s\n", DTSmoduleID);
				if( arc_shot < 56221) {
					fprintf(fp,"TriggerCh,%d\n", DTStriggerChannel);
					fprintf(fp,"ClockCh,%d\n", DTSclockChannel);
				} else {
					fprintf(fp,"TriggerCh,%s\n", strDTStriggerChannel);
					fprintf(fp,"CLKsource,%s\n", CLKsource);
					fprintf(fp,"CLKhostID,%s\n", CLKhostID);
					fprintf(fp,"CLKmoduleID,%s\n", CLKmoduleID);
					fprintf(fp,"ClockCh,%s\n", strDTSclockChannel);
				}
				fprintf(fp,"ClockSource,%s\n", ExtOrInt);
				fprintf(fp,"TimeArraySize,%d\n", DTStimeArraySize);
				if(SwDoubleTimeArray){
					fprintf(fp,"ArrayDataType,double\n");
					fprintf(fp,"ClockCycle,%fsec\n", *(double*)tmp_clock );
					fprintf(fp,"StartTiming,%fsec\n", *(double*)tmp_start);
				}else{
					fprintf(fp,"ArrayDataType,float\n");
					fprintf(fp,"ClockCycle,%fsec\n", *(float*)tmp_clock);
					fprintf(fp,"StartTiming,%fsec\n", *(float*)tmp_start);
				}

				fclose(fp);
				if(bSilent == 0)		printf("OK\n");
				break;	// --exit for loop write time data
			}	// --end for loop write time data
			if(SwDoubleTimeArray) {
				delete (double*)tmp_clock;
				delete (double*)tmp_start;
				if(timeArray) delete[] (double*)timeArray;
			}
			else {
				delete (float*)tmp_clock;
				delete (float*)tmp_start;
				if(timeArray) delete[] (float*)timeArray;
			}

		}	// -- end if SwMakeTimeArray
		for(j=0;j<param_count;j++) {
			delete[] Key[j];
			delete[] Val[j];
		}
		if( NULL != Typ )	delete[] Typ;
		if( NULL != Val )	delete[] Val;
		if( NULL != Key )	delete[] Key;

		if( 0 != ret ) {
			exit_status = ret;
			break;	// -- exit for loop ch.
		}

		if( 0 == is_nframe && 0 == sframe ) {
			filepath[0] = '\0';
			if( SwOriginalName ) {
				if( '\0' != orgname[0] ) {
					sprintf(filepath,"%s%s",path_for_orgname,orgname);
				}
				else {
					if( SwVerbose ) printf("Warning undefined original file name.\n");
				}
			}
			if( '\0' == filepath[0] ) {
				if( NULL != signal_name && SwSignalName ) {
					sprintf(filepath,"%s-%s.dat", Filename, signal_name);
				}
				else {
					sprintf(filepath,"%s-%d.dat", Filename, chno);
				}
			}
			if( SwVerbose ) printf("-----Save to data file           [%s ]...", filepath);
			if((fp=fopen(filepath,"wb"))==NULL)		{
				if( SwVerbose ) printf("Error %d Can't open the file[%s]\n",errno,filepath);
				retrieveClose( rd );
				exit_status = ERROR_FILE_OPEN;
				break;	// -- exit for loop ch
			}
			ret = 0;
			if(0 < data_length)		{
				for(;;) {
					if( SwVoltsDbl ) {
						if( 0== volts_buf_size ) {
							if( SwVerbose ) 
								printf("Error No-supported Resoution(bit=%d) (CH=%d)!\n" ,value_len,chno);
							ret = -1;
							break;	// --exit for loop write-data
						}
						unsigned int ret_length;
						double* buf=new double[volts_buf_size];
						func_name = "retrieveChVoltsDbl";
						ret=retrieveChVoltsDbl(rd, chno, buf, volts_buf_size, &ret_length);
						if(0 > ret) {
							if( SwVerbose ) printf("Error %d [ %s ] from DLL(%s CH=%d)!\n" ,
										ret, retrieveErrorMessage(ret), func_name, chno);
						}
						else {
							fwrite(buf,sizeof(double),ret_length,fp);
						}
						delete[] buf;
					}
					else if( SwVoltsFlt ) {
						if( 0== volts_buf_size ) {
							if( SwVerbose ) 
								printf("Error No-supported Resoution(bit=%d) (CH=%d)!\n" ,value_len,chno);
							ret = -1;
							break;	// --exit for loop write-data
						}
						unsigned int ret_length;
						float* buf=new float[volts_buf_size];
						func_name = "retrieveChVolts";
						ret=retrieveChVolts(rd, chno, buf, volts_buf_size, &ret_length);
						if(0 > ret) {
							if( SwVerbose ) printf("Error %d [ %s ] from DLL(%s CH=%d)!\n" ,
										ret, retrieveErrorMessage(ret), func_name, chno);
						}
						else {
							fwrite(buf,sizeof(float),ret_length,fp);
						}
						delete[] buf;
					}
					else {
						unsigned int ret_length;
						char* data=new char[data_length];
						func_name = "retrieveChData";
						ret=retrieveChData(rd,chno,data,data_length,&ret_length);
						if(0 > ret) {
							if( SwVerbose ) printf("Error %d [ %s ] from DLL(%s CH=%d)!\n" ,
										ret, retrieveErrorMessage(ret), func_name, chno);
						}
						else {
							fwrite(data,1,ret_length,fp);
						}
						delete[] data;
					}
					break;	// --exit for loop write-data
				}	// --end for loop write-data
			}	// --end if
			fclose(fp);
			if( -1 < ret ) {
				if( SwVerbose ) printf("OK\n");
			}
			else {
				exit_status = ret;
			}
		}
		else {
			int start_no = sframe;
			int end_no	= eframe;
			if( 0 == start_no ) start_no = 1;
			if( 0 == end_no || is_nframe < end_no ) end_no = is_nframe;
			unsigned int	data_length;
			unsigned int	frame_x;
			unsigned int	frame_y;
			unsigned int	ret_length;
			char*			frame_data = NULL;
			for( int frno=start_no; frno<=end_no ; frno++ ) {
				func_name = "retrieveFrameInfo";
				ret = retrieveFrameInfo( rd ,chno ,frno ,&data_length ,image_type ,&frame_x ,&frame_y ); 
				if(0 > ret) {
					if( SwVerbose ) 
						printf("Error %d [ %s ] from DLL(%s CH=%d FRAME=%d)!\n" ,
							ret, retrieveErrorMessage(ret), func_name, chno, frno);
					exit_status = ret;
					continue;
				}
				if( NULL != signal_name && SwSignalName ) {
					sprintf(filepath,"%s-%s-%d.dat", Filename, signal_name, frno);
				}
				else {		
					sprintf(filepath,"%s-%d-%d.dat",Filename, chno, frno);
				}
				if((fp=fopen(filepath,"wb"))==NULL)		{
					if( SwVerbose ) printf("Error %d Can't open the file[%s]!\n",errno, filepath);
					exit_status = ERROR_FILE_OPEN;
					break;	// --exit loop frame
				}
				if( SwVerbose ) printf("-----Save to data file           [%s ]...", filepath);
				ret = 0;
				if(0 < data_length)		{
					frame_data=new char[data_length];
					func_name = "retrieveFrameData";
					ret=retrieveFrameData(rd,chno,frno,frame_data,data_length,&ret_length);
					if(0 > ret) {
						if( SwVerbose ) 
							printf("Error %d [ %s ] from DLL(%s CH=%d FRAME=%d)!\n" ,
									ret, retrieveErrorMessage(ret), func_name, chno, frno);
					}
					else {
							fwrite(frame_data,1,ret_length,fp);
					}
					delete[] frame_data;
				}	
				fclose(fp);
				if( -1 < ret ) {
					if( SwVerbose ) printf("OK\n");
				}
				else {
					exit_status = ret;
				}
				if( 0 != exit_status ) break;	// --exit loop frame
			}	// -- end for loop frame
		}	// -- end if 
		if( 0 != exit_status ) break;	// --exit loop ch
	}	// -- end for loop ch
		delete[] filepath;
		if( 0 != exit_status ) break;	// --exit loop chs_cnt
	}	// -- end for loop chs_cnt
	if( 0 < rd ) retrieveClose( rd );
	endIndexNoClosing();
	return exit_status;
}
