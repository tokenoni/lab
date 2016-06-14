/**
 *
 * @date 2012-07-13
 * @note SHOT#56220以前はDTSリンク情報がチャネルパラメータに保存されていなかった。
 *       そのため、dtsinfoテーブルより取得する。　dtsinfoテーブルには、その他の
 *       情報(サンプリング数等)もふくまれているが、必要十分でないため、利用しない。
 *       SHOT#56220以前はDTSホストはトリガとクロック別には定義できなかった。
 *       つまり　トリガホスト=クロックホストである。
 */
#if defined(WIN32) || defined(WIN64)
 #pragma warning(disable:4786)
  #include <windows.h>
  #include <string.h>
  typedef unsigned int uint32_t;
  typedef unsigned short uint16_t;
  typedef int int32_t;
  typedef short int16_t;
#else
 #include <stdio.h>
 #include <string>
 #define 	_MAX_PATH	256
#endif // WIN32

#include "stdio.h"
#include "retrieve.h"

#define	PATH_SEP_CHAR	NAME_SEPARATOR 

int	ret_t_error_print( int ret , const char* opt24, const char* opt25);

#define	startIndexNoClosing()	IndexNoClosing(1)
#define	endIndexNoClosing()		IndexNoClosing(0)

int main(int argc, char* argv[])
{
	if(argc < 5 || argc > 10)
	{
		const char *exec_name=strrchr(argv[0], NAME_SEPARATOR);
		if( NULL == exec_name ) exec_name = argv[0];
		else					exec_name++;
		printf(" Version %s, Link Library Version %s\n\n",RETRIEVE_T_VERSION,retrieveVersion());
		printf("<< The call must have between 4 - 9 parameters! >>\n");
		printf("	(1) diag name (Magnetics,Langmuir,.., etc.)\n");
		printf("	    or DTS name (DTS:133.75.175.100,DTS:VME01:DMOD1,.., etc.)\n");
		printf("	(2) shot number\n");
		printf("	(3) sub shot number\n");
//		printf("	(4) channel number ( 0 -> all ch, start_ch:end_ch -> ranged ch )\n");
		printf("	(4) channel number(s) (e.g. 1,3,5-8  1,3,5:8 0=all )\n");
		printf("	         or a signal tag name (e.g. sig-tag1 )\n");

		printf("	(5) [filename or pathname] option\n");
		printf("	    ( Create file name to \"filename+'-'+shotNumber+'-'+subShotNumber+'-'+channelNumber+.time\" )\n");
		printf("	    ( Create file name to \"pathname+DiagName+'-'+shotNumber+'-'+subShotNumber+'-'+channelNumber+.time\" )\n");
		printf("\n");
		printf("== Options ==\n");
		printf("	[-S | --silent]\n");
		printf("        No working messages output.\n");
		printf("	[-D | --double]\n");
		printf("        Make Time Array in double type.\n");
		printf("	[-w | --wait [-t timeout (sec)] ]\n");
		printf("		Data waiting, i.e. blocking, mode. Default timeout: 180.\n");
		printf("The usage of %s\n",exec_name);
		printf("> %s DiagName shotNumber subShotNumber channelNumber [filename] [options...]\n",exec_name);
		return -1;
	}
	char* IndexServer=getenv("INDEXSERVERNAME");
	const char* Diag_Name=argv[1];
	char Server[32]={'\0'};
	char* workDiagName;
	char* filePrefix;
	filePrefix = new char[strlen(Diag_Name)+1];
	strcpy(filePrefix, Diag_Name);
	workDiagName = new char[strlen(Diag_Name)+1];
	strcpy(workDiagName, Diag_Name);

	char* indicatedFileName = NULL;
	uint32_t AShot_Number=atoi(argv[2]);
	uint16_t ASub_Shot_Number=atoi(argv[3]);
	char* Channel=argv[4];
	short sch;
	short ech;
	int ret;
	char* parg;
	short bSilent, bDouble, bFilenameIndicated, bWait;
	bSilent = bDouble = bFilenameIndicated = bWait = 0;
	int		TimeoutSec	=180;

	if(argc > 5) {
		for(int i = 5; i < argc; i++)
		{
			parg = argv[i];

			if(*parg=='-') {
				switch(*(parg+1))
				{
					case 'S':
						bSilent = 1;
						break;
					case 'D':
						bDouble = 1;
						break;
					case 'w':
						bWait = 1;
						break;
					case 't':
						i++;
						if( i < argc ) TimeoutSec=atoi(argv[i]);
						break;
					case '-': // 2009-09-14
						if( 0 == strcmp(parg+1,"-silent") ) {
							bSilent = 1;
							break;
						}
						else if( 0 == strcmp(parg+1,"-double") ) {
							bDouble = 1;
							break;
						}
						else if( 0 == strcmp(parg+1,"-wait") ) {
							bWait = 1;
							break;
						}
						break;
					case 'X': // for TEST
						retrieveVersion();
						break;

				}
			}
			else {
				size_t len = strlen(parg);
				indicatedFileName= new char[len+1];
				strcpy(indicatedFileName, parg);
				if(filePrefix) delete [] filePrefix; 
				bFilenameIndicated = 1;

				if( parg[len-1] == PATH_SEP_CHAR ) {
					len += strlen(Diag_Name);
					filePrefix = new char[len+1];
					strcpy(filePrefix, parg);
					strcat(filePrefix, Diag_Name);
				} else {
					filePrefix = new char[len+1];
					strcpy(filePrefix, parg);
				}
			}
		}
	}
	int chs_cnt = 0;
	short* st_chs = NULL;
	short* en_chs = NULL;
	const char* signal_name =NULL;
	bool	bDTS = FALSE;

	char	*DTShost_name=NULL, *DTSmodule_name=NULL;

	startIndexNoClosing();
	ret = 0;
	for(;;) {
		if(!IndexServer)	{
			if(bSilent == 0)		printf("Environment variable INDEXSERVERNAME unknown!\n");
			ret = INDEX_ENV_UNDEFINED;
			break;
		}
	
		chs_cnt = ChannelDecode(Channel, NULL, NULL);
		if( 0 > chs_cnt ) {
			if(bSilent == 0)		printf("Illegal channel number!\n");
			ret = WRONG_PARAMETER;
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

		char	*str[3], *head, *next;
		int		str_num;

		head = workDiagName;
		next = strtok(head, ":");
		for(str_num=0 ; next != NULL && str_num < 3; str_num++){
			str[str_num] = next;
			next = strtok(NULL, ":");
		}

		if( 0 == strcmp(str[0], "DTS") ){
			if(str_num == 2){
				bDTS = TRUE;
				DTShost_name = str[1];
				DTSmodule_name = "NONAME";
			}
			else if(str_num == 3){
				bDTS = TRUE;
				DTShost_name = str[1];
				DTSmodule_name = str[2];
			}
		}
		if( NULL != signal_name && ( bDTS || AShot_Number < 56221)  ) {
			if(bSilent == 0)		printf("Signal name is not supported!\n");
			ret = WRONG_PARAMETER;
			break;
		}
		break;
	}
	if( 0 != ret ) {
		if(st_chs)				delete[] st_chs;
		if(en_chs)				delete[] en_chs;
		if(workDiagName)		delete[] workDiagName;
		if(filePrefix)			delete[] filePrefix;
		if(indicatedFileName)	delete[] indicatedFileName;
		endIndexNoClosing();
		return ret;
	}
	uint32_t	retShot;
	uint16_t	retSubShot;
	short		retParameterCount;

	char		DTSsource[64], DTShostID[64], DTSmoduleID[64], CLKsource[64], CLKhostID[64], CLKmoduleID[64];
	short		DTStriggerChannel, DTSclockChannel;
	char		strDTStriggerChannel[64], strDTSclockChannel[64], ExtOrInt[64], InternalClock[64], SamplingInterval[64], PreSampling[64];
	int32_t		DTSuserDefine, DTStimeArraySize;
	int 		iPreSamplings, iSamplingInterval, ClkDTSid;
	char		ClkDTSmodule[128];

	if( bDTS ) {
		//------------------------------- start DTS
		char	tmp_module_name[32];
		int	s=st_chs[0];
		int	e=en_chs[0];
		for(int i=0;i<chs_cnt;i++) {
			if( 0 == en_chs[i] )	{
				s = 1;
				e = 0;
				break;
			}
			if( s > st_chs[i] )	s = st_chs[i];
			if( e < en_chs[i] )	e = en_chs[i];
		}
		// Access DTSTBL db-Open Close
		ret=retrieveGetDTSParametersCount(IndexServer, DTShost_name, DTSmodule_name, AShot_Number, ASub_Shot_Number,
										&retParameterCount, &s, &e, tmp_module_name);

		if(ret) {
			if(filePrefix)			delete[]	filePrefix;
			if(indicatedFileName)	delete[]	indicatedFileName;

			if( 0 ==  bSilent ) {
				if( WRONG_PARAMETER == ret ) {
					printf("Channel parameter is illegal.\n");
				}
				else {
					printf("Error from DLL( retrieveGetDTSParametersCount )!\n");
					ret_t_error_print( ret, IndexServer, Server);
				}
			}
			endIndexNoClosing();
			return ret;
		}
		DTSmodule_name =tmp_module_name;
		for(int i=0;i<chs_cnt;i++) {
			if( 0 == en_chs[i] )	en_chs[i] = e;
			if( 0 == st_chs[i] )	st_chs[i] = 1;
		}

		char**	p_name	=NULL;
		char**	p_val	=NULL;
		int*	p_type	=NULL;

		if(bSilent == 0) {
			printf("---DTS Infomation---\n");
			printf("Shot Number        = %10d\n",AShot_Number);
			printf("Sub Shot Number    = %10d\n",ASub_Shot_Number);
			printf("DTS Source         = %s\n",DTShost_name);
			printf("Channel Range      = [ %3d - %3d ]\n", s, e);
		}

		char *file_name=new char[_MAX_PATH];
		p_name = new char*[retParameterCount];
		p_val = new char*[retParameterCount];
		p_type = new int[retParameterCount];
		for(int i = 0; i < retParameterCount; i++) {
			p_name[i] = new char[64];
			p_val[i] = new char[128];
		}
		FILE *fp;
		for(int l=0;l<chs_cnt;l++) {
			sch = st_chs[l];
			ech = en_chs[l];
		for(short t_ch = sch; t_ch <= ech; t_ch++) {
	
			for(int i = 0; i < retParameterCount; i++) {
				*p_name[i] = '\0';
				*p_val[i] = '\0';
			}
			fp = NULL;
			file_name[0]= '\0';
		
			ret = retrieveGetDTSParameters(IndexServer, DTShost_name, DTSmodule_name, AShot_Number, ASub_Shot_Number, t_ch,
										retParameterCount, p_name, p_val, p_type);

			if(ret) {
				if( 0 ==  bSilent ) {
					printf("Error from DLL( retrieveGetDTSParameters )! %d\n",ret);
					ret_t_error_print( ret, IndexServer, Server);
				}
				break;
			}
		

			if(!bFilenameIndicated) {
				sprintf(file_name,"%s-%s-%d-%d-%d.tprm",
					DTShost_name, DTSmodule_name, AShot_Number, ASub_Shot_Number, t_ch);
			}
			else {
				size_t len = strlen(indicatedFileName);
				if( indicatedFileName[len-1] == PATH_SEP_CHAR ) {
					sprintf(file_name, "%s%s-%s-%d-%d-%d.tprm",
						indicatedFileName, DTShost_name, DTSmodule_name, AShot_Number, ASub_Shot_Number, t_ch);
				} else {
					sprintf(file_name, "%s-%d-%d-%d.tprm",
						filePrefix, AShot_Number, ASub_Shot_Number, t_ch);
				}
			}

			if((fp=fopen(file_name,"w"))==NULL)	{
				if(bSilent == 0)	printf("Error!! : Can't open the file[%s]\n", file_name);
				ret = ERROR_FILE_OPEN;
				break;
			}
			if(bSilent == 0) {
				printf("-----Save to parameter file[%s]...", file_name);
			}
			for(int i = 0; i < retParameterCount; i++) {
				fprintf(fp,"%s:%s,%s,%s,%d\n", DTShost_name, DTSmodule_name, p_name[i], p_val[i], p_type[i]);
			}
			fclose(fp);

			if(bSilent == 0)		printf("OK\n");
			ret = 0;

		}// --end for loop ech
		if( 0 != ret) break;
		}// --end for loop chs_cnt
		for(int i = 0; i < retParameterCount; i++) {
			delete[]	p_name[i];
			delete[]	p_val[i];
		}
		delete[]	p_name;
		delete[]	p_val;
		delete[]	p_type;	
		if(file_name)			delete[]	file_name;

	}//------------------------------- end DTS 
	else {

		ret = 0;
		if(bWait == 0)			TimeoutSec = 0;

		unsigned int last_ch = 0;
		int rd = -1;
		rd=retrieveOpenWait(Diag_Name,IndexServer,AShot_Number,ASub_Shot_Number,TimeoutSec);
		if(0 > rd ) {
			if(bSilent == 0) {
				printf("Error %d [ %s ] !\n" ,rd ,retrieveErrorMessage(rd));
			}
			ret = rd;
		}
		else {
			char	management[32], comment[128], cserver[128];
			int		comment_size=128;
			short	retYear, retMonth, retDay, retHour, retMinute, retSecond;
		
			ret = retrieveShotInfo( rd ,&last_ch ,&retYear ,&retMonth ,&retDay ,&retHour ,&retMinute ,&retSecond,
								   management ,comment ,comment_size ,cserver ); 
			if( 0 > ret ) {
				if(bSilent == 0) {
					printf("Error %d [ %s ] !\n" ,rd ,retrieveErrorMessage(rd));
				}
			}			
		}
		if(ret) {
			if(filePrefix)			delete[]	filePrefix;
			if(indicatedFileName)	delete[]	indicatedFileName;
			if( 0 < rd ) 			retrieveClose(rd);
			endIndexNoClosing();
			return ret;
		}
		for(int i=0;i<chs_cnt;i++) {
			if( 0 == st_chs[i] )	st_chs[i] = 1;
			if( 0 == en_chs[i] )	en_chs[i] = last_ch;
		}

		void* timeArray=NULL;
		void* tmp_clock=NULL;
		void* tmp_start=NULL;

		if(bSilent == 0) {
			printf("---DTS Link Infomation---\n");
			printf("Shot Number        = %10d\n",AShot_Number);
			printf("Sub Shot Number    = %10d\n",ASub_Shot_Number);
		}
		int bytes_per_data;
		if(bDouble) {
			tmp_clock = new double;
			tmp_start = new double;
			bytes_per_data = sizeof(double);
		}
		else {
			tmp_clock = new float;
			tmp_start = new float;
			bytes_per_data = sizeof(float);
		}
		char* file_name=new char[_MAX_PATH];

		short t_ch;
		int64_t	LastCount;
		const char* p_channel;
		int	type_channel;
		int real_channel;
		if( NULL == signal_name ) {
			type_channel = sizeof(t_ch);
			p_channel = (char*)&(t_ch); 
		}
		else {
			type_channel = 0;
			p_channel = signal_name;
		}
		int needDTSlink = 1;
		if(AShot_Number < 56221) {// DTS LINK INFO. from DTSinfo table.
			needDTSlink = 0;
		}
		for(int l=0;l<chs_cnt;l++) {
			sch = st_chs[l];
			ech = en_chs[l];
			t_ch = sch;
			int inc = 1;
			if( sch > ech ) {
				inc = -1;
				sch = ech;
				ech = t_ch;
			}
		for(int cloop = sch; cloop <= ech; cloop++, t_ch += inc) {
			if( NULL != timeArray ) {
				if(bDouble) 	delete[]	(double *)timeArray;
				else 			delete[]	(float *)timeArray;
				timeArray = NULL;
			}
			ret = 0;
			if(AShot_Number < 56221) {// DTS LINK INFO. from DTSinfo table.
				ret=retrieveGetDTSinfox(Diag_Name, AShot_Number, ASub_Shot_Number, t_ch, t_ch, IndexServer,
									DTSsource, DTShostID, DTSmoduleID,
									&DTStriggerChannel, &DTSclockChannel, &DTSuserDefine, &DTStimeArraySize,
									ExtOrInt, &iPreSamplings, &iSamplingInterval, &ClkDTSid, ClkDTSmodule);
				if(ret) {
					if( 0 ==  bSilent ) {
						printf("Error from DLL( retrieveGetDTSinfox )! %d\n",ret);
					}
					break;
				}
				// DTStriigerホストとDTSclockホストは同一である。
				strcpy(CLKsource, DTSsource);
				strcpy(CLKhostID, DTShostID);
				strcpy(CLKmoduleID, DTSmoduleID);
				sprintf(strDTStriggerChannel, "%d", DTStriggerChannel);
				sprintf(strDTSclockChannel, "%d", DTSclockChannel);
			} 
			//	Retrieveからパラメータ文字列を取得	パラメータからDTSリンク情報を取り出す
			ret = retrieveGetDTSinfoFromRetrieve(rd, type_channel, p_channel, needDTSlink, &real_channel,
								&retShot, &retSubShot, 
								DTSsource, DTShostID, DTSmoduleID, strDTStriggerChannel,  
								CLKsource, CLKhostID, CLKmoduleID, strDTSclockChannel,
								&DTSuserDefine, &DTStimeArraySize, ExtOrInt,
								InternalClock, SamplingInterval, PreSampling, &LastCount);
			//	DTSリンク情報から時間軸配列を作成
			if(0 != ret ) {
				if(bSilent == 0) printf("Error from DLL( retrieveGetDTSinfoFromRetrieve )! %d\n",ret);
				break;
			}

			if(bDouble) timeArray = new double[DTStimeArraySize];
			else 		timeArray = new float[DTStimeArraySize];

			ret = retrieveGetDTSdatax2(IndexServer,
								DTSsource, DTShostID, DTSmoduleID, strDTStriggerChannel, 
								CLKsource, CLKhostID, CLKmoduleID, strDTSclockChannel,
								&DTSuserDefine, DTStimeArraySize, retShot, retSubShot,
								ExtOrInt, InternalClock, SamplingInterval, PreSampling, LastCount, bSilent,
								bDouble, timeArray, tmp_clock, tmp_start);
			if(ret) {
				if(bSilent == 0) printf("Error from DLL( retrieveGetDTSdatax2 )! %d\n",ret);
				break;
			}
		
			
			FILE *fp = NULL;

			sprintf(file_name,"%s-%d-%d-%d.time", filePrefix, AShot_Number, ASub_Shot_Number, real_channel);

			if((fp=fopen(file_name,"wb")) == NULL) {
				if(bSilent == 0)	printf("Error!! : Can't open the file[%s]\n", file_name);
				ret= ERROR_FILE_OPEN;
				break;
			}
			if(bSilent == 0)		printf("-----Save to data      file[%s]...", file_name);

			fwrite(timeArray, bytes_per_data, DTStimeArraySize, fp);
			fclose(fp);

			if(bSilent == 0)		printf("OK\n");

			sprintf(file_name,"%s-%d-%d-%d.tprm", filePrefix, AShot_Number, ASub_Shot_Number, real_channel);
			
			if((fp=fopen(file_name,"w"))==NULL)	{
				if(bSilent == 0)	printf("Error!! : Can't open the file[%s]\n", file_name);
				ret = ERROR_FILE_OPEN;
				break;
			}
			if(bSilent == 0)		printf("-----Save to parameter file[%s]...", file_name);

			fprintf(fp,"DTSsource,%s\n", DTSsource);
			fprintf(fp,"DTShostID,%s\n", DTShostID);
			fprintf(fp,"DTSmoduleID,%s\n", DTSmoduleID);
			if( AShot_Number < 56221) {
				fprintf(fp,"TriggerCh,%s\n", strDTStriggerChannel);
				fprintf(fp,"ClockCh,%s\n", strDTSclockChannel);
			} else {
				fprintf(fp,"TriggerCh,%s\n", strDTStriggerChannel);
				fprintf(fp,"CLKsource,%s\n", CLKsource);
				fprintf(fp,"CLKhostID,%s\n", CLKhostID);
				fprintf(fp,"CLKmoduleID,%s\n", CLKmoduleID);
				fprintf(fp,"ClockCh,%s\n", strDTSclockChannel);
			}
			fprintf(fp,"ClockSource,%s\n", ExtOrInt);
			fprintf(fp,"TimeArraySize,%d\n", DTStimeArraySize);
			if(bDouble){
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
				
			ret = 0;
		}// --end for loop ech
		if( 0 != ret) break;
		}// --end for loop chs_cnt
		if(rd >= 0)	retrieveClose( rd );

		if(ret) {
			if( 0 ==  bSilent ) {
				ret_t_error_print( ret, IndexServer, Server);
			}
		}
		if( NULL != timeArray ) {
			if(bDouble) 	delete[]	(double *)timeArray;
			else 			delete[]	(float *)timeArray;
		}
		if(bDouble) {
			delete (double*)tmp_clock;
			delete (double*)tmp_start;
		}
		else {
			delete (float*)tmp_clock;
			delete (float*)tmp_start;
		}
		if(file_name)			delete[]	file_name;
	}

	if(st_chs)				delete[] st_chs;
	if(en_chs)				delete[] en_chs;
	if(workDiagName)		delete[] workDiagName;
	if(filePrefix)			delete[] filePrefix;
	if(indicatedFileName)	delete[] indicatedFileName;
	endIndexNoClosing();
	return ret;

}
