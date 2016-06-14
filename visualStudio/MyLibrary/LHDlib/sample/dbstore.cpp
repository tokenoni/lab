/**
 * @file  dbStore.cpp
 * @brief Implemention dbStore.exe and CdbStore Class. 
 * 
 * @author s.imazu
 * @date   2007-03
 * @version 1.0
 *
 * @date   2012-04-26
 * @version 1.1
 * @note    add option : --accept_sizes,--cndb_timeout,--cndb_retry
 * @note    set sigal_name in shot file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "dbstore.h"
#include "param_item_name.h"
#include <vector>

#if defined(__linux__) || defined(__MACH__)
 #include <strings.h>
 #include <getopt.h>
 #define _stricmp( x1 ,x2) strcasecmp( x1 ,x2)
#else
 #include "getopt_win32.h"
#endif

class CdbStore
{
	char*			mailAddress;
	char*			diagName;
	unsigned int	shotNumber;
	short			subShot;
	int			startChNumber;
	int			endChNumber;
	short			dataType;
	char*			dataPath;
	const char*	defImageType;
	bool			swAcceptSizes;
	int			cndbTimeout;
	int			cndbRetry;
public:
	CdbStore( void );
	~CdbStore( void);

	int setALL( int ac, char** av);
	
	int	run( void );
private:
	void makeChParamFileName( char* filename, int ch, int is_old);
	void makeChDataFileName( char* filename, int ch,int is_old);
	void makeFrameDataFileName( char* filename, int ch, int frame);
	int	 getParameters( char* pfname,unsigned long* p_cnt,unsigned long* p_size,unsigned char** p_list,
				int* is_frame_num,char** image_type,int *x_size,int *y_size,int *frame_size);
	int getData( char* dfname, unsigned long*data_length,unsigned char **data_data);
	
	bool getKeywordValue( const char* key, const char* value,
				int* framenum,int* framesize,int *xsize,int *ysize,int *valuelen,
				char*image_type );
	bool isExceptParam( const char* key);
	bool isExceptFrameParam( const char* key);
	bool existFrameDat( char *pfname );
};

void usage( char *av );

int main( int ac, char **av)
{
	if( ac< 7  ) {
		usage(av[0]);
		return(-1);
	}
	
	int	ret;
	CdbStore	store;
	
	ret = store.setALL( ac , av );
	if( 0 != ret ) {
		usage(av[0]);
		return(-1);
	}
	
	ret = store.run();
	
	if( 0 == ret )	printf("Success\n");
	else  			printf("Error(%d:%s)\n",ret,dbsErrorMessage(ret));
	return(ret);	
}

void usage( char *av)
{
	printf(" Version %s, Link Library Version %s, zlib %s\n\n",DBSTORE_VERSION,dbsVersion(),dbsZversion());
	printf("Usage : %s mail_addr diag shot sub_shot ch data_type [path [image_type]] [Options]\n\n",av);
	printf(
"    mail_addr : your mail address.     e.g.) LABCOM@xxxx.xxxx.xxxx\n"
"    diag      : diagnostics name.      e.g.) Magnetics\n"
"    shot      : shot number.           e.g.) 990000\n"
"    sub_shot  : sub shot number.       e.g.) 1\n"
"    ch        : orignal start channel number and colon(:) and end number.\n"
"                                       e.g.) 1:10\n"
"    data_type : data type. RAW or ANA  e.g.) RAW\n"
"    path      : data saved path. defalut is curent.    e.g.) TestData \n"
"    image_type: data binary type. INT8 ,INT16 ,INT32 ,FLT32 or FLT64 \n"
"\n"
"List of Options:\n"
"\n"
"    --accept_sizes (-z) : accept the frame data of illegal size.\n"
"    --cndb_timeout=n (-t n) : set timeout(sec) at db connection.\n"
"    --cndb_retry=m   (-r m) : set retry count at db conncection.\n"
"\n"
	);
}
/**　コンストラクタ
*
*/
CdbStore::CdbStore(void)
{
	mailAddress	= NULL;
	diagName		= NULL;
	dataPath		= NULL;
	defImageType	= NULL;
	startChNumber	= 0; 
	endChNumber	= 0;
	swAcceptSizes	= false;		// 2012-04-20
	cndbTimeout	= -1;
	cndbRetry		= -1;
}
/**　デストラクタ
*
*/
CdbStore::~CdbStore(void)
{
	if( NULL != mailAddress	) delete [] mailAddress;
	if( NULL != diagName		) delete [] diagName;
	if( NULL != dataPath		) delete [] dataPath;
}


/**　属性のSET
*
* @param[in] ac 実行時の引数の数 
* @param[in] av 実行時の引数 
*
* @retval 0 成功  
* @retval <0 エラー発生
*   
*/
static char* IMAGE_TYPES[] = 
{ "INT8","INT16","INT32","FLT32","FLT64","INT64",NULL};
int CdbStore::setALL( int ac, char** av)
{
	int len;
	struct option lngopt[] = {
			{"accept_sizes",	0, NULL, 'z'},
			{"cndb_timeout",	1, NULL, 't'},
			{"cndb_retry"  ,	1, NULL, 'r'},
			{0, 0, 0, 0},
	};
	const char* shtopt = "-z?t:r:";

	int	opt = 0;
	int	opt_idx=0;
	int arg_no =0;
	char* sep = NULL;
	int	tmp;
	
	while( (opt=getopt_long(ac, av, shtopt, lngopt, &opt_idx)) != -1 )
	{
		switch( opt ) {
			case 1 :
				arg_no++;
				len = (int)strlen( optarg);
				if( 1 > len ) return( 0 - arg_no );
				switch( arg_no ) {
					case 1:
						mailAddress = new char[len+1];
						strcpy(mailAddress, optarg);
						break;
					case 2:
						diagName = new char[len+1];
						strcpy(diagName, optarg);
						break;
					case 3:
						shotNumber = atol(optarg);
						if( 1 > shotNumber ) 						return(-3);
						break;
					case 4:
						subShot = (short)atoi(optarg);
						if( 0 > subShot ) 							return(-4);	// 2009-10-27
						break;
					case 5:
						sep = strchr(optarg, ':');
						startChNumber = atoi(optarg);
						if( 1 > startChNumber )						return(-5);
						if( NULL == sep ) {
							endChNumber = startChNumber;
						}
						else {
							endChNumber = atoi(sep+1);
							if( startChNumber > endChNumber ) 		return(-5);
						}
						break;
					case 6:
						if     ( 0 == strcmp("RAW", optarg) )	dataType =1 ;
						else if( 0 == strcmp("ANA", optarg) )	dataType =2 ;
						else										return(-6);
						break;
					case 7:
						dataPath = new char[len+1];
						strcpy(dataPath, optarg);
						break;
					case 8:
						defImageType = NULL;
						for(int i=0;;i++) {
							if( NULL == IMAGE_TYPES[i] ) break; 
							if( 0 == strcmp(optarg, IMAGE_TYPES[i]) ) {
								defImageType = IMAGE_TYPES[i];
								break;
							}
						}
						if( NULL == defImageType ) return(-8);
						break;
					default:

						break;
				}
				break;
			case 'z' :
				swAcceptSizes = true;
				break;
			case 't' :
				tmp = atoi(optarg);
				if( -1 < tmp ) cndbTimeout = tmp;
				break;
			case 'r' :
				tmp = atoi(optarg);
				if( -1 < tmp ) cndbRetry = tmp;
				break;
			case '?' :
			default :
				return -9;
				break;
		}
	}
	if( NULL == dataPath ) {
		dataPath = new char[2];
		strcpy(dataPath,".");
	}
	return 0;
}

/**　dbStore処理
*
* @retval 0　成功
* @retval <0 エラー 
*   
*/
int CdbStore::run( void )
{
	long	dbs;
	dbs = dbsOpen2( mailAddress, diagName, (unsigned int)shotNumber, (short)subShot,
			dataType, cndbTimeout, cndbRetry);
	if( 0 >dbs ) return (int)dbs;
	
	unsigned long	param_cnt;
	unsigned long	param_size;
	unsigned char*	param_list=NULL;
	int 	is_frame_num;
	char* 	image_type=NULL;
	int		frame_x;
	int		frame_y;
	int		frame_s;
	unsigned long	data_length;
	unsigned char*	data_data=NULL;
	char	ch_fname[256];
	int		ret=0;
	int		dst_ch = 1;
	int		is_old = 0;
	for( int src_ch=startChNumber ; src_ch<=endChNumber ; src_ch++,dst_ch++) {
		
		if( NULL != data_data ) {	delete [] data_data; data_data=NULL;	}
		param_cnt = param_size = is_frame_num = 0;
		for( is_old=0;is_old<2;is_old++) {
			if( NULL != param_list ) {	delete [] param_list; param_list=NULL;	}
			if( NULL != image_type ) {	delete [] image_type; image_type=NULL;	}
			makeChParamFileName( ch_fname, src_ch,is_old);
			ret = getParameters( ch_fname, &param_cnt, &param_size, &param_list,
						&is_frame_num, &image_type, &frame_x, &frame_y, &frame_s);
			if( -1 != ret ) break;  
		}
		if( -2 == ret ) { ret = PARAM_FILE_READ_ERROR;	break;}
		if( 1  == ret ) { ret = WRONG_PARAMETER_ITEM;	break;}
		ret = 0;
		if( 0 != is_frame_num ) {
			if( NULL  == image_type ) { ret = WRONG_PARAMETER_ITEM;	break;}
			ret = dbsOpenFrame2( dbs, dst_ch, frame_x,frame_y, image_type, swAcceptSizes );
			if( 0 != ret ) break;
			for(int fr=1;fr <= is_frame_num ;fr++) {
				data_length = 0;
				if( NULL != data_data ) {	delete [] data_data; data_data=NULL;	}
				makeFrameDataFileName( ch_fname, src_ch, fr);
				ret = getData( ch_fname, &data_length, &data_data);
				if( -2 == ret ) {
					ret = DATA_FILE_READ_ERROR;
					break;
				}
				else if( -1 == ret ) {		// 2012-04-16 
					ret = FILE_NOT_FOUND;
					break;
				}
				ret = dbsWriteFrame( dbs, dst_ch, fr, data_length, data_data);
				if( 0 != ret ) break;
			}
			if( 0 == ret ) { 
				ret = dbsCloseFrame( dbs, dst_ch, param_cnt, param_size, param_list);
			}
		}
		else {
			
			data_length = 0;
			for( is_old=0;is_old<2;is_old++) {
				makeChDataFileName( ch_fname, src_ch,is_old);
				ret = getData( ch_fname, &data_length, &data_data);
				if( -1 != ret ) break;
			}
			if( -2 == ret ) {
				ret = DATA_FILE_READ_ERROR;
				break;
			}
//			if( 0 < param_cnt || 0 < data_length ) {
			if( 0 < param_cnt || 0 == ret ) {
				if( NULL != image_type ) {
					ret = dbsWrite( dbs, dst_ch,
								param_cnt, param_size, param_list,
								data_length, data_data, image_type );
				}
				else {
					if( NULL != defImageType || 0 == data_length ) {
						ret = dbsWrite( dbs, dst_ch,
									param_cnt, param_size, param_list,
									data_length, data_data, (char*)defImageType );
					}
					else {
						ret = ILLEGAL_IMAGE_TYPE;
					}
				}
			}
			else {
				ret = FILE_NOT_FOUND;
			}
		}
		if( 0 != ret ) break;
		
	}
	if( NULL != param_list ) {	delete [] param_list; param_list=NULL;	}
	if( NULL != image_type ) {	delete [] image_type; image_type=NULL;	}
	if( NULL != data_data ) {	delete [] data_data; data_data=NULL;	}
	
	if( 0 == ret ) {
		ret = dbsClose(dbs);
	}
	else {
		dbsAbort(dbs);
	}
	
	return(ret);
}
/**　チャネルパラメータファイル名の作成
*
* @param[out] fname チャネルパラメータファイル名
* @param[in]  ch	チャネル番号
* @param[in] is_old 旧式ファイル名の指定　0:旧　1:新
*   
*/
void CdbStore::makeChParamFileName( char* fname, int ch, int is_old)
{
	if( 0 == is_old ) { 
		sprintf(fname,"%s/%s-%u-%hd-%d.prm",dataPath,diagName,shotNumber,subShot,ch);
	}
	else {
		sprintf(fname,"%s/%s%u-%hd-%d.prm",dataPath,diagName,shotNumber,subShot,ch);
	}
}
/**　チャネルデータファイル名の作成
*
* @param[out] fname チャネルデータファイル名
* @param[in]  ch	チャネル番号
* @param[in] is_old 旧式ファイル名の指定　0:旧　1:新
*   
*/
void CdbStore::makeChDataFileName( char* fname, int ch, int is_old)
{
	if( 0 == is_old ) { 
		sprintf(fname,"%s/%s-%u-%hd-%d.dat",dataPath,diagName,shotNumber,subShot,ch);
	}
	else {
		sprintf(fname,"%s/%s%u-%hd-%d.dat",dataPath,diagName,shotNumber,subShot,ch);
	}
}
/**　フレームデータファイル名の作成
*
* @param[out] fname フレームデータファイル名
* @param[in]  ch	チャネル番号
* @param[in]  frame フレーム番号
*   
*/
void CdbStore::makeFrameDataFileName( char* fname, int ch, int frame)
{
	sprintf(fname,"%s/%s-%u-%hd-%d-%d.dat",dataPath,diagName,shotNumber,subShot,ch,frame);
}
/**　チャンネルパラメータファイルから、パラメータを取得する。
*
* @param[in]  pfname		チャンネルパラメータファイル名
* @param[out] p_cnt			パラメータ数
* @param[out] p_size		パラメータのサイズ
* @param[out] p_list		パラメータの格納領域
* @param[out] isframenum	フレーム数（0の場合、フレーム構造でない）
* @param[out] imagetype		データタイプ
* @param[out] xsize			フレームＸサイズ
* @param[out] ysize			フレームＹサイズ
* @param[out] framesize		フレームサイズ
*
* @retval 0 成功
* @retval -1 ファイルオープンエラー
* @retval -2 ファイルリードエラー  
* 
*/
int	CdbStore::getParameters( char* pfname,
				unsigned long* p_cnt,unsigned long* p_size,unsigned char** p_list,
				int* isframenum,char** imagetype,int *xsize,int *ysize,int *framesize)
{
	FILE *fp;
	fp=fopen(pfname,"rt");
	if( NULL == fp ) return(-1);
	char	buffer[256];
	
	int		i;
	char	*item1;
	char	*item2;
	char	*item3;
	int		len;
	int		max_len=0;
	int		x_size=0,y_size=0,frame_size=0,value_len=0,frame_num=0;
	char	img_type[32];
	int		status	= -2;
	int		params	= 0;
	std::vector<char*> item_key;	
	std::vector<char*> item_val;	
	std::vector<int>   item_type;
	char	*p_key;
	char	*p_val;
	int		p_type;	
	memset(img_type,'\0',sizeof(img_type));
	for(;;) {
		if( NULL == fgets(buffer,sizeof(buffer),fp) ) {
			if( 0 != feof(fp) ) status =0; ///< normal eof;
			break;
		}
		if( 4>strlen(buffer) ) continue;
		
		item1 = strchr(buffer,',');
		if( NULL == item1 ) 	break;
		item2 = strchr(++item1,',');
		if( NULL == item2 ) 	break;
		*item2='\0';
		item3 = strchr(++item2,',');
		if( NULL == item3 ) 	break;
		*item3='\0';item3++;
		
		if( !getKeywordValue(item1,item2,&frame_num,&frame_size,&x_size,&y_size,
					&value_len,img_type) ) {
			if( isExceptParam(item1) ) continue;
		}		
		len = (int)(item2-item1);
		if( len > max_len ) max_len=len;
		p_key = new char[len];
		strcpy(p_key,item1);
		item_key.push_back(p_key);
				
		len = (int)(item3-item2);
		if( len > max_len ) max_len=len;
		p_val = new char[len];
		strcpy(p_val,item2);
		item_val.push_back(p_val);
		
		p_type = atoi(item3);
		item_type.push_back(p_type);
		
		params++;
	}
	fclose(fp);
	
	if( 0 == status ) {


		if( 0 < frame_num && existFrameDat(pfname) ) {
			*isframenum	= frame_num;
			*xsize		= x_size;
			*ysize		= y_size;
			*framesize	= frame_size;
			len = (int)strlen(img_type)+1;
			*imagetype = new char[len];
			strcpy(*imagetype,img_type);
		}
		else {
			*isframenum	= 0;
			if( 1 > strlen(img_type) ) {
				if      ( 1  > value_len )	{	;	}
				else if( 64 < value_len )	{	;	}
				else if( 32 < value_len )	{	strcpy(img_type,"INT64");	}
				else if( 16 < value_len )	{	strcpy(img_type,"INT32");	}
				else if(  8 < value_len )	{	strcpy(img_type,"INT16");	}
				else						{	strcpy(img_type,"INT8");	}
			}
			
			len = (int)strlen(img_type)+1;
			if( 1 < len ) {
				*imagetype = new char[len];
				strcpy(*imagetype,img_type);
			}
		}
	}
	if( 0 == status ) {	
		unsigned char *param_list = new unsigned char[params*3*max_len];
		*p_list =param_list;
		int num = params;
		for( i=0;i<num;i++) {
			if( 0 != *isframenum &&  isExceptFrameParam(item_key[i]) ) {
				params--;
				continue;
			}
			strcpy((char *)param_list,item_key[i]);			param_list += max_len;
			strcpy((char *)param_list,item_val[i]);			param_list += max_len;
			sprintf((char *)param_list,"%d",item_type[i]);	param_list += max_len;
		}
		*p_cnt	= params*3;
		*p_size	= max_len; 
	}
	int num =(int)item_key.size();
	for( i=0;i<num;i++) {
		delete [] item_key[i];
		delete [] item_val[i];
	}
	
	return(status);
}
bool CdbStore::existFrameDat( char *pfname )
{
	struct stat buf;
	char	dfname[256];
	char	*pt;
	strcpy( dfname, pfname);
	pt = strrchr(dfname,'.');
	if( NULL == pt ) return (false);
	strcpy(pt,"-1.dat");
	int ret = stat(dfname,&buf);
	if( 0 != ret  ) return(false);
	
	return( true ); 
}
/**　パラメータを必要な属性をバイナリ値に変換する。
*
* @param[in]  key			パラメータのキー（項目名）
* @param[in] value			パラメータの値
* @param[out] p_size		パラメータのサイズ
* @param[out] p_list		パラメータの格納領域
* @param[out] framenum		フレーム数（0の場合、フレーム構造でない）
* @param[out] framesize		フレームサイズ
* @param[out] xsize			フレームＸサイズ
* @param[out] ysize			フレームＹサイズ
* @param[out] valuelen		データの有効ビット長
* @param[out] image_type	データタイプ
*
* @retval true　 予約された属性である。（除外するパラメータである。）
* @retval false　予約された属性でない。（除外するパラメータでない。）
* 
*/
bool CdbStore::getKeywordValue( const char* key,const char* value,
		int* framenum,int* framesize,int *xsize,int *ysize,int *valuelen,
		char*image_type)
{
	if( 0 == strcmp(key,ARC_PARAM_FRAME) ) {
		*framenum = atoi(value);
		return( true );
	}
	else if(  0 == strcmp(key,ARC_PARAM_FRAME_SIZE) ) {
		*framesize = atoi(value);
		return( true );
	}
	else if(  0 == strcmp(key,ARC_PARAM_FRAME_X) ) {
		*xsize = atoi(value);
		return( true );
	}
	else if(  0 == strcmp(key,ARC_PARAM_FRAME_Y) ) {
		*ysize = atoi(value);
		return( true );
	}
	else if(  0 == strcmp(key,OTH_PARAM_VALUE_LEN) ) {
		*valuelen = atoi(value);
		return( true );
	}
	else if(  0 == strcmp(key,ARC_PARAM_IMAGE_TYPE) ) {
		strcpy(image_type,value);
		return( false );
	}
	return( false );
}
static char* exceptParams[]={
	ARC_PARAM_DIAG_NAME, ARC_PARAM_SHOT, ARC_PARAM_SUB_SHOT,
	ARC_PARAM_DATATYPE, ARC_PARAM_MOD_GROUP, ARC_PARAM_MOD_TYPE,
	ARC_PARAM_MANAGE,
	ARC_PARAM_CH,  ARC_PARAM_IMAGE_TYPE, ARC_PARAM_DATA_LEN,
	ARC_PARAM_COMP_METHOD,
	NULL
};
/**　除外するパラメータを判断する。　（旧パラメータ用）
*
* @param[in]  key			パラメータのキー（項目名）
*
* @retval true　 除外するパラメータである。
* @retval false　除外するパラメータでない。
* 
*/
bool CdbStore::isExceptParam( const char* key )
{
	for(int i=0;;i++) {
		if( NULL == exceptParams[i] ) break;
		if( 0 == _stricmp(key,exceptParams[i]) ) return true;
	}
	return false;	
}
static char* exceptFrameParams[]={
	ARC_PARAM_FRAME, ARC_PARAM_FRAME_SIZE,
	ARC_PARAM_FRAME_X, ARC_PARAM_FRAME_Y,ARC_PARAM_CRC32,
	NULL
};
/**　除外するパラメータを判断する。　（新パラメータ用）
*
* @param[in]  key			パラメータのキー（項目名）
*
* @retval true　 除外するパラメータである。
* @retval false　除外するパラメータでない。
* 
*/
bool CdbStore::isExceptFrameParam( const char* key )
{
	for(int i=0;;i++) {
		if( NULL == exceptFrameParams[i] ) break;
		if( 0 == _stricmp(key,exceptFrameParams[i]) ) return true;
	}
	return false;	
}
/**　データファイルから、データを取得する。
*
* @param[in]  pfname		データファイル名
* @param[out] data_length	データバイト長
* @param[out] data_data		データ(格納領域)
*
* @retval 0 成功
* @retval -1 ファイルオープンエラー
* @retval -2 ファイルリードエラー  
* 
* @note データ格納領域は、使用側で管理すること。
* 
*/
int CdbStore::getData( char* dfname, unsigned long*data_length,unsigned char **data_data)
{
	struct stat buf;
	int ret = stat(dfname,&buf);
	if( 0 != ret  ) return(-1); 

	unsigned long length = buf.st_size;
	
	FILE *fp=fopen( dfname,"rb");
	if( NULL == fp ) return(-1);

	unsigned char *dbuf= new unsigned char[length];
	size_t 	nread	= 0;
	size_t 	rlen	= 0;
	ret	= 0;
	for(;rlen < length;) {
		nread = fread(dbuf+rlen,1,length-rlen,fp);
		if( 0 != ferror(fp) ) {
			ret = -2;
			break;
		}
		rlen += nread;
	}
	fclose(fp);
	if( 0 == ret ) { 
		*data_length= length;
		*data_data 	= dbuf; 
	}
	else {
		delete [] dbuf;
	}
	return ret;	
}
