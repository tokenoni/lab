#if defined(WIN32) || defined(WIN64)
#include <string.h>
#include <stdio.h>
struct option {
	const char *name;
	int has_arg;
	int *flag;
	int val;
};

/* スレッドセーフではありません。 */
char*	optarg = NULL; 
int		optind = 1; 
int		opterr = 0;
static const struct option* LONG_OPTIONS= NULL;
static int*	LONG_OPTIONS_NAME_LEN		= NULL;

/* 以下の仕様を厳密にサポートしていません。

デフォルトでは getopt() は argv をスキャンする際に順序を変更し、オプション
以外の要素を最後に移動する。  他にも 2 つのモードが実装されている。

optstring の先頭文字が '+' であるか、環境変数 POSIXLY_CORRECT が設定されて
いる場合には、オプションを対象とする動作は、非オプションの引き数が現れた段
階で終了する。 

optstring の先頭文字が '-' である場合には、オプションでない argv 要素は、文
字コード 1 のオプションであるかのように扱われる  (これを用いるプログラムは、
オプションや argv 要素を任意の順序で受け入れ、かつそれらの順序が意味を持つよ
うに書かれている必要がある)。

"--" は特殊な引き数で、スキャンのモードによらず、オプションのスキャンを強制
的に終了させる。 

*/

int  getopt (int argc, char **argv, char *optstring) 
{
	bool opt_1_ret = false;
	if( '-' == optstring[0] ) opt_1_ret = true;

	for(int i=optind;i<argc;i++) {
		if( '-' == argv[i][0] ) {
			optarg = NULL;
			char *opt=strchr(optstring,argv[i][1]);
			if( NULL == opt ) {
				opt = "?";
			}
			else {
				if( ':' == *(opt+1) ){
					if( 2 < strlen(argv[i]) ) {
						optarg = &argv[i][2];
					}
					else {
						if( argc>(i+1) ) {
							if( '-' != argv[i+1][0] ) {
								optarg = argv[++i];
							}
						}
					}
				}
			}
			optind = ++i;
			return *opt;
		}
		else {
			if( opt_1_ret ) {
				optind = ++i;
				return 1;
			}
		}
	}
	return(-1);
}
/* 厳密にサポートしていません。 */
int getopt_long(int argc, char * const argv[],
                const char *optstring,
                const struct option *longopts, int *longindex)
{
	bool opt_1_ret = false;
	if( '-' == optstring[0] ) opt_1_ret = true;

	if( LONG_OPTIONS != longopts ) {
		LONG_OPTIONS = longopts;
		if( NULL != LONG_OPTIONS_NAME_LEN ) delete [] LONG_OPTIONS_NAME_LEN;
		int l_cnt = 0;
		for(;;l_cnt++ ) {
			if( NULL == longopts[l_cnt].name ) break;
		}
		LONG_OPTIONS_NAME_LEN = new int[l_cnt];
		for(int n=0;n<l_cnt;n++ ) {
			LONG_OPTIONS_NAME_LEN[n] = (int)strlen(longopts[n].name);
		}
	}
	for(int i=optind;i<argc;i++) {
		if( '-' == argv[i][0] ) {
			optarg = NULL;
			if( '-' == argv[i][1] ) {
				if( NULL == longopts ) continue;
				for(int ll=0;;ll++ ) {
					if( NULL == longopts[ll].name ) break;
					int n_len = LONG_OPTIONS_NAME_LEN[ll];
					if( 0 == strncmp(&argv[i][2], longopts[ll].name, n_len) ) {
						if( '\0' == argv[i][2+n_len] ) {
							if( longopts[ll].has_arg ) {
								optarg = argv[++i];
							}
						}
						else if( '=' == argv[i][2+n_len] ) {
							if( longopts[ll].has_arg ) {
								optarg = &argv[i][2+n_len];
							}
						}
						optind = ++i;
						if( NULL != longindex ) *longindex = ll;
						if( NULL == longopts[ll].flag )	return longopts[ll].val;
						else	{
							*longopts[ll].flag = longopts[ll].val;
							return 0;
						}
					}
				}
				continue;
			}
			const char *opt=strchr(optstring,argv[i][1]);
			if( NULL == opt ) {
				opt = "?";
			}
			else {
				if( ':' == *(opt+1) ){
					if( 2 < strlen(argv[i]) ) {
						optarg = &argv[i][2];
					}
					else {
						if( argc>(i+1) ) {
							if( '-' != argv[i+1][0] || '\0' == argv[i+1][1] ) {
								optarg = argv[++i];
							}
						}
					}
				}
			}
			optind = ++i;
			return *opt;
		}
		else {
			if( opt_1_ret ) {
				optarg =  argv[i];
				optind = ++i;
				return 1;
			}
		}
	}
	return(-1);
}
#endif

