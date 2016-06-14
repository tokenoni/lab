#if defined(WIN32) || defined(WIN64)
#include <string.h>
#include <stdio.h>
struct option {
	const char *name;
	int has_arg;
	int *flag;
	int val;
};

/* �X���b�h�Z�[�t�ł͂���܂���B */
char*	optarg = NULL; 
int		optind = 1; 
int		opterr = 0;
static const struct option* LONG_OPTIONS= NULL;
static int*	LONG_OPTIONS_NAME_LEN		= NULL;

/* �ȉ��̎d�l�������ɃT�|�[�g���Ă��܂���B

�f�t�H���g�ł� getopt() �� argv ���X�L��������ۂɏ�����ύX���A�I�v�V����
�ȊO�̗v�f���Ō�Ɉړ�����B  ���ɂ� 2 �̃��[�h����������Ă���B

optstring �̐擪������ '+' �ł��邩�A���ϐ� POSIXLY_CORRECT ���ݒ肳���
����ꍇ�ɂ́A�I�v�V������ΏۂƂ��铮��́A��I�v�V�����̈����������ꂽ�i
�K�ŏI������B 

optstring �̐擪������ '-' �ł���ꍇ�ɂ́A�I�v�V�����łȂ� argv �v�f�́A��
���R�[�h 1 �̃I�v�V�����ł��邩�̂悤�Ɉ�����  (�����p����v���O�����́A
�I�v�V������ argv �v�f��C�ӂ̏����Ŏ󂯓���A�������̏������Ӗ�������
���ɏ�����Ă���K�v������)�B

"--" �͓���Ȉ������ŁA�X�L�����̃��[�h�ɂ�炸�A�I�v�V�����̃X�L����������
�I�ɏI��������B 

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
/* �����ɃT�|�[�g���Ă��܂���B */
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

