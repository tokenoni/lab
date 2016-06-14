#ifndef _GETOPT_WIN32_H_
#define _GETOPT_WIN32_H_
#if defined(WIN32) || defined(WIN64)
 extern char*	optarg; 
 extern int		optind; 
 extern int		opterr; 
 int getopt (int argc, char **argv, char *optstring);

 struct option {
	const char *name;
	int has_arg;
	int *flag;
	int val;
 };
 int getopt_long(int argc, char * const argv[],
                const char *optstring,
                const struct option *longopts, int *longindex);
#endif
#endif