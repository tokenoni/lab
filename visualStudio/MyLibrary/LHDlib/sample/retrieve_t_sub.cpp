#include "stdio.h"
#include "retrieve.h"

int	ret_t_error_print( int ret , const char* opt24, const char* opt25)
{
	switch(ret)	{
		case 10:
			printf("Parameter is illegal.\n");
			break;
		case 20:
			printf("Memory allocation error!\n");
			break;
		case 24:
			printf("Fail to connect Index Server from %s!\n", opt24);
			break;
		case 25:
			printf("Fail to connect Trans Server from %s!\n", opt25);
			break;
		case 29:
			printf("Index protocole version disagreement!\n");
			break;
		case 30:
			printf("Arc shot data retrieve on DATABASE but not found!\n");
			break;
		case 31:
			printf("Socket time out!\n");
			break;
		case 32:
			printf("Trans protocole version disagreement!\n");
			break;
		case 33:
			printf("Fail to get Data Location from Index Server!\n");
			break;
		case 34:
			printf("Fast channel number within extent!\n");
			break;
		case 35:
			printf("Last channel number within extent!\n");
			break;
		case 36:
			printf("Can not extract of compressed data!\n");
		case 37:
			printf("Arc shot data retrieve on O2 data base but not found!(Abnomal porosess of Transd)\n");
			break;
		case 40:
			printf("Arc shot channel data transed error!\n");
			break;
		case 42:
			printf("Unmatch Attribute\n");
			break;
		case 53:
			printf("Not supported DTS info\n");
			break;
		case 54:
			printf("No DTS link Data. Please check setup for Diagnostics\n");
			break;
		case 55:
			printf("No Diagnostics name\n");
			break;
		case 56:
			printf("No Host ID\n");
			break;
		case 57:
			printf("No DTS Data.\n");
			break;
		default:
			printf("Error Code:%d",ret);
			if(ret < 0) {
				printf("\t[ %s ]", retrieveErrorMessage(ret));
			}
			printf("\n");
			break;
	}

	return ret;
}
