#include "common.h"

#ifdef USING_R
# include <R.h>
# include <Rinternals.h>
void oops(char *name, char *msg)
{
	error("%s: Error: %s\n",name, msg);
}
void smsg(char *name, char *msg)
{
	Rprintf("%s: %s\n", name, msg);
}
#else
/* OOPS - print a string and exit */
void oops(char *name, char *msg)
{
	fprintf(stderr,"%s: Error: %s\n",name, msg);
	exit(1);
}
void smsg(char *name, char *msg)
{
	fprintf(stderr,"%s: %s\n", name, msg);
}
#endif


