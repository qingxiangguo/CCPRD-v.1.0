#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "INLINE.h"
    /* Strip all new line (\n) and carriage return (\r) characters
       from string str
    */
    char* _strip_crnl(char* str) {
        char *s;
        char *s2 = str;
        for (s = str; *s; *s++) {
            if (*s != '\n' && *s != '\r') {
              *s2++ = *s;
            }
        }
        *s2 = '\0';
        return str;
    }

MODULE = Bio::DB::IndexedBase_168b  PACKAGE = Bio::DB::IndexedBase  

PROTOTYPES: DISABLE


char *
_strip_crnl (str)
	char *	str

