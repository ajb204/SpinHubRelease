
/***/
/* Include file for file inquire:
/***/

#define SORT_NONE  0
#define SORT_ALPHA 1
#define SORT_DATE  2

int dirExists();
int fileExists();
int fileExists0();
int getDirName();
int getFileName();
int isTTY();

char *fileDir();
char *fileExt();
char *fileRoot();
char *fileTail();

