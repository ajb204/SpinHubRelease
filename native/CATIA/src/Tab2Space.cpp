#include <StringMethods.h>
/*
  Clear a character buffer.
*/
void Tab2Space(char* buf, int n) {
  for (int i = 0; i < n; i++) {
    if (buf[i]=='\t'){
      buf[i]=' ';
    };
  };
};
