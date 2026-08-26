#include <StringMethods.h>
/*
  Clear a character buffer.
*/
void ClearBuf(char* buf, int n) {
  for (int i = 0; i < n; i++) {
    buf[i] = 0;
  };
};
