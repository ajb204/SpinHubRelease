
#if defined (WINNT) || defined (WIN95)
#ifdef min 
#undef min
#endif

#ifdef max
#undef max
#endif
#endif

#define min( I1, I2 ) (((I1) < (I2)) ? (I1) : (I2))
#define max( I1, I2 ) (((I1) > (I2)) ? (I1) : (I2))
#define sq( I1 ) ((I1)*(I1))
