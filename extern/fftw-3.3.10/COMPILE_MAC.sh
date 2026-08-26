## note we need fftw3f_theads and fftw3.
# for this we need to include openMP.
# for mac, this seemed to require gfortran, for nonobvious reasons.
# enable-float is essential for our applications in decon.

#this is the line that homebrew dumps if set to compile from source
#./configure --enable-single --enable-shared --enable-threads --enable-mpi --enable-openmp --enable-sse2 --enable-avx --enable-avx2

#played with a few permutations:
#on my intel mac, best performance by far is with system compiler (not homebrew g++).


#./configure CC=gfortran --enable-threads CFLAGS="-fopenmp" --enable-float
#./configure CC=gcc-13 --enable-threads CFLAGS="-fopenmp" --enable-float
#./configure CC=gcc --enable-threads CFLAGS="-fopenmp" --enable-float
#--enable-i386-hacks
# --enable-openmp

#if( -e ../libs) rm -rf ../libs
#mkdir ../libs

#do floating point libraries
#./configure --enable-single --enable-shared --enable-threads --enable-mpi  --enable-sse2 --enable-avx --enable-avx2
./configure --enable-single --enable-threads 
make -j8
mv .libs/libfftw3f.a ../libs
mv libfftw3f.la ../libs
mv threads/.libs/libfftw3f_threads.a ../libs
mv threads/libfftw3f_threads.la ../libs
##cleanup.
make clean
rm config.log
rm libtool
find . -name ".deps" -exec rm -rf {} \;
find . -name "Makefile" -exec rm  {} \;


#do 'normal' double libraries
./configure --enable-threads
make -j8
mv .libs/libfftw3.a ../libs
mv libfftw3.la ../libs
mv threads/.libs/libfftw3_threads.a ../libs
mv threads/libfftw3_threads.la ../libs
#cleanup.
make clean
rm config.log
rm libtool
find . -name ".deps" -exec rm -rf {} \;
find . -name "Makefile" -exec rm  {} \;




#mv .libs/libfftw3f.a ../libs
#mv libfftw3f.la ../libs
#mv threads/.libs/libfftw3f_threads.a ../libs
#mv threads/libfftw3f_threads.la ../libs
#cleanup.
#make clean
#rm config.log
#rm libtool
#find . -name ".deps" -exec rm -rf {} \;
#find . -name "Makefile" -exec rm  {} \;



#mv .libs/libfftw3f.3.dylib ../libs/
#mv .libs/libfftw3f.dylib ../libs/
#mv .libs/libfftw3f.lai ../libs



#mv threads/.libs/libfftw3f_threads.3.dylib ../libs/
#mv threads/.libs/libfftw3f_threads.dylib ../libs/
#mv threads/.libs/libfftw3f_threads.lai ../libs



#cleanup.
#make clean
#rm config.log
#rm libtool
#find . -name ".deps" -exec rm -rf {} \;
#find . -name "Makefile" -exec rm  {} \;
