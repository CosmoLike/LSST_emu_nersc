
opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -O0 -g -O3 -std=gnu99 -ffast-math -funroll-loops 
opt_ocelote := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.1/include -L/cm/shared/uaapps/gsl/2.1/lib \
-lfftw3 -lgsl -lgslcblas -lm -O0 -g -O3 \
-ffast-math -funroll-loops -std=gnu99 

home:
	make home_lib
	make home_test

home_lib:
	gcc $(opt_home) -shared -o like_fourier.so -fPIC like_fourier.c $(opt_home) 


home_test:
	gcc  like_fourier.c -o ./like_fourier $(opt_home) 


ocelote: 
	make ocelote_lib


ocelote_lib:
	gcc -shared -o like_fourier.so -fPIC like_fourier.c $(opt_ocelote)  
