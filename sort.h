#ifndef JLIB_SORT_H
#define JLIB_SORT_H

#include <stdlib.h>

void sort(void *arr, size_t bytes, size_t n, int (*cmp)(void*, void*, long, long, void*), void *arg);

#define JLIB_SORT_IMPL
#ifdef JLIB_SORT_IMPL

void sort(void *arr, size_t bytes, size_t n, int (*cmp)(void*, void*, long, long, void*), void *arg) {
	if(!arr || n == 0) return;

    size_t BUFFER[n&63]; 
	size_t *stack = 0;
	size_t l = 0;
	size_t h = n - 1;

	if(n >= 64)
		stack = malloc(sizeof(size_t) * n);
	else
		stack = BUFFER;
  
    long top = -1; 
  
    stack[++top] = l; 
    stack[++top] = h; 
  
    while(top >= 0) { 
        h = stack[top--]; 
        l = stack[top--]; 
  
        long p = 0; 

		{
			unsigned char TMP[bytes];
			void *x = (unsigned char*)arr + h*bytes; 
			long i = l - 1; 

			for(long j = l; j < h; ++j) { 
				if(cmp((unsigned char*)arr + j*bytes, x, j, h, arg) <= 0) { 
					i++; 
					unsigned char *ip = (unsigned char*)arr + i*bytes;
					unsigned char *jp = (unsigned char*)arr + j*bytes;
					int index;
					for(index = 0; index < bytes; ++index)
						TMP[index] = ip[index];
					for(index = 0; index < bytes; ++index)
						ip[index] = jp[index];
					for(index = 0; index < bytes; ++index)
						jp[index] = TMP[index];
				} 
			} 

			unsigned char *ip = (unsigned char*)arr + i*bytes + bytes;
			unsigned char *hp = (unsigned char*)arr + h*bytes;
			int index;
			for(index = 0; index < bytes; ++index)
				TMP[index] = ip[index];
			for(index = 0; index < bytes; ++index)
				ip[index] = hp[index];
			for(index = 0; index < bytes; ++index)
				hp[index] = TMP[index];

			p = i + 1; 
		}
  
        if(p - 1 > l) { 
            stack[++top] = l; 
            stack[++top] = p - 1; 
        } 
  
        if(p + 1 < h) { 
            stack[++top] = p + 1; 
            stack[++top] = h; 
        } 
    } 

	if(stack != BUFFER)
		free(stack);
}

#endif
#endif
