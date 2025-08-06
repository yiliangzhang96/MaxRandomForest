#ifndef TYPE_H
#define TYPE_H

#include "stl.h"

#if defined(_MSC_VER) || defined(__BORLANDC__)
	#define WINDOWS
#endif

#ifdef WINDOWS
	typedef __int64 int64;
	typedef unsigned __int64 uint64;
#else
	#include <stdint.h>
	typedef int64_t int64;
	typedef uint64_t uint64;
#endif

typedef unsigned char uchar;
typedef unsigned int uint;
		
#define __DEBUG

#ifdef __ASSERT
	#undef __ASSERT
#endif

#ifdef __DEBUG
	//#define __ASSERT(b, msg) if (!(b)) panic(msg);
	#define __ASSERT(b, msg)  assert(b);
#else
	#define __ASSERT
#endif

inline void panic(const string errormsg){
	cout << errormsg.c_str() << endl;
	cin.get();
	exit(1);
}

#endif //TYPE_H
