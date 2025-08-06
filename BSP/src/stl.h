#ifndef STL_H
#define STL_H

//include all the necessary stl header files
#include <vector>
#include <stack>
#include <queue>
//#include <list>
//#include <deque>
#include <map>
//#include <set>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <assert.h>
#include <sstream>
//use namespace std to make the the source codes clear
using namespace std;

//+ by CHTsai
#define isnan(x) ::_isnan(x)
#define isinf(x) (!::_finite(x))

inline bool compareinterval (const pair<double, double> &i, const pair<double, double> &j) {
	return (i.first<j.first);
}

#endif //STL_H
