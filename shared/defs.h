#ifndef DEFS_H
#define DEFS_H 1

#ifdef WIN32
#include <windows.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <tcl/tcl.h>

const int STRLEN = 100;
typedef char STR[STRLEN];

enum GlobalStates
{ 
	STATE_PLAYBACK, 
	STATE_ANIMATE,
	STATE_SIMULATE,
	STATE_IDLE
};

struct GLMouseButtonInfo
{
	int button;
	int state;
};

#endif