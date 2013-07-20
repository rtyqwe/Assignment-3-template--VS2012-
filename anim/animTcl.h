
#ifndef	_ANIMTCL_H_
#define	_ANIMTCL_H_ 1

//#include <shared/GlobalResourceManager.h>

#include <tcl/tcl.h>
#include <tcl/tk.h>

#include <shared/defs.h>

#include <stdarg.h>

typedef	struct {
		char name[STRLEN];
		int type;
		char *ptr;
} SETVAR;


class animTcl {
	private:
	
	static Tcl_Interp * m_interpreter;
    static int UsingTk;
	static SETVAR setlist[];

	static void appendMessage(Tcl_Interp *interp, char *message);

	public:

	// Utility routines
	static void checkEvents(void);
    static void ExecuteCommand(char *command) ;
	static void OutputMessage(char *format, ...);
	static void OutputResult(char *format, ...);
	static void OutputListElement(char *format, ...);


	// Initialization routines
	static void InitTclTk(int argc,	char **argv);
	static int Tcl_AppInit(Tcl_Interp *interp);
	static void InitUserTclCommands(Tcl_Interp *interp);
	static void LinkUserVariables(Tcl_Interp *interp);

	static int Command(ClientData clientData, Tcl_Interp *interp, int argc,	char **argv);
	static int quit(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);
	static int startGlut(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

};
#endif



