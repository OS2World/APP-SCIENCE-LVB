 
#define DRAW
#include "phylip.h"

#ifdef X
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#ifdef MAC
#include "interface.h"
#endif

/* Added by Dan F. for bmp code */
#include "math.h"  
#define  DEFAULT_STRIPE_HEIGHT 20

#define maxnodes        1200
#define minus           '-'
#define stripewidth     3000L
#define maxstripedepth  3500
#define fontsize        3800
#define pi              3.1415926535897932384626433
#define epsilond        0.00001
#define ebcdic          EBCDIC
#define segments        40
#define xstart          10
#define ystart          35
#define LF              10
#define CR              13
#define escape  (ebcdic ?  '\'' :  '\033')
#define null  '\000'
#define AFMDIR "/usr/lib/transcript/" /* note trailing slash */

typedef unsigned char byte;
/*typedef char byte; */
typedef enum {treepen, labelpen} pentype;
typedef enum {lw,hp,tek,ibm,mac,houston,decregis,epson,oki,fig,
		citoh,toshiba,pcx,pcl,pict,ray,pov,xpreview,xbm,bmp,
		gif,idraw,vrml,winpreview,other} plottertype;
typedef enum {vertical, horizontal} growth;
typedef enum {cladogram,phenogram,curvogram,
              eurogram,swoopogram,circular} treestyle;
typedef enum {penup,pendown} pensttstype;
typedef enum {plotnow, changeparms, quitnow} winactiontype;
typedef short fonttype[fontsize];
typedef Char *striparray;
typedef striparray striptype[maxstripedepth];

struct LOC_plottext {              /* Local variables for plottext: */
  double height, compress;
  short *font;
  short coord;
  double heightfont, xfactor, yfactor, xfont, yfont, xplot, yplot, sinslope,
         cosslope, xx, yy;
  pensttstype penstatus;
} ;

typedef struct colortype {
  Char *name;
  double red, green, blue;
} colortype;

double lengthtext();
double heighttext();

/* For povray, added by Dan F. */
#define TREE_TEXTURE "T_Tree\0"
#define NAME_TEXTURE "T_Name\0"

#ifdef X
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#ifdef X
Display *display;               /* the X display */
Window  mainwin;                /* the main display window */
int x, y;                       /* the corner of the window */
int width, height;              /* the width and height of the window */
#define FONT "-*-new century schoolbook-medium-r-*-*-14-*"
XFontStruct *fontst;            /* the font strcture for the font */
char *fontrsc;                  /* the font resource */
XFontStruct *fontst;             /* the font strcture for the font */
XGCValues gcv;                   /* graphics context values */
GC gc1;                          /* a graphics context */
#define DEFGEOMETRY "600x400+20+50"
#endif
#define LARGE_BUF_LENGTH 500
extern char fontname[LARGE_BUF_LENGTH]; /* the font name to use */

#ifdef WIN32
#define DEFPLOTTER lw
#define DEFPREV winpreview
#endif
#ifdef DOS
#define DEFPLOTTER lw
#define DEFPREV ibm
#endif   
#ifdef MAC 
#define DEFPLOTTER pict
#define DEFPREV mac 
#endif   
#ifdef VMS
#define DEFPLOTTER lw
#define DEFPREV decregis
#endif
#ifndef DOS
#ifndef MAC 
#ifndef VMS
#ifndef WIN32
#define DEFPLOTTER lw
#ifdef X
#define DEFPREV xpreview
#endif
#ifndef X
#define DEFPREV tek
#endif
#endif
#endif   
#endif   
#endif   

/* Define SEEK_SET (needed for fseek()) for machines that haven't 
   got it already, */
#ifndef SEEK_SET
#define SEEK_SET 0
#endif


/* prototypes should be here */

void   plotdot(long, long);
void   circlepoints(int, int, int, int);
void   drawpen(long, long, long);
void   drawfatline(long, long, long, long, long);
void   plot(pensttstype, double, double);
void   idellipse(double, double); 
void   splyne(double,double,double,double,boolean,long,boolean,boolean);
void   swoopspline(double,double,double,double,double,double,boolean,long);
void   curvespline(double, double, double, double, boolean, long);
static void   putshort(FILE *, int);

static void   putint(FILE *, int);
void   write_bmp_header(FILE *, int, int);
void   reverse_bits (byte *, int); 
void   turn_rows(byte *, int, int);
void   translate_stripe_to_bmp(striptype *, byte *, int, int, int, int *); 
void   write_full_pic(byte *, int);
void   makebox_no_interaction(char *, double *, double *, double *, long);
boolean plot_without_preview(char *, double *, double *, double *, long,
		node *);
void   void_func(void);
double computeAngle(double, double, double, double);
void plottree(node*,node*);
void plotlabels(char*);
