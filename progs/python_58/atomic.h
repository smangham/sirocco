
#define MAXRAND 2147486748.
#define INFINITY 1e50
#define TRUE		1
#define FALSE		0

#define H   				6.6262e-27
#define HC  				1.98587e-16
#define HEV				4.13620e-15	/* Planck's constant in eV */
#define C   				2.997925e10
#define G				6.670e-8
#define BOLTZMANN 			1.38062e-16
#define H_OVER_K			4.799437e-11
#define STEFAN_BOLTZMANN 		5.6696e-5
#define THOMPSON			0.66524e-24
#define PI  				3.1415927
#define MELEC 				9.10956e-28
#define E				4.8035e-10	/* Electric charge in esu */
#define MPROT 				1.672661e-24
#define MSOL 				1.989e33
#define PC				3.08e18
#define YR   				3.1556925e7
#define PI_E2_OVER_MC			0.02655103	/* Classical cross-section */
#define PI_E2_OVER_M  			7.96e8
#define ALPHA				7.297351e-3	/* Fine structure constant */
#define BOHR				0.529175e-8	/* Bohr radius */
#define	CR				3.288051e15	/*Rydberg frequency for H != Ryd freq for infinite mass */
#define ANGSTROM                        1.e-8           /*Definition of an Angstrom in units of this code, e.g. cm */

#define EV2ERGS   			1.602192e-12
#define RADIAN				57.29578
#define RYD2ERGS                        2.1798741e-11	/* Rydberg in units of ergs */


/* The next term attempts gloabally to define a minimum density to prevent zero devides in some routines */
#define DENSITY_MIN		1.e-20
/* 

   Note that the structure ele array may not be completely filled.  In this structure, the dimension is simply
   the order in which elements are read in, and one may skip elements (which are not of interest).

   So that one can easily find, say carbon, however, *elz contains pointers to ele such that elz[6] is the
   pointer to carbon.  Elements for which there is no data, point to the dummy element structure foo_element

   If you want to cycle through the elements one will act on ele directly.  If you want to know something
   about a specific element *elz may be better, at least that is the idea

 */

#define NELEMENTS		50	/* Maximum number of elements to consider */
int nelements;			/* The actual number of ions read from the data file */
#define NIONS		500	/* Maximum number of ions to consider */
int nions;			/*The actual number of ions read from the datafile */
#define NLEVELS 	5000	/* Maximum number of levels for all elements and ions */
int nlevels;			/*These are the actual number of levels which were read in */
#define NLTE_LEVELS	200     /* Maximum number of levels to treat explicitly */
int nlte_levels;		/* Actual number of levels to treat explicityly */
#define NLEVELS_MACRO   50     /* Maximum number of macro atom levels. (SS, June 04) */
int nlevels_macro;              /* Actual number of macro atom levels. (SS, June 04) */
#define NLINES 10000            /* Maximum number of lines to be read */
int nlines;                     /* Actual number of lines that were read in */
int nlines_macro;                /* Actual number of Macro Atom lines that were read in.  New version of get_atomic
				data assumes that macro lines are read in before non-macro lines */

#define NBBJUMPS         20      /* Maximum number of Macro Atom bound-bound jumps from any one configuration (SS) */

#define NBFJUMPS         20      /* Maximum number of Macro Atom Bound-free jumps from any one configuration (SS) */

#define MAXJUMPS          1000000      /* The maximum number of Macro Atom jumps before emission (if this is exceeded
                                     it gives up (SS) */
typedef struct elements
  {				/* Element contains physical parameters that apply to element as a whole and 
				   provides and index to the atomic data */
    char name[20];		/* Element name */
    int z;			/* Atomic number */
    int firstion;		/*Index into struct ions  ion[firstion] is the lowest ionization state of this ion */
    int nions;			/*The number of ions actually read in from the data file for this element */
    double abun;		/*Abundance */
  }
ele_dummy,*ElemPtr; 

ElemPtr ele;

double rho2nh;			/* The conversion constant from rho to nh the total number of H atoms */

/* Note that ion is the basic structure.  It is filled from 0 up to nions.  *ionzi is a set of pointers which 
   can be accessed via z and i (the ionization state).  But it is important to know that the pointer is really
   not pointed at foo_ion */


typedef struct ions
  {
    int z;			/* defines the element for this ion */
    int istate;			/*1=neutral, 2 = once ionized, etc. */
    double ip;			/* ionization potential of this ion (converted from eV to ergs by get_atomic) */
    double g;			/* multiplicity of ground state, note that this is not totally consistent
				   with energy levels and this needs to be reconciled */
    int nmax;			/* The maximum number of allowed lte configurations for this ion. 
				Used only be get_atomic */
    int n_lte_max;		/* The maximum number of allowed nlte configurations for this ion. 
				Used only be get_atomic */
    int firstlevel;		/* Index into config struc; should also refer to ground state for this ion */
    int nlevels;		/* Actual number of "lte" levels for this ion. "lte" here refers to any level 
				in which densities are to be calculated on the fly, instead of the beginning 
				or end of an ionization cycle */
    int first_nlte_level;       /* Index into config structure for the nlte_levels. There will be nlte 
				levels associated with this ion*/
    int nlte_first;             /* Index into the levden array in wind structure */
    int nlte;                   /* Actual number of nlte levels for this ion */
    int phot_info;		/*-1 means there are no photoionization cross sections for this ion, 
				 0  means the photoionization is given on an ion basis (e.g. for the 
					ground state using VFKY and the xphot structure
				   1 means photoionization is given on a level basis, using topbase value
                                   Topbase photoinization trumps VFKY photoionization
				  */
    int macro_info;             /* Identifies whether ion is to be treated using a Macro Atom approach.
				   set to -1 initially (means not known) 
				   set to 1 if macro atom method is used
				   set to 0 if macro atom method is not used  (SS) for this ion (ksl)
				Note: program will exit if -1 before leaving get_atomicdata
				 */
    int ntop_first;             /* Internal index into topbase photionization structure */
    int ntop;                   /* Number of topbase photoionization cross sections for this ion */
    int nxphot;			/* Internal index into VFKY photionionization structure.  There
				is only one of these per ion */
  }
ion_dummy,*IonPtr;

IonPtr ion;


/* And now for the arrays which describe the energy levels.  In the Topbase data, g is float (although
I have never seen it as anything else.  I have kept g an integer here, because I am pretty sure it is
an integer everywhere else in the code....but this is something to watch ??? */

typedef struct configurations
  {
    int z;                      /*The element associated with this configuration*/
    int istate;                 /*The ion associated with the configuration */
    int nion;			/*Internal index to ion structure */
    int nden;                   /*Internal index to levden array in wind structure. -1 if no such level */
    int isp;			/*Topbase description of angular momentum (2s+1)*100+l*10+p */
    int ilv;			/* Unique Level number (within an ion. Topbase ilv when Topbase record*/
    int macro_info;             /* Unambigous way to determine this is part of a macro-level ksl 
				   set to -1 initially (means not known)
				   set to 1 if macro atom method is used
				   set to 0 if macro atom method is not used for this configuration 
				Note: program will exit if -1 before leaving get_atomicdata
				 */
    int n_bbu_jump, n_bbd_jump; /* No. of bound-bound upward/downward jumps available to this configuration (SS) */
    int n_bfu_jump, n_bfd_jump; /* No. of bound-free upward/downward jumps available to this configuration (SS) */
    double g;			/* multiplicity for level */
    double q_num;		/* principal quantum number.  In Topbase this has non-integer values */
    double ex;			/*excitation energy of level */
    double rad_rate;		/* Total spontaneous radiative de-excitation rate for level */
    int bbu_jump[NBBJUMPS];         /* list of upwards bb jumps from this configuration (SS) */
    int bbd_jump[NBBJUMPS];         /* list of downwards bb jumps from this configuration (SS) */
    int bfu_jump[NBFJUMPS];         /* list of upwards bf jumps from this configuration (SS) */
    int bfd_jump[NBFJUMPS];         /* list of downwards bf jumps from this configuration (SS) */
  }
config_dummy,*ConfigPtr;

ConfigPtr config;

/* So what is the energy of the first level CIV 
   ex[ion[6][0].index]
 */

/* Structure to describe the lines */

typedef struct lines
  {
    int nion;			/*The ion no (in python) of the transition */
    int z, istate;		/*element and ion associated with the line */
    double gl, gu;		/*multiplicity of lower and upper level respectively */
    int nconfigl,nconfigu;      /*The configuration no (in python) of the transition */
    int levl,levu;		/*level no of transition..parallel/redundant with el,eu hopefully */
    int macro_info;             /* Identifies whether line is to be treated using a Macro Atom approach.
				   set to -1 (not known initially) 
				   set to 0 if not a macro atom line  
				   set to 1 if a macro atom line  (ksl 04 apr)
				Note: program will exit if -1 before leaving get_atomicdata
				 */
    double freq;		/* The frequency of the resonance line */
    double f;			/*oscillator strength.  Note: it might be better to keep PI_E2_OVER_MEC flambda times this.
				   Could do that by initializing */
    double el,eu;               /* The energy of the lower and upper levels for the transition */
    double pow;			/*The power in the lines as last calculated in total_line_emission */
    int where_in_list;          /* Position of line in the line list: i.e. lin_ptr[line[n].where_in_list] points
				   to the line. Added by SS for use in macro atom method. */
    int down_index;             /* This is to map from the line to knowing which macro atom jump it is (and therefore find
				   the estimator with which it is associated. The estimator is identified by the
				   upper configuration (nconfigu) and then down_index (for deexcitation) or the lower
				   configuration (nconfigl) and then up_index. (SS) */
    int up_index;
  }
line_dummy,*LinePtr;


LinePtr line,lin_ptr[NLINES];	/* line[] is the actual structure array that contains all the data, *lin_ptr
				   is an array which contains a frequency ordered set of ptrs to line */
                                /* fast_line (added by SS August 05) is going to be a hypothetical
                                   rapid transition used in the macro atoms to stabilise level populations */
struct lines fast_line;

int nline_min, nline_max, nline_delt;	/*For calculating which line are likely to be in resonance, see
					   the routine limit_lines */



/*structure containing photoionization data */

/* Photoionization crossections from Verner, Ferland, Korista & Yakovlev */
typedef struct photoionization
  {
    int nion;			/* index to the appropriate ion in the structure ions, so for example, ion would
				   normally be 0 for neutral H, 1 for H+, 1 for He, 2 for He+ etc */
    int z, istate;
    double freq_t;		/*frequency of threshold derived from ionization potential */
    double freq_max;		/* maximum frequency for which fit formula apllies */
    double freq0;		/* fit parameter */
    double sigma;		/*cross section at freq0 */
    double ya, p, yw, y0, y1;	/* Fit prarameters */
    double f_last,sigma_last;            /*last freq, last x-section */

  } Photoionization,*PhotoionizationPtr;

Photoionization xphot[NIONS];
PhotoionizationPtr xphot_ptr[NIONS]; /* Pointers to xphot in threshold frequency order --57h -- ksl*/

int nxphot;			/*The actual number of ions for which there are VFKY photoionization x-sections */
double phot_freq_min;		/*The lowest frequency for which photoionization can occur */

#define NCROSS 150
#define NTOP_PHOT 250           /* Maximum number of photoionisation processes. (SS) */
int ntop_phot;			/* The actual number of TopBase photoionzation x-sections */

typedef struct topbase_phot {           /* If the old topbase treatment is to be replaced by Macro Atoms perhaps this
				   can be dropped - need for backward compatibility? (SS)*/
  int nlev;			/* Internal index to the config structure for the lower state */
  int uplev;                    /* Internal index to the config structure for the upper state (SS) */
  int nion;			/* Internal index to the    ion structure for this x-section */
  int z, istate;
  int np;
  int nlast;			/* nlast is an index into the arrays freq and x.  It allows 
				   a quick check to see whether one needs to do a search
				   for the correct array elements, or whether this process
				   can be short circuited */
  int macro_info;             /* Identifies whether line is to be treated using a Macro Atom approach.
				   set to -1 initially
				   set to 0 if not a macro atom line  
				   set to 1 if a macro atom line  (ksl 04 apr)
				   Note: Program will exit before leaving get_atomicdata if not initallized to 0 or 1
				 */
  int down_index;             /* This is to map from the line to knowing which macro atom jump it is (and therefore find
				 the estimator with which it is associated. The estimator is identified by the
				 upper configuration (uplev) and then down_index (for deexcitation) or the lower
				 configuration (nlev) and then up_index. (SS) */
  int up_index;
  double freq[NCROSS],x[NCROSS];
  double f,sigma;            /*last freq, last x-section */
} Topbase_phot,*TopPhotPtr;

Topbase_phot phot_top[NLEVELS];
TopPhotPtr phot_top_ptr[NLEVELS];  /* Pointers to phot_top in threshold frequency order */

/* now do the table of fractions of recombinations going straight to the ground state */
struct ground_fracs
  {
    int z, istate;		/*element and ion */
    double frac[20];		/* this is a table containing the fraction of recombinations going directly
				   to the ground state as a function of temperature. ground_frac[0] is or t=5000
				   and then we go in steps of 5000 to ground_frac[19] which is for t=1e5. these
				   fractions must have been computed elsewhere */
  }
ground_frac[NIONS];


/* Data needed for collisionally excited lines in a thin plasma.   The
   data is taken in atomic.dat from Gaetz and Salpeter.  Note that .tm
   is converted to degrees (from log degrees) in atomic.c */

#define NTRANS 200
int nxcol;			/*number of transition for collisional excitation */
struct collision_strength
  {
    int nion;			/*The index into the ion array */
    int z, istate;
    double freq;		/* frequency of the line */
    double ex;			/* excitation energy of line (converted from eV 2 ergs in get_atomic) */
    double alpha;
    double beta;
    double tm;			/* Temp max, expanded from log in get_atomic */
    double pow;			/*power associated with this line */
  }
xcol[NTRANS], *xcol_ptr[NTRANS];	/* xcol[] is the actual structure array that contains all the data, *xcol_ptr
					   is an array which contains a frequency ordered set of ptrs to line */

int nxcol_min, nxcol_max, nxcol_delt;	/*For calculating a frequency range within which collisions should
					   be included, see the routine limit_collisions */

//#include "recipes.h"
