

/* Note that there are some limits on how large NMODS*NWAVES can be.  When you exceed this limit, you will get a segmentation
 * fault on startup */

/* Note that there is some overlap between these definitions and those in python.h and so mixing these two
 * header files can be dangerous.  
 */

#define NWAVES 	8000
#define NDIM	10	// The maximum number of free parameters
#define NCOMPS	10     //The maximum number of separate components
#define NPARS	10     //The maximum number of parameters in a component (not all variable)
#define NMODS   80000    //The total number of models read in in all components
#define NXMODS	10      //The number of models used in a complex models.   
#define LINELEN 132



int ncomps; // The number of components that have been read.  
int nmods_tot;  // The total number of models that have been read in 

    
/* This is just a structure to make one can determine what sets of wavelengths
 * are stored for the models
 */

struct Wavelengths 
{
	float *w;
	int nwaves;
} xwav[NMODS];

int n_xwav;   /* The number of different wavelength arrays that exit */
int n_xflux;  /* Tells you how many of the models have populated fluxes */


/* This is the structure that describes an individual continuum model. 
 * mods is the set of all models.  Note that reading a model consists
 * of one stetp to read the file that lists the model and the parameters.  
 * and a second step to actually read the model.  The model may not actually
 * be read in until it is needed; space is allocated at that time
 *
 * The name of a model is not necessarily unique, e. g. in the case
 * of python models.  However the combination of name and nspec should
 * be unique.  
 *
 * Normally, however One should access models through ModSum.  
 */
typedef struct Model
  {
    char name[LINELEN];
    int ncomp;  // Component number for this model
    int nspec;		//Usually 0, unless multiple spectra are read from the 
    			// same ascii file
    double par[NPARS];
    int ix,iy;   // Column number where wavelength(freq) and flux is found
    		 // Implemented so that one can work directly with python
		 // output spectra
    float *w;    // Floats to save space
    float *f;
    int nwaves;
  } *ModelPtr;



struct Model mods[NMODS];

/* There is one element of comp for each set of models of the same type, i.e. if
one reads in a list of WD atmosphers this will occupy one componenet here 

There is a built-in assumption that the header information for all of the
models will be read in at once.  Otherwise modstart and modstop would not
be a useful way to identify the models.
*/

# define SIMPLE  0
# define PYTHON  1  

/* Used by kslfit */
#define T_MULT  0
#define T_ADD   1


typedef struct Comp  
  {
    char var[LINELEN];  //Abbreviation used locate this component for arithmetic
    char name[LINELEN]; //Name of master file describint models; this for historical reasons
    int filetype;       //Type of model (Uniquely defines format for master file and spectral files)
    			// The two types are SIMPLE and PYTHON at present
    int nmods, modstart, modstop;    //number of models of this type and the index into the model struct
    int npars;          // Number of parameters that describe the model
    double min[NPARS];           // The minimum and maximum for each "free" paratmenter of the model
    double max[NPARS];
    int nwaves;         //All models in each comp should have same wavelengths;
    struct Model xmod;  //The current interpolated model of this type 
//Python    struct Pdf	xpdf;	//The current cumulative distribution function for this component

    /* Items below here are specifically for kslfit */
    int type;		// Indicate additive or multiplicative model 
    float new[NPARS];       // Current set of parameters being used -- not necessarily the best fit set
    float orig[NPARS];  // An original set (guess) at the parameters...used by init_fits.
//OLD    float f[NWAVES];    //The flux associated with this component in the wavelengnth bins of the data
  }
*CompPtr;

struct Comp comp[NCOMPS];

// Possible alternative approach -- struct Model comp_model[NCOMPS];  // A parallel model structure containing the last model created by 
//				// routine called model (whose name we need to change as it is too
//				// confusing to something like ModelPtr make_model(CompPtr comp)::

int ncomps;  //Actual number of components to the current model.  N.B. This is not the number of variables

/* In modsum, current[0] often refers to a normalization.  Therefore for parallism, a set of
models which really only have one parameter, e.g. T, will store that parmeter in .par[1] */

/* The Spectrum structure has been left as it is in kslfit, except that the wavelength
 * boundaries have been omitted
 */


typedef struct Spectrum 
	  {                             /* need both w and f for irregularly spaced grids */
	char name[LINELEN];
	float w[NWAVES];
	float f[NWAVES];
	float e[NWAVES];
	float fit[NWAVES];
	int good[NWAVES];           // 0 = nogood, 1 = good   
	int nboot[NWAVES];          // Bootstrap weight  

	int nw, ng;
	int regrid_spectra;       // 0 is a request to leave model spectra alone, 1 to
				// regrid to the wavelength scale here
	int npars;
	float par[NPARS];       // A place to store pars (e.g pyfit4) when dealing with multiple spectral fits

	/* Used by kslfit for something */
//OLD	float f_orig[NWAVES];
}
*SpecPtr;

struct Spectrum spec;


#define MAX_LINES 10000

char ascii_store[MAX_LINES][LINELEN];
int ascii_nlines;


// #include "model_prototypes.h"
