/* get_models.c */
int get_models(char modellist[], char modelname[], int ichoice);
int model_locate(char *modelname);
int get_one_model(char filename[]);
int model(char *modelname, double par[]);
int init_models(void);
int assign_wavelengths(struct Model *onemod, int nwaves, float waves[]);
int parse_xy(char *line, int xcol, int ycol, float *w, float *f);
/* get_python_models.c */
int get_python_models(char modellist[], char modelname[], int ichoice);
int get_python_model(char filename[]);
/* ascii_read.c */
int ascii_read(char *filename);
int check_comment(char *xstring, char *answer);
