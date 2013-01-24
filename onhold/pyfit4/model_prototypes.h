/* get_data.c */
int get_data(char filename[]);
int regrid(struct Model *mod_in, struct Model *mod_out);
/* get_models.c */
int get_models(char modellist[], char modelname[], int ichoice);
int model_locate(char *modelname);
int get_one_model(char filename[]);
ModelPtr make_model(CompPtr xcomp);
int init_models(void);
int assign_wavelengths(struct Model *onemod, int nwaves, float waves[]);
int parse_xy(char *line, int xcol, int ycol, float *w, float *f);
/* get_python_models.c */
int get_python_models(char modellist[], char modelname[], int ichoice);
int get_python_model(char modelname[]);
int get_python_filename(char *modelname, char *filename);
/* ascii_read.c */
int ascii_read(char *filename);
int check_comment(char *xstring, char *answer);
int get_comment_line(char *xstring, char *answer);
/* model_arith.c */
CompPtr get_component(char *model_name);
