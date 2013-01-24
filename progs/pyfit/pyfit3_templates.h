/* pyfit3.c */
int main(int argc, char *argv[]);
int help(void);
int update_unique(float xx[], int nnn);
int write_unique(FILE *optr, int n);
int write_record(FILE *optr, int nrec, int nreg);
int write_output_spectrum(char filename[], int nreg);
int get_data(char filename[]);
int get_abs(char filename[]);
int model_multiply(struct spect *model);
int get_python_model(char filename[]);
int regrid(struct spect *modin, struct spect *modout);
float chisq_pure(struct spect *mod);
float metric_abs(struct spect *mod);
int parse_waverange(char waverange[]);
