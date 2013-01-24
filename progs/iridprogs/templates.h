/* bb.c */
double planck(double t, double freqmin, double freqmax);
double get_rand_pow(double x1, double x2, double alpha);
double get_rand_exp(double alpha_min, double alpha_max);
double integ_planck_d(double alphamin, double alphamax);
int init_integ_planck_d(void);
double planck_d(double alpha);
double emittance_bb(double freqmin, double freqmax, double t);
/* get_atomicdata.c */
int get_atomic_data(char masterfile[]);
int index_lines(void);
int index_phot_top(void);
int index_phot_verner(void);
int index_collisions(void);
void indexx(int n, float arrin[], int indx[]);
int limit_lines(double freqmin, double freqmax);
/* python.c */
int main(int argc, char *argv[]);
int init_geo(void);
int photon_checks(PhotPtr p, double freqmin, double freqmax);
/* photon2d.c */
int translate(WindPtr w, PhotPtr pp, double tau_scat, double *tau, int *nres);
int translate_in_space(PhotPtr pp);
double ds_to_wind(PhotPtr pp);
int translate_in_wind(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres);
int walls(PhotPtr p, PhotPtr pold);
/* photon_gen.c */
int define_phot(PhotPtr p, WindPtr w, double f1, double f2, int nphot_tot, int ioniz_or_final, int iwind, int freq_sampling);
double populate_bands(WindPtr w, double f1, double f2, int ioniz_or_final, int iwind, struct xbands *band);
int xdefine_phot(WindPtr w, double f1, double f2, int ioniz_or_final, int iwind);
int xmake_phot(PhotPtr p, WindPtr w, double f1, double f2, int ioniz_or_final, int iwind, double weight, int iphot_start, int nphotons);
double star_init(double r, double tstar, double freqmin, double freqmax, int ioniz_or_final, double *f);
int photo_gen_star(PhotPtr p, double r, double t, double weight, double f1, double f2, int spectype, int istart, int nphot);
double disk_init(double rmin, double rmax, double m, double mdot, double freqmin, double freqmax, int ioniz_or_final, double *ftot);
int photo_gen_disk(PhotPtr p, double weight, double f1, double f2, int spectype, int istart, int nphot);
int phot_gen_sum(char filename[], char mode[]);
double bl_init(double lum_bl, double t_bl, double freqmin, double freqmax, int ioniz_or_final, double *f);
/* saha.c */
int nebular_concentrations(PlasmaPtr xplasma, int mode);
int concentrations(PlasmaPtr xplasma, int mode);
int saha(double ne, double nh, double t, int nelem, double density[]);
int lucy(PlasmaPtr xplasma);
int lucy_mazzali1(double nh, double t_r, double t_e, double www, int nelem, double ne, double density[], double xne, double newden[]);
int fix_concentrations(double nh, double density[], double *ne);
double get_ne(double density[]);
int do_partitions(PlasmaPtr xplasma, int mode);
double partition(int nion, double ne, double t, double w, int pmode);
int dlucy(PlasmaPtr xplasma);
/* spectra.c */
int spectrum_init(double f1, double f2, int nangle, double angle[], double phase[], int scat_select[], int top_bot_select[], int select_extract, double rho_select[], double z_select[], double az_select[], double r_select[]);
int spectrum_create(PhotPtr p, double f1, double f2, int nangle, int select_extract);
int spectrum_summary(char filename[], char mode[], int nspecmin, int nspecmax, int select_spectype, double renorm);
/* wind2d.c */
int define_wind(void);
int where_in_grid(double x[]);
int vwind_xyz(PhotPtr p, double v[]);
int wind_div_v(WindPtr w);
double rho(WindPtr w, double x[]);
int mdot_wind(WindPtr w, double z, double rmax);
int get_random_location(int n, double x[]);
/* wind.c */
int where_in_wind(double x[]);
int wind_check(WindPtr www, int n);
double model_velocity(double x[], double v[]);
int model_vgrad(double x[], double v_grad[][3]);
double model_rho(double x[]);
/* vector.c */
double dot(double a[], double b[]);
double length(double a[]);
int renorm(double a[], double scalar);
int cross(double a[], double b[], double c[]);
int vmove(double u[], double lmn[], double s, double result[]);
int vsub(double u[], double v[], double result[]);
int vadd(double u[], double v[], double result[]);
int stuff_v(double vin[], double vout[]);
double dot_tensor_vec(double tensor[3][3], double vin[3], double vout[3]);
int project_from_xyz_cyl(double a[], double b[], double result[]);
int project_from_cyl_xyz(double a[], double b[], double result[]);
int create_basis(double u[], double v[], struct basis *basis_new);
int project_from(struct basis *basis_from, double v_in[], double v_out[]);
int project_to(struct basis *basis_from, double v_in[], double v_out[]);
int reorient(struct basis *basis_from, struct basis *basis_to, double v_from[], double v_to[]);
/* debug.c */
int DebugStr(char *string);
/* recipes.c */
double qromb(double (*func)(double), double a, double b, double eps);
double trapzd(double (*func)(double), double a, double b, int n);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double zbrent(double (*func)(double), double x1, double x2, double tol);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
double expint(int n, double x);
double *vector(int i, int j);
void free_vector(double *a, int i, int j);
double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2, double xacc);
double golden(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);
/* trans_phot.c */
int trans_phot(WindPtr w, PhotPtr p, int iextract);
/* phot_util.c */
int stuff_phot(PhotPtr pin, PhotPtr pout);
int move_phot(PhotPtr pp, double ds);
int comp_phot(PhotPtr p1, PhotPtr p2);
double ds_to_cone(ConePtr cc, struct photon *p);
double ds_to_sphere(double r, struct photon *p);
double ds_to_sphere2(double x[], double r, struct photon *p);
int quadratic(double a, double b, double c, double r[]);
double ds_to_plane(struct plane *pl, struct photon *p);
double ds_to_closest_approach(double x[], struct photon *p, double *impact_parameter);
/* resonate.c */
double calculate_ds(WindPtr w, PhotPtr p, double tau_scat, double *tau, int *nres, double smax, int *istat);
int select_continuum_scattering_process(double kap_cont, double kap_es, double kap_ff, PlasmaPtr xplasma);
double kappa_bf(PlasmaPtr xplasma, double freq, int macro_all);
int kbf_need(double fmin, double fmax);
double sobolev(WindPtr one, PhotPtr p, double den_ion, struct lines *lptr, double dvds);
int doppler(PhotPtr pin, PhotPtr pout, double v[], int nres);
int scatter(PhotPtr p, int *nres, int *nnscat);
/* radiation.c */
int radiation(PhotPtr p, double ds);
double kappa_ff(PlasmaPtr xplasma, double freq);
double sigma_phot(struct photoionization *x_ptr, double freq);
double sigma_phot_topbase(struct topbase_phot *x_ptr, double freq);
double den_config(PlasmaPtr xplasma, int nconf);
/* wind_updates2d.c */
int wind_update(WindPtr (w));
int wind_rad_init(void);
int wind_rad_summary(WindPtr w, char filename[], char mode[]);
/* windsave.c */
int wind_save(WindPtr w, char filename[]);
int wind_read(char filename[]);
int wind_complete(WindPtr w);
/* extract.c */
int extract(WindPtr w, PhotPtr p, int itype);
int extract_one(WindPtr w, PhotPtr pp, int itype, int nspec);
/* pdf.c */
int pdf_init(void);
int pdf_gen_from_func(PdfPtr pdf, double (*func)(double), double xmin, double xmax, int njumps, double jump[]);
int pdf_gen_from_array(PdfPtr pdf, double x[], double y[], int n_xy, double xmin, double xmax, int njumps, double jump[]);
double pdf_get_rand(PdfPtr pdf);
int pdf_limit(PdfPtr pdf, double xmin, double xmax);
double pdf_get_rand_limit(PdfPtr pdf);
int pdf_to_file(PdfPtr pdf, char filename[]);
int pdf_check(PdfPtr pdf);
int recalc_pdf_from_cdf(PdfPtr pdf);
/* roche.c */
int binary_basics(void);
double ds_to_roche_2(PhotPtr p);
int hit_secondary(PhotPtr p);
double pillbox(PhotPtr p, double *smin, double *smax);
double phi(double s);
double dphi_ds(double s);
double d2phi_ds2(double s);
double roche_width(double x);
double roche2_width_max(void);
void roche(double s, double *value, double *derivative);
void roche_deriv(double s, double *value, double *derivative);
/* random.c */
int randvec(double a[], double r);
int randvcos(double lmn[], double north[]);
double vcos(double x);
/* stellar_wind.c */
int get_stellar_wind_params(void);
double stellar_velocity(double x[], double v[]);
double stellar_rho(double x[]);
int stellar_vel_grad(double x[], double velgrad[][3]);
/* proga.c */
int get_proga_wind_params(void);
int get_proga(void);
double proga_velocity(double x[], double v[]);
double proga_rho(double x[]);
/* corona.c */
int get_corona_params(void);
double corona_velocity(double x[], double v[]);
double corona_rho(double x[]);
/* knigge.c */
int get_knigge_wind_params(void);
double kn_velocity(double x[], double v[]);
double kn_rho(double x[]);
double kn_vzero(double r);
double kn_wind_mdot_integral(double r);
double kn_rho_zero(double r);
/* thierry.c */
int get_thierry_params(void);
int read_thierry(char infile[]);
int thierry_read_set(FILE *fptr, double array[], int nd);
double xthierry_velocity(double x[], double v[]);
double thierry_velocity(double x[], double v[]);
int thierry_vel_grad(double x[], double velgrad[][3]);
double thierry_rho(double x[]);
double xthierry_rho(double x[]);
/* disk.c */
double tdisk(double m, double mdot, double r);
double teff(double t, double x);
double gdisk(double mass, double mdot, double rmin);
double geff(double g0, double x);
double vdisk(double x[], double v[]);
double zdisk(double r);
double ds_to_disk(struct photon *p, int miss_return);
void disk_deriv(double s, double *value, double *derivative);
/* lines.c */
double total_line_emission(WindPtr one, double f1, double f2);
double lum_lines(WindPtr one, int nmin, int nmax);
int lum_pdf(PlasmaPtr xplasma, double lumlines);
double q21(struct lines *line_ptr, double t);
double q12(struct lines *line_ptr, double t);
double a21(struct lines *line_ptr);
double two_level_atom(struct lines *line_ptr, PlasmaPtr xplasma, double *d1, double *d2);
double line_nsigma(struct lines *line_ptr, PlasmaPtr xplasma);
double scattering_fraction(struct lines *line_ptr, PlasmaPtr xplasma);
double p_escape(struct lines *line_ptr, PlasmaPtr xplasma);
int line_heat(PlasmaPtr xplasma, PhotPtr pp, int nres);
/* continuum.c */
double one_continuum(int spectype, double t, double g, double freqmin, double freqmax);
double emittance_continuum(int spectype, double freqmin, double freqmax, double t, double g);
/* emission.c */
double wind_luminosity(double f1, double f2);
double total_emission(WindPtr one, double f1, double f2);
double adiabatic_cooling(WindPtr one, double t);
int photo_gen_wind(PhotPtr p, double weight, double freqmin, double freqmax, int photstart, int nphot);
double one_line(WindPtr one, double freqmin, double freqmax, int *nres);
double total_free(WindPtr one, double t_e, double f1, double f2);
double ff(WindPtr one, double t_e, double freq);
double one_ff(WindPtr one, double f1, double f2);
/* recomb.c */
double fb_verner_partial(double freq);
double fb_topbase_partial(double freq);
double integ_fb(double t, double f1, double f2, int nion, int fb_choice);
double total_fb(WindPtr one, double t, double f1, double f2);
double one_fb(WindPtr one, double f1, double f2);
int num_recomb(PlasmaPtr xplasma, double t_e);
double fb(PlasmaPtr xplasma, double t, double freq, int ion_choice, int fb_choice);
int init_freebound(double t1, double t2, double f1, double f2);
double get_nrecomb(double t, int nion);
double get_fb(double t, int nion, int narray);
double xinteg_fb(double t, double f1, double f2, int nion, int fb_choice);
int fb_save(char filename[]);
int fb_read(char filename[]);
/* diag.c */
int open_diagfile(void);
/* sv.c */
int get_sv_wind_params(void);
double sv_velocity(double x[], double v[]);
double sv_rho(double x[]);
double sv_find_wind_rzero(double p[]);
int sv_zero_init(double p[]);
double sv_zero_r(double r);
double sv_theta_wind(double r);
double sv_wind_mdot_integral(double r);
/* ionization.c */
int ion_abundances(PlasmaPtr w, int mode);
int convergence(PlasmaPtr xplasma);
int check_convergence(void);
int one_shot(PlasmaPtr xplasma, int mode);
double calc_te(PlasmaPtr xplasma, double tmin, double tmax);
double zero_emit(double t);
/* ispy.c */
int ispy_init(char filename[], int icycle);
int ispy_close(void);
int ispy(PhotPtr p, int n);
/* levels.c */
int levels(PlasmaPtr xplasma, int mode);
/* gradv.c */
double dvwind_ds(PhotPtr p);
int dvds_ave(void);
/* reposition.c */
int reposition(WindPtr w, PhotPtr p);
/* anisowind.c */
int randwind(PhotPtr p, double lmn[3], double north[3]);
double vrandwind(double x);
double reweightwind(PhotPtr p);
int make_pdf_randwind(double tau);
/* util.c */
int fraction(double value, double array[], int npts, int *ival, double *f);
int linterp(double x, double xarray[], double yarray[], int xdim, double *y);
int coord_fraction(int ichoice, double x[], int ii[], double frac[], int *nelem);
int where_in_2dcell(int ichoice, double x[], int n, double *fx, double *fz);
int wind_n_to_ij(int n, int *i, int *j);
int wind_ij_to_n(int i, int j, int *n);
/* density.c */
double get_ion_density(PhotPtr p, int nion);
/* detail.c */
int detailed_balance(PlasmaPtr xplasma, int nelem, double newden[]);
int rebalance(double rates_up[], double rates_down[], double fraction[], int ntot);
int wind_update_after_detailed_balance(PlasmaPtr xplasma, int nelem, double newden[]);
/* bands.c */
int init_bands(double t, double f1, double f2, int mode, struct xbands *band);
/* time.c */
double timer(void);
/* matom.c */
int matom(PhotPtr p, int *nres, int *escape);
double b12(struct lines *line_ptr);
double alpha_sp(struct topbase_phot *cont_ptr, PlasmaPtr xplasma, int ichoice);
double alpha_sp_integrand(double freq);
int kpkt(PhotPtr p, int *nres, int *escape);
int fake_matom_bb(PhotPtr p, int *nres, int *escape);
int fake_matom_bf(PhotPtr p, int *nres, int *escape);
int macro_pops(PlasmaPtr xplasma, double xne);
int macro_gov(PhotPtr p, int *nres, int matom_or_kpkt, int *which_out);
double get_kpkt_f(void);
double get_matom_f(void);
int photo_gen_kpkt(PhotPtr p, double weight, int photstart, int nphot);
int photo_gen_matom(PhotPtr p, WindPtr w, double weight, int photstart, int nphot);
int emit_matom(WindPtr w, PhotPtr p, int *nres, int upper);
double q_ioniz(struct topbase_phot *cont_ptr, double electron_temperature);
double q_recomb(struct topbase_phot *cont_ptr, double electron_temperature);
/* estimators.c */
int bf_estimators_increment(WindPtr one, PhotPtr p, double ds);
int bb_estimators_increment(WindPtr one, PhotPtr p, double tau_sobolev, double dvds, int nn);
int mc_estimator_normalise(int n);
double total_fb_matoms(PlasmaPtr xplasma, double t_e, double f1, double f2);
double total_bb_cooling(PlasmaPtr xplasma, double t_e);
double macro_bb_heating(PlasmaPtr xplasma, double t_e);
double macro_bf_heating(PlasmaPtr xplasma, double t_e);
int bb_simple_heat(PlasmaPtr xplasma, PhotPtr p, double tau_sobolev, double dvds, int nn);
double get_gamma(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double gamma_integrand(double freq);
double get_gamma_e(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double gamma_e_integrand(double freq);
double get_alpha_st(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double alpha_st_integrand(double freq);
double get_alpha_st_e(struct topbase_phot *cont_ptr, PlasmaPtr xplasma);
double alpha_st_e_integrand(double freq);
/* wind_sum.c */
int xtemp_rad(WindPtr w);
/* yso.c */
int get_yso_wind_params(void);
double yso_velocity(double x[], double v[]);
double yso_rho(double x[]);
/* cylindrical.c */
double cylind_ds_in_cell(PhotPtr p);
int cylind_make_grid(WindPtr w);
int cylind_wind_complete(WindPtr w);
int cylind_volumes(WindPtr w);
int cylind_where_in_grid(double x[]);
int cylind_get_random_location(int n, double x[]);
int cylind_extend_density(WindPtr w);
/* rtheta.c */
double rtheta_ds_in_cell(PhotPtr p);
int rtheta_make_grid(WindPtr w);
int rtheta_wind_complete(WindPtr w);
int rtheta_volumes(WindPtr w);
int rtheta_where_in_grid(double x[]);
int rtheta_get_random_location(int n, double x[]);
int rtheta_extend_density(WindPtr w);
/* spherical.c */
double spherical_ds_in_cell(PhotPtr p);
int spherical_make_grid(WindPtr w);
int spherical_wind_complete(WindPtr w);
int spherical_volumes(WindPtr w);
int spherical_where_in_grid(double x[]);
int spherical_get_random_location(int n, double x[]);
int spherical_extend_density(WindPtr w);
/* cylind_var.c */
double cylvar_ds_in_cell(PhotPtr p);
int cylvar_make_grid(WindPtr w);
int cylvar_wind_complete(WindPtr w);
int cylvar_volumes(WindPtr w);
int cylvar_where_in_grid(double x[], int ichoice, double *fx, double *fz);
int cylvar_get_random_location(int n, double x[]);
int cylvar_extend_density(WindPtr w);
int cylvar_coord_fraction(int ichoice, double x[], int ii[], double frac[], int *nelem);
/* bilinear.c */
int bilin(double x[], double x00[], double x01[], double x10[], double x11[], double *f, double *g);
int xquadratic(double a, double b, double c, double r[]);
/* gridwind.c */
int create_maps(int ichoice);
int calloc_wind(int nelem);
int calloc_plasma(int nelem);
int check_plasma(PlasmaPtr xplasma, char message[]);
int calloc_macro(int nelem);
/* py_wind_sub.c */
int zoom(int direction);
int overview(WindPtr w, char rootname[]);
int position_summary(WindPtr w);
int abs_summary(WindPtr w, char rootname[], int ochoice);
int adiabatic_cooling_summary(WindPtr w, char rootname[], int ochoice);
int lum_summary(WindPtr w, char rootname[], int ochoice);
int photo_summary(WindPtr w, char rootname[], int ochoice);
int recomb_summary(WindPtr w, char rootname[], int ochoice);
int electron_summary(WindPtr w, char rootname[], int ochoice);
int freq_summary(WindPtr w, char rootname[], int ochoice);
int nphot_summary(WindPtr w, char rootname[], int ochoice);
int temp_summary(WindPtr w, char rootname[], int ochoice);
int temp_rad(WindPtr w, char rootname[], int ochoice);
int weight_summary(WindPtr w, char rootname[], int ochoice);
int velocity_summary(WindPtr w, char rootname[], int ochoice);
int mo_summary(WindPtr w, char rootname[], int ochoice);
int vol_summary(WindPtr w, char rootname[], int ochoice);
int wind_element(WindPtr w);
int tau_h_summary(WindPtr w, char rootname[], int ochoice);
int coolheat_summary(WindPtr w, char rootname[], int ochoice);
int complete_file_summary(WindPtr w, char root[], int ochoice);
int wind_reg_summary(WindPtr w, char rootname[], int ochoice);
int dvds_summary(WindPtr w, char rootname[], int ochoice);
/* py_ion.c */
int ion_summary(WindPtr w, int element, int istate, int iswitch, char rootname[], int ochoice);
int tau_ave_summary(WindPtr w, int element, int istate, double freq, char rootname[], int ochoice);
int line_summary(WindPtr w, int element, int istate, char rootname[], int ochoice);
int total_emission_summary(WindPtr w, char rootname[], int ochoice);
int modify_te(WindPtr w, char rootname[], int ochoice);
int partial_measure_summary(WindPtr w, int element, int istate, char rootname[], int ochoice);
/* py_write.c */
int write_array(char filename[], int choice);
int display(char name[]);
