

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"

/* This is a routine to calculate some of the constants that are used in python */

int
main ()
{
  double x;
  double n;
  printf ("hello world\n");

  printf
    ("This is a constant which relates when collisions and radiation compete in photoionzation\n");

  x = MELEC / (6. * PI * BOLTZMANN);
  x = sqrt (x);
  x *= (E * E * 1e-15 / H);

  printf ("at 1000 Angstroms the constant is %e\n", x);





  printf (" This is the bremsstrahlung constant \n");


  x = pow (E, 6.0) / (C * C * C) / (MELEC * MELEC) * sqrt (MELEC / BOLTZMANN);

  x *= 4. * PI * 8. / 3. * sqrt (2. * PI / 3);

  printf ("x %e\n", x);

  printf ("This is the KL crossection for H at threshold\n");

  x = 64 / pow (3., 1.5) * PI * BOHR * BOHR * ALPHA;

  printf ("x %e\n", x);


  printf ("This is the fb constant\n");

  /* H **4/ m_e **3 /c**2  * 4 / sqrt(PI) *(m_e/2 k) ** 1.5 */

  x = H * H * H * H / (MELEC * MELEC * MELEC) / (C * C);
  x *= 4. / sqrt (PI);
  x *= pow (MELEC / (2. * BOLTZMANN), 1.5);


  printf ("x %e\n", x);



  printf ("The fb constant a second try\n");

  x = 2. * PI * PI * E * E * E * E * MELEC / (H * H);

  printf ("Ionization energy of H %g\n", x);

  x =
    64 * sqrt (PI) * E * E * E * E * H / (pow (3., 1.5) * MELEC * MELEC * C *
					  C * C) * pow (x, 1.5);

  printf ("A new fb constant %g\n", x);


  x = (2 * PI * MELEC * BOLTZMANN);
  x = pow (x, 1.5) / (H * H * H);

  printf ("xsaha %g 2x %g \n", x, 2.0 * x);


// Calculate a new FB emissivity constant for use with Ferland emissivity

  x = 2. * PI * MELEC * BOLTZMANN / (H * H);
  x = pow (x, -1.5);
  x *= 8. * PI / (C * C);
  x *= H;


  printf (" Ferland FB emissivity constant %g\n", x);

  for (n = 20; n < 50; n++)
    printf ("exp %g\n", exp (-n));
  return EXIT_SUCCESS;
}
