#include <stdio.h>

main ()
{

  char hello[] =
    "2.8019E-093.2268E-091.2423E-092.5084E-091.6263E-095.5754E-099.9628E-099.6235E-09";
  float foo, goo;

  sscanf (hello, "%10e%10e", &foo, &goo);

  printf ("foo %e %e\n", foo, goo);

}
