gsl_rng *random_setup(void)
{
const gsl_rng_type *TT;
long seedy;
srand((unsigned) time(NULL));
//seedy = -rand();
//seedy = time(NULL)*(int)getpid();  // This is the different line!
seedy = time(NULL)*getpid();  // This is the different line!

gsl_rng_env_setup ();
gsl_rng_default_seed = seedy;
TT = gsl_rng_default;

return gsl_rng_alloc(TT);
}

// -------------------------------------------------------------------- //

int ipow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = ipow(x, p/2);
  if (p%2 == 0)             return tmp * tmp;
  else                      return x * tmp * tmp;
}
