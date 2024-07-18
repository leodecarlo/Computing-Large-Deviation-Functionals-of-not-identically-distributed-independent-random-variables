/* stat.hd: wick's elementary statistics routines */
/* include wdefs.h, linalg.h */
/* doubled version */

long random();
double ran0(void)
{
  double f;long r;
  r = random(); /* make sure to give srandom(SEED) in main */
  if(r == 0){
    printf("WARNING: ran0() delivered exact 0! Trying again ...");
    r = random(); 
    r = random(); 
    r = random(); 
    if(r == 0) nrerror("ERROR in ran0(): still delivered exact 0."); else;
   } else;
  f = r/2147483649.0;
  return f;
}

long random();
void srandom(int SEED)
{ int k;long dum;
  for(k=1;k<=SEED;k++) dum = random();
} 

int ransign()
{ int ret;
  double r;
  r = ran0();
  if(r < 0.5) ret = 1;
  else ret = -1;
  return ret;
}

/* randint: generates random integers in [first,last] */
int randint(int first,int last)
{
   double r;int ret;
   r = ran0();
   ret = first + r*(last - first + 1.0);
   if( (ret<first)||(ret>last) ) nrerror("error in randint:outside range");else;
   return ret;
}

double ran1(void) /*NumR in C: supposed to break up local correlations*/
{
  double f;int i,k;static int first=1;
  static double vec[100];
   k = randint(0,99);
  if( first==1) for(i=0;i<=99;i++) vec[i] = ran0();
   first++;
  f = vec[k];vec[k] = ran0();
  return f;
}

/*double ran1(void)
{
  double f;int i,k;
  char rt[] = "ran1";
  struct vector vec;
  invec(&vec,1,100,"vec");
   k = randint(1,100);
  for(i=1;i<=100;i++)
   Avec(&vec,i,ran0(),rt);
  f = Evec(&vec,k,rt);
  freevec(&vec);
  return f;
}*/
     
double uniform(double low,double high)
{ double ret,r;
  r = ran0();
ret = r*low + (1.0-r)*high;
return ret;
}

double normal(void)
  {
    static int iset=0;
    static float gset;
    double fac,r,v1,v2;

    if(iset==0){
    do{
    v1 = 2.0*ran0() -1.0;
    v2 = 2.0*ran0() -1.0;
    r = v1*v1 + v2*v2;
    } while( r >= 1.0);

    if(r <= 0.0) nrerror("error in normal");

    fac = sqrt(-2.0*log(r)/r);
    gset = v1*fac;
    iset=1;
    return v2*fac;
    }

    else{
      iset=0;
      return gset;
      }

 }

void normalvec(struct vector *V)
{int i;
  if(V->type != 'v') nrerror("type error in normalvec");else;
   for(i=V->in;i<=V->fin;i++) V->vec[i] = normal();
}

/* gammavec: creats a vector of independent gamma's of given mean and index. */
void gammavec(struct vector *V,int n,double mean)
{
 int i,j;double mn,r,dum,sum;
 char rt[] = "gammavec";
 mn = mean/n;
 if(V->type != 'v') nrerror("Type error in gammavec.");else;
 for(i=V->in;i<=V->fin;i++){
   sum = 0.0;
   for(j=1;j<=n;j++){
     r = ran0(); 
     if(r > 0.0) dum = - log(r);
     else dum = 1.0e+30;
     sum = sum + dum;
   }
   Avec(V,i,sum,rt);
 }
}


void uniform_vec(struct vector *V,struct vector *low,struct vector *high)
{
  int i;
  if( (V->type != 'v')||(high->type != 'v')||(low->type != 'v') )
   nrerror("type error in normalvec");else;
  if( (V->in != low->in)||(V->fin != low->fin) )
     nrerror("index error 1 in normalvec");else;
   if( (V->in != high->in)||(V->fin != high->fin) )
     nrerror("index error 2 in normalvec");else;

    for(i=V->in;i<=V->fin;i++){
     double r;r = ran0();
     V->vec[i] = r*low->vec[i] + (1.0-r)*high->vec[i];
    }
}

void u_vec(struct vector *V,double a,double b)
{
 int i;
  if(V->type != 'v') nrerror("type error in u_vec");else;
   for(i=V->in;i<=V->fin;i++) V->vec[i] = a + (b - a)*ran0();
}

double stdev(struct vector *Vec)
{
  int i,length;
  double mean,std,fac;
  char rt[] = "stdev";
  if(Vec->type != 'v') nrerror("type error in stdev");else;
  length = Vec->fin - Vec->in + 1;
  if(length <= 1) nrerror("Error in stdev: only one entry??");else;
  mean = vecsum(Vec,Vec->in,Vec->fin)/length;
  std=0.0;
  for(i=Vec->in;i<=Vec->fin;i++){
   std = std + Evec(Vec,i,rt)*Evec(Vec,i,rt);
   }
   fac = length/(length - 1.0);
   std = std/length - mean*mean;
   if(std < 0.0) nrerror("Error in stdev: std < 0??"); else;
   std = sqrt(fac*std);
   return std;
 }

double t_test(struct vector *y1,struct vector *y2)
{
  int n1,n2;
  double y1bar,y2bar,sp,s1,s2,ret;
  if(y1->type != 'v') nrerror("type error in y1; routine t_test");else;
  if(y2->type != 'v') nrerror("type error in y2; routine t_test");else;
  n1 = y1->fin - y1->in + 1;
  n2 = y2->fin - y2->in + 1;
  y1bar = vecsum(y1,y1->in,y1->fin)/n1;
  y2bar = vecsum(y2,y2->in,y2->fin)/n2;
  s1 = stdev(y1);
  s2 = stdev(y2);
  if(n1 + n2 -2 > 2) sp = sqrt( ( (n1 - 1)*s1*s1 + (n2 - 1)*s2*s2 )/( n1 + n2 - 2) );
  else nrerror("Error in t_test: d.f. < 3??"); 
  ret = (y1bar - y2bar)/(sp*sqrt( 1.0/n1 + 1.0/n2 ) );
  return ret;
}


void stats(struct vector *Vec,double *Mean,double *standard_dev)
{
  double mean,std;int i,length;
  if(Vec->type != 'v') nrerror("type error in stdev");else;
  length = Vec->fin - Vec->in + 1;
  if(length <= 0) nrerror("error1 in stdev");
  mean = vecsum(Vec,Vec->in,Vec->fin)/length;
  std=0.0;
  for(i=Vec->in;i<=Vec->fin;i++){
   std = std + Vec->vec[i]*Vec->vec[i];
   }
   std = std/length - mean*mean;
   if(std < 0.0) nrerror("Error in stats:error 2 in stdev");else;
   std = sqrt(std);
   (*Mean) = mean;
   (*standard_dev) = std;
}

void stats_alt(struct vector *Vec,double *Mean,double *Var)
{
  double mean,std;int i,length;
  if(Vec->type != 'v') nrerror("type error in stdev");else;
  length = Vec->fin - Vec->in + 1;
  if(length <= 0) nrerror("error1 in stdev");
  mean = vecsum(Vec,Vec->in,Vec->fin)/length;
  std=0.0;
  for(i=Vec->in;i<=Vec->fin;i++){
   std = std + Vec->vec[i]*Vec->vec[i];
   }
   std = std/length - mean*mean;
   if(std < 0.0) nrerror("Error in stats_alt:error 2 in stdev"); else;
   (*Mean) = mean;
   (*Var) = std;
}


double variance(struct vector *v)
{
  double mean,ret;int length;
  char rt[] = "variance";
  if( v->type != 'v' )  nrerror("type error in variance");else;
  length = v->fin - v->in + 1;
  if(length <= 0) nrerror("Error in variance: wrong length.");else;
  mean = vecsum(v,v->in,v->fin)/length;
  ret = dotprod(v,v)/length - mean*mean;
  return ret;
 }


/* "pvar": variance of a vector, v, relative to a probability
vector, p:*/
double pvar(struct vector *p,struct vector *v)
{
  double dum,mean,ret;int i;
  struct vector vv;
  char rt[] = "pvar";
  if( (v->type != 'v')||(p->type != 'v') )  nrerror("type error in pvar");else;
  invec(&vv,v->in,v->fin,"vv");
  mean = dotprod(p,v);
  for(i=v->in;i<=v->fin;i++){
   dum = Evec(v,i,rt);
   Avec(&vv,i,dum*dum,rt);
  }
  ret = dotprod(p,&vv) - mean*mean;
  freevec(&vv);
  return ret;
 }


/* irandlist: stores random integers in [first,last] */
 
int irandlist(int first,int last,struct ivector *rlist)
 {double r;int j;

  if(rlist->type != 'V') nrerror("type error in irandlist");else;

   for(j=rlist->in;j<=rlist->fin;j++){
   r = ran0();
   rlist->vec[j] = first + r*(last - first + 1.0);
   if( (rlist->vec[j]<first)||(rlist->vec[j]>last) ) nrerror("error in irandlist: outside range");else;
   }

   return 1;
}

double correlation(struct vector *V,struct vector *W)
{
   double ret,x,y,z,sv,sw,denom;int N;
   if( (V->in != W->in)||(V->fin != W->fin) ) 
   nrerror("index error in correlation");else;
   N = V->fin - V->in + 1;
   sv = vecsum(V,V->in,V->fin);
   sw = vecsum(W,W->in,W->fin);
   x = dotprod(V,V) - (sv*sv)/N;
   y = dotprod(W,W) -  (sw*sw)/N; 
   z = dotprod(W,V) -  (sv*sw)/N; 
   denom = x*y;
   if( (x > 0.0)&&(y > 0.0) ) denom = sqrt(denom); 
   else nrerror("error in correlation: a vector is constant");
   ret = z/denom;
return ret;
}

double entropy(struct vector *rho)
{
 int i;double dum,sum;
 char rt[] = "entropy";
  if(rho->type != 'v') nrerror("Type error in entropy.");else;
  sum = 0.0;
  for(i = rho->in;i <= rho->fin;i++){
      dum = Evec(rho,i,rt);
      if( dum > 1.0e-30 ) sum = sum - dum*log( dum );  
      else{
       if( dum > 0.0);
       else nrerror("Error in entropy: a component is negative??");
      }
  }
 return sum;
}

double iProb(struct ivector *data,int n)
{
 int i,count,denom;double ret;
 char rt[] = "iProb";
  if(data->type != 'V') nrerror("Type error in iProb.");else;
  count=0;denom=0;
  for(i=data->in;i<=data->fin;i++){
   if( Eivec(data,i,rt) >= n) count++;else;
   denom++;
  }
  ret = (1.0*count)/denom;
  return ret;
}

int rand_int_vec(int RNG,struct ivector *kvec,struct vector *pvec)
{
  int ret,i;
  double s,rr;
  char rt[] = "rand_int_vec";
  if(pvec->type != 'v') nrerror("Type error in vector pvec in routine rand_int_vec");else;
  if(kvec->type != 'V') nrerror("Type error in vector kvec in routine rand_int_vec");else;
  if( (pvec->in != kvec->in)||(pvec->fin != kvec->fin) ) nrerror("Index error in rand_int_vec."); else;
  s = Evec(pvec,pvec->in,rt);
  if(RNG == 0) rr = ran0(); else rr = ran1();
  i = pvec->in;
  ret = Eivec(kvec,i,rt);
  while(i < pvec->fin){
   if(rr < s ){
     ret = Eivec(kvec,i,rt);
     i = pvec->fin + 1;
   }
   else{
     i++;
     s = s + Evec(pvec,i,rt);
     if(i == pvec->fin) ret = Eivec(kvec,i,rt);
   }
  }
 return ret;
}

/* binomial: executes N choices of a 'type' from some range, according to probabilities in vector p_types,
and stores the numbers of each hit in ivector n_types. */
void binomial(int RNG,struct ivector *n_types,struct vector *p_types,int N)
{
  int k,i,dum,kch;
  char rt[] = "binomial";
  struct ivector kvec;
  if(p_types->type != 'v') nrerror("Type error in vector p_types in routine binomial");else;
  if(n_types->type != 'V') nrerror("Type error in vector n_types in routine binomial");else;
  if( (p_types->in != n_types->in)||(p_types->fin != n_types->fin) ) nrerror("Index error in binomial."); else;
  inivec(&kvec,n_types->in,n_types->fin,"kvec");
  for(k = n_types->in;k <= n_types->fin;k++) Aivec(&kvec,k,k,rt);
  setivec(n_types,n_types->in,n_types->fin,0);
  for(i = 1;i <= N;i++){
   kch = rand_int_vec(RNG,&kvec,p_types);
   dum = Eivec(n_types,kch,rt);
   Aivec(n_types,kch,dum + 1,rt);
  }
  freeivec(&kvec);
}

/* B delivers a Bernouli (no. of successes in N trials) */
int B(int N, double p)
{ 
   int i,ret;double rr;
   ret = 0;
  for(i = 1;i <= N;i++){
   rr = ran0();
   if(rr < p) ret++;else;
  }
  return ret;
}


/* Poisson delivers a random integer */
int Poisson(double mean)
{
  int ret,i;
  double fac,dum,sum;
  ret = 200;
  if(mean < 40.0){
    fac = exp( mean )*ran0();
    dum = 1.0;sum = 1.0;i = 0;
    while(i <= 200){
     if(fac < sum){
       ret = i;i = 201;
     }
     else{
      dum = (dum*mean)/(i+1);
      sum = sum + dum;
      i++;
     }
    }
    if(ret >= 200) nrerror("Error in Poisson: mean < 40, yet return >= 200??");else;
  }
  else{
   fac = mean + sqrt(mean)*normal();
   ret = fac;
   if( ret < 0 ) nrerror("Error in Poisson: in Gaussian approximation, ret would be negative."); else;
  }
  return ret;
}


int geometric(double inf,int stop)
{
   int j,ret,hit;
   double rr;
   j = 1;hit = 0;
   while((j < stop)&&(hit == 0)){
    rr = ran0();
    if(rr <= inf){
     ret = j;hit++;
    }
    else j++;
   }
  if(j >= stop) nrerror("Error in geometric: exceeded stop."); else;
  return ret;
}
    

 /* stores a random list of exactly K ones and N-K zeroes, for use e.g. in constructing a random Kauffman graph */
/* void rand_indicator(int N,int K,struct ivector *vec)
 {
char rt[] = "rand_indicator";
 int i,j,dum,dum1;
  if(vec->type != 'V') nrerror("Type error in routine rand_indicator, in ivector vec.");else;
  if( (vec->in != 1)||(vec->fin != N) )
   nrerror("Index error in routine rand_indicator.");else;

 setivec(vec,1,N,0);
 for(i=1;i<=K;i++){
 Aivec(vec,i,1,rt);
 }
 for(i=1;i<=N;i++){
  j = randint(1,N);
  dum = Eivec(vec,i,rt);
  dum1 = Eivec(vec,j,rt);
  Aivec(vec,j,dum,rt); 
  Aivec(vec,i,dum1,rt); 
 }
}*/

/*
void sampleworeplacement(struct vector *source, struct vector *sample)
{
 char rt[] = "sampleworeplacement";
 int n,k;double dum;

  if( (source->type != 'v')||(sample->type != 'v') )  nrerror("type error in sampleworeplacement");else;
  if( (source->in != sample->in)||(sample->fin != source->fin) ) nrerror("index error in sampleworeplacement");else;

  for(n = source->in;n <= source->fin;n++){
   k = randint(source->in,source->fin);
   dum = Evec(source,k,rt);
   Avec(sample,n,dum,rt);
  }
} */
