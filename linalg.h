/* new_linalg.hd; contains wick's linear algebra routines */
/* Assumes inclusion of wdefs.h */
/* "Short" version, containing frequently-used routines only */
/* doubled version */
/* 5/1/01: partly changed-over to new defs. */


 #define swap(a,b) {double temp; temp = (a);(a) = (b);(b) = temp; } 
 #define iswap(a,b) {int temp; temp = (a);(a) = (b);(b) = temp; } 

void vswap(struct vector *v,int i,int j)
{
  double dum,dum1;
   char rt[] = "vecswap";
     dum = Evec(v,i,rt);
     dum1 = Evec(v,j,rt);
     Avec(v,i,dum1,rt);
     Avec(v,j,dum,rt);
}

void ivswap(struct ivector *v,int i,int j)
{
  int dum,dum1;
   char rt[] = "vecswap";
     dum = Eivec(v,i,rt);
     dum1 = Eivec(v,j,rt);
     Aivec(v,i,dum1,rt);
     Aivec(v,j,dum,rt);
}


void vec_inc(struct vector *V,int index,double value)
{  
   char rt[] = "vec_inc";
   if(V->type != 'v') nrerror("Type error in vec_inc.");else;
   Avec(V,index,Evec(V,index,rt) + value,rt);
}

void ivec_inc(struct ivector *V,int index,int value)
{  
   char rt[] = "ivec_inc";
   if(V->type != 'V') nrerror("Type error in ivec_inc.");else;
   Aivec(V,index,Eivec(V,index,rt) + value,rt);
}

int ivec_check(struct ivector *V,int value)
{  int i,ret; 
   char rt[] = "ivec_check";
   if(V->type != 'V') nrerror("Type error in ivec_check.");else;
   ret = 1;
   for(i = V->in;i <= V->fin;i++){
     if(Eivec(V,i,rt) >= value); else ret--;
   }
   return ret;
}

/* vec_over tests for over-flow dangers: */
void vec_over(struct vector *V,double Over,char routine[],char message[])
{ 
 int j;
 char rt[] = "vec_over";
 if(V->type != 'v') nrerror("Type error in vec_over.");
 for(j=V->in;j<=V->fin;j++){
    if( fabs(Evec(V,j,rt)) > Over ){
      printf("Possible overflow in routine %s;\n",routine);
      printf("In vector %s, entry %i;\n",V->ident,j);
      printf("|component| > %g.\n",Over);
      printf("%s.\n",message);
      nrerror("Error reported by vec_over.");
     } else;
 }
}

double vecfabs(struct vector *V)
{ 
 int j;double sum;
 char rt[] = "vecfabs";
   if(V->type != 'v') nrerror("type error in vecfabs");
   sum = 0.0;
   for(j=V->in;j<=V->fin;j++){
    sum = sum + fabs( Evec(V,j,rt) );
   }
  return sum;
}

double vecmax(struct vector *V)
{ 
 int j;double ret,dum;
 char rt[] = "vecmax";
   if(V->type != 'v') nrerror("type error in vecfabs");else;
   ret = fabs( Evec(V,V->in,rt) );
   for(j = V->in + 1;j <= V->fin;j++){
    dum = fabs( Evec(V,j,rt) );
     if(dum > ret) ret = dum;else;
   }
  return ret;
}
double vecsize(struct vector *U,struct vector *scale)
{
  double m,temp,denom;int k;
 char rt[] = "vecsize";
 if( (U->in != scale->in)||
  (U->fin != scale->fin) ) nrerror("index error in vecsize");else;
   denom = Evec(scale,U->in,rt);
   if( denom > 0.0) temp  = fabs( Evec(U,U->in,rt)/denom );
   else nrerror("error in vecsize: scale <=0");
  for(k=U->in+1;k<=U->fin;k++){
   denom = Evec(scale,k,rt);
   if( denom > 0.0) temp  = fabs( Evec(U,k,rt)/denom );
   else nrerror("error in vecsize: scale <=0");
   if(temp > m) m = temp;else;
  }
  return m;
}

double vecerror(struct vector *U,struct vector *V,struct vector *scale)
{
  double m,temp,fac;int k;
  char rt[] = "vecerror";
  if( (U->in != V->in)||(U->fin != V->fin) )
   nrerror("index error in vecerror");else;
  if( (U->in != scale->in)||(U->fin != scale->fin) )
   nrerror("index error 2 in vecerror");else;
   if( fabs(Evec(scale,1,rt)) > 0.0) fac = 1.0/(Evec(scale,1,rt));
   else nrerror("error in vecerror: scale[] = 0?");
   m = fabs( ( Evec(U,1,rt) - Evec(V,1,rt) )*fac );
  for(k = U->in + 1;k <= U->fin;k++){
   if( fabs(Evec(scale,k,rt)) > 0.0) fac = 1.0/(Evec(scale,k,rt));
   else nrerror("error in vecerror: scale = 0?");
   temp = fabs( ( Evec(U,k,rt) - Evec(V,k,rt) )*fac );
   if(temp > m) m = temp;else;
  }
  return m;
}

int ivecerror(struct ivector *U,struct ivector *V)
{int k;int sum=0;
   if( V->type != 'V' ) nrerror("type error in ivecerror");else
   if( U->type != 'V' ) nrerror("type error in ivecerror");else
  if( (U->in != V->in)||(U->fin != V->fin) )
  nrerror("index error in ivecerror");else;
  for(k=U->in;k<=U->fin;k++){
  sum = sum + abs(U->vec[k] - V->vec[k] );
  }
  return sum;
}


double vecsum(struct vector *ptr,int first,int last)
{ int j;double sum =0.0;
  char rt[] = "vecsum";
   if(ptr->type != 'v') nrerror("type error in vecsum");
   if( (first < ptr->in)||(last > ptr->fin) ){
    printf("Index error in vecsum, on vector %s\n",ptr->ident); 
    nrerror("Index error in vecsum"); 
   } else;
   for(j=first;j<=last;j++){
    sum = sum + Evec(ptr,j,rt);
   }
  return sum;
}

int ivecsum(struct ivector *ptr,int first,int last)
{ int j;int sum =0;
   if(ptr->type != 'V') nrerror("type error in ivecsum");
   if( (first < ptr->in)||(last > ptr->fin) )
    nrerror("Index error in ivecsum");else;
   for(j=first;j<=last;j++){
    sum = sum + ptr->vec[j];
    }
  return sum;
}

double vecfrac(int greater,int equal,double level,struct vector *V)
{ int rep,count;double ret;
    count = 0;
   for(rep=V->in;rep<=V->fin;rep++){
    if(greater==1){
      if(equal==1){
         if(V->vec[rep] >= level) count++;else;
      }
      else{
         if(V->vec[rep] > level) count++;else;
      }
    }
    else{
      if(equal==1){
         if(V->vec[rep] <= level) count++;else;
      }
      else{
         if(V->vec[rep] < level) count++;else;
      }
    }

   }

  ret = (1.0*count)/(V->fin - V->in + 1);
  return ret;

}

void vecassign(struct vector *V,struct vector *V1)
{ int k; double dum;
 char rt[] = "vecassign";
  if( (V->type != 'v')||(V1->type != 'v') )
  nrerror("type error in vecassign");else;
  if( (V->in != V1->in)||(V->fin != V1->fin) )
    nrerror("index error in vecassign");else;
  for(k=V->in;k<=V->fin;k++){ 
    dum = Evec(V1,k,rt);
    Avec(V,k,dum,rt);
  }
}

void ivecassign(struct ivector *V,struct ivector *V1)
{ int k;int dum;
 char rt[] = "ivecassign";
  if( (V->type != 'V')||(V1->type != 'V') )
  nrerror("type error in ivecassign");else;
  if( (V->in != V1->in)||(V->fin != V1->fin) )
    nrerror("index error in ivecassign ");else;
  for(k=V->in;k<=V->fin;k++){ 
    dum = Eivec(V1,k,rt);
    Aivec(V,k,dum,rt);
  }
}

int setvec(struct vector *V,int first,int last,double x)
{int i;
 char rt[] = "setvec";
 if(V->type != 'v')  nrerror("type error in setvec");
  if( (first < V->in)||(last > V->fin) )
  nrerror("Index error in setvec"); else;
  for(i=first;i<=last;i++) Avec(V,i,x,rt);
  return 1;
}

int setivec(struct ivector *V,int first,int last,int x)
{int i;
  char rt[] = "setivec";
 if(V->type != 'V')  nrerror("type error in setivec");
  if( (first < V->in)||(last > V->fin) )
  nrerror("Index error in setivec"); else;
  for(i=first;i<=last;i++) Aivec(V,i,x,rt);
  return 1;
}

/* vecext: extracts Source[Sfirst],...,Source[Sfirst + num - 1]
and copies them into Target[Tfirst],...,Target[Tfirst + num - 1] */

void vecext(int Tfirst,struct vector *Target,int Sfirst,
struct vector *Source, int num)
{
    int i; double dum;
    char rt[] = "vecext";
    if( (Target->type != 'v')||(Source->type != 'v') )
      nrerror("type error in vecext");else;
    if( (Tfirst < Target->in)||(Tfirst+num-1 > Target->fin) )
      nrerror("Index error in Target, in vecext");else;
    if( (Sfirst < Source->in)||(Sfirst+num-1 > Source->fin) )
      nrerror("Index error in Source, in vecext");else;

    for(i=0;i <= num-1;i++){
      dum =  Evec(Source,Sfirst+i,rt);
      Avec(Target,Tfirst + i,dum,rt);
    }
}

void Vecext(int Tfirst,struct vector *Target,int Sfirst,
struct vector *Source, int num,char rt[],int dum)
{
    int i;
    if( (Target->type != 'v')||(Source->type != 'v') )
      nrerror("type error in vecext");else;
    if( (Tfirst < Target->in)||(Tfirst+num-1 > Target->fin) ){
      printf("Index error in target in vecext,in routine %s, number %i\n",rt,dum);
      nrerror("Index error in target in vecext");
    } else;
    if( (Sfirst < Source->in)||(Sfirst+num-1 > Source->fin) ){
      printf("Index error in source in vecext,in routine %s, number %i\n",rt,dum);
      nrerror("Index error in source in vecext");
     }else;

    for(i=0;i <= num-1;i++) Target->vec[Tfirst + i] = Source->vec[Sfirst+i];
}


void ivecext(int Tfirst,struct ivector *Target,int Sfirst,
struct ivector *Source, int num)
{
    int i;
    if( (Target->type != 'V')||(Source->type != 'V') )
      nrerror("type error in vecext");else;
    if( (Tfirst < Target->in)||(Tfirst+num-1 > Target->fin) )
      nrerror("Index error 1 in vecext");else;
    if( (Sfirst < Source->in)||(Sfirst+num-1 > Source->fin) )
      nrerror("Index error 2 in vecext");else;

    for(i=0;i <= num-1;i++) Target->vec[Tfirst + i] = Source->vec[Sfirst+i];
}


/* matext: extracts Source[Sfirstrow],...,Source[Sfirstrow + numrow - 1],
Source[Sfirstcol],...,Source[Sfirstcol + numcol - 1],
and copies them into Target[(same)] */

void matext(int Tfirstrow,int Tfirstcol,struct matrix *Target,
int Sfirstrow, int Sfirstcol,
struct matrix *Source, int numrow,int numcol)
{
    int i,j;double dum;
    char rt[] = "matext"; 
    if( (Target->type != 'm')||(Source->type != 'm') )
      nrerror("type error in matext");else;
     if( (Tfirstrow < Target->inrow)||(Tfirstrow+numrow-1 > Target->finrow) )
      nrerror("Index error 1 in matext");else;
    if( (Tfirstcol < Target->incol)||(Tfirstcol+numcol-1 > Target->fincol) )
      nrerror("Index error 2 in matext");else;
    if( (Sfirstrow < Source->inrow)||(Sfirstrow+numrow-1 > Source->finrow) )
      nrerror("Index error 3 in matext");else;
    if( (Sfirstcol < Source->incol)||(Sfirstcol+numcol-1 > Source->fincol) )
      nrerror("Index error 4 in matext");else;

    for(i=0;i <= numrow-1;i++){
    for(j=0;j <= numcol-1;j++){
      dum = Emat(Source,Sfirstrow+i,Sfirstcol+j,rt);
      Amat(Target,Tfirstrow+i,Tfirstcol+j,dum,rt);  
     }
     }
}


int setmat(struct matrix *M,double x)
{
  int i,j;
  char rt[] = "setmat";
  if(M->type != 'm')  nrerror("type error in setmat");
   for(i=M->inrow;i<=M->finrow;i++){
     for(j=M->incol;j<=M->fincol;j++){
      Amat(M,i,j,x,rt);
     }
   }
   return 1;
}

int setimat(struct imatrix *M,int x)
{
  int i,j;
  char rt[] = "setimat";
  if(M->type != 'M')  nrerror("type error in setmat");
   for(i=M->inrow;i<=M->finrow;i++){
     for(j=M->incol;j<=M->fincol;j++){
      Aimat(M,i,j,x,rt);
     }
   }
   return 1;
}

int Identity(struct matrix *M)
{
  int i;
  if(M->type != 'm')  nrerror("type error in setmat");else;
  if( (M->inrow != M->incol)||(M->finrow != M->fincol) )
   nrerror("index error in Identity");else;
   setmat(M,0.0);
   for(i=M->inrow;i<=M->finrow;i++){
      M->mat[i][i] = 1.0;
     }
     return 1;
}

void matassign(struct matrix *A,struct matrix *B)
{ int i,j;
  if( (A->type != 'm')||(B->type != 'm') )
  nrerror("type error in matasign");
  if( (A->inrow != B->inrow)||(A->incol != B->incol) )
   nrerror("index error 1 in mattassign");
  if( ( A->finrow != B->finrow)||(A->fincol != B->fincol) )
   nrerror("index error 2 in mattassign");
    for(i=A->inrow;i<=A->finrow;i++){
    for(j=A->incol;j<=A->fincol;j++){
      A->mat[i][j] = B->mat[i][j];
    }
    }
}

/* "make_mat" makes a square matrix from entries in a vector,
reading successively across rows;
if upper==1 assumes V stores diagonal plus upper entries,
and fills in matrix to be symmetric */

void make_mat(struct vector *V,struct matrix *M,int first,int last,int upper)
{ 
  char rt[] = "make_mat";
  int i,j,d,pl;double dum;
    if( (M->inrow != 1)||(M->incol != 1)||(M->finrow != M->fincol)||(M->type != 'm') )
     nrerror("error in make_mat:type or index error");else;
    if( (V->in > first)||(V->fin < last) )
    nrerror("error in make_mat:index error 2");else;
    d = M->finrow;
    pl = 0;
     for(i=1;i<=d;i++){
      for(j=1;j<=d;j++){
	 if( j < i){
	  if(upper==0){
	   dum = Evec(V,first + pl,rt);
	   Amat(M,i,j,dum,rt);
	   pl++;
	  }else;
	 }
	 else{
	  dum = Evec(V,first + pl,rt);
	  Amat(M,i,j,dum,rt);
	  if(upper==1) Amat(M,j,i,dum,rt); 
	  pl++;
	}
      }
     }
  if(first + pl - 1 == last); else nrerror("Error in make_mat: indices didn't agree.");
}

void Make_mat(struct vector *V,struct matrix *M,int first,int last,int upper,char st[])
{ 
  char rt[] = "make_mat";
  int i,j,d,pl;double dum;
    if( (M->inrow != 1)||(M->incol != 1)||(M->finrow != M->fincol)||(M->type != 'm') )
     nrerror("error in make_mat:type or index error");else;
    if( (V->in > first)||(V->fin < last) )
    nrerror("error in make_mat:index error 2");else;
    d = M->finrow;
    pl = 0;
     for(i=1;i<=d;i++){
      for(j=1;j<=d;j++){
	 if( j < i){
	  if(upper==0){
	   dum = Evec(V,first + pl,rt);
	   Amat(M,i,j,dum,rt);
	   pl++;
	  }else;
	 }
	 else{
	  dum = Evec(V,first + pl,rt);
	  Amat(M,i,j,dum,rt);
	  if(upper==1) Amat(M,j,i,dum,rt); 
	  pl++;
	}
      }
     }
 if(first + pl - 1 == last); 
  else werror("Error in make_mat: indices didn't agree.\n",st);
}

/* codes entries of a matrix into a vector, reading across rows;
set upper==1 if matrix is symmetric */
void mat_incode(struct vector *V,struct matrix *M,
int first,int last,int upper)
{ 
 char rt[] = "mat_incode";
   int i,j,d,pl;double dum;
    if( (M->inrow != 1)||(M->incol != 1)||(M->finrow != M->fincol)||(M->type != 'm') )
     nrerror("error in mat_incode:type or index error");else;
    if( (V->in > first)||(V->fin < last) )
    nrerror("error in mat_incode:index error 2");else;
    d = M->finrow;
     pl = 0;
     for(i=1;i<=d;i++){
      for(j=1;j<=d;j++){
	 if( j < i){
	   if(upper==0){
	    dum = Emat(M,i,j,rt);
	    Avec(V,first + pl,dum,rt);
	    pl++;
	   } else;
	 }
	 else{
	  dum = Emat(M,i,j,rt);
	  Avec(V,first + pl,dum,rt);
	  pl++;
	}
      }
     }
  if(first + pl - 1 == last); else nrerror("Error in mat_incode: indices didn't agree.");
 }

void Mat_incode(struct vector *V,struct matrix *M,
int first,int last,int upper,char st[])
{ 
 char rt[] = "mat_incode";
   int i,j,d,pl;double dum;
    if( (M->inrow != 1)||(M->incol != 1)||(M->finrow != M->fincol)||(M->type != 'm') )
     nrerror("error in mat_incode:type or index error");else;
    if( (V->in > first)||(V->fin < last) )
    nrerror("error in mat_incode:index error 2");else;
    d = M->finrow;
     pl = 0;
     for(i=1;i<=d;i++){
      for(j=1;j<=d;j++){
	 if( j < i){
	   if(upper==0){
	    dum = Emat(M,i,j,rt);
	    Avec(V,first + pl,dum,rt);
	    pl++;
	   } else;
	 }
	 else{
	  dum = Emat(M,i,j,rt);
	  Avec(V,first + pl,dum,rt);
	  pl++;
	}
      }
     }
  if(first + pl - 1 == last); else werror("Error in mat_incode: indices didn't agree.",st);
 }


void make_imat(struct ivector *V,struct imatrix *M,int first,int last,int upper)
{ 

  int i,j,d,sig,pl;
    if( (M->inrow != 1)||(M->incol != 1)||(M->finrow != M->fincol)||(M->type != 'M') )
     nrerror("error in make_mat:type or index error");else;
    if( (V->in > first)||(V->fin < last) )
    nrerror("error in make_mat:index error 2");else;
    sig = last - first + 1;d = M->finrow;
   if(upper==1){
     if(sig != (d*(d+1))/2) nrerror("error in make_imat:index error 3");else;
    }
    else{
    if(sig != d*d) nrerror("error in make_imat:index error 4");else;
    }
     pl = 0;

     for(i=1;i<=d;i++){
      for(j=1;j<=d;j++){

	 if( j < i){
	   if(upper==0){
	    M->mat[i][j] = V->vec[first + pl];
	    pl++;
	    }else;
	 }

	 else{
	   M->mat[i][j] = V->vec[first + pl];
	   if(upper==1) M->mat[j][i] = M->mat[i][j];
	   pl++;
	 }

       }
      }

      if(pl != sig) nrerror("error 4 in make_imat");else;
}

void imat_incode(struct ivector *V,struct imatrix *M,
int first,int last,int upper)
{ 
   int i,j,d,sig,pl;
    if( (M->inrow != 1)||(M->incol != 1)||(M->finrow != M->fincol)||(M->type != 'M') )
     nrerror("error in mat_incode:type or index error");else;
    if( (V->in > first)||(V->fin < last) )
    nrerror("error in :index error 2");else;
    sig = last - first + 1;d = M->finrow;
    if(upper==1){
     if(sig != d*(d+1)/2) nrerror("error in imat_incode:index error 3");else;
    }
    else{ if(sig != d*d) nrerror("error in imat_incode:index error 4");else;
    }
     pl = 0;

     for(i=1;i<=d;i++){
      for(j=1;j<=d;j++){
	 if( j < i){
	   if(upper==0){
	     V->vec[first + pl] = M->mat[i][j];
	    pl++;
	    }else;
	 }
	 else{
	   V->vec[first + pl] = M->mat[i][j];
	   pl++;
	 }
       }
      }

      if(pl != sig) nrerror("error 4 in make_mat");else;
 }


void col_extract(struct matrix *A,int n,struct vector *V)
{int j;double dum;
 char rt[] = "col_extract";

   if( A->type != 'm') nrerror("type error in A in col_extract");else;
   if( V->type != 'v') nrerror("type error in V in col_extract");else;
   if( A->inrow != V->in ) nrerror("index error in col_extract");else;
   if( A->finrow != V->fin ) nrerror("index error 2 in col_extract");else;

   for(j=A->inrow;j<=A->finrow;j++){
    dum = Emat(A,j,n,rt);
    Avec(V,j,dum,rt); 
   }
}

void col_inject(struct matrix *A,int n,struct vector *V)
{
 int j;double dum;
 char rt[] = "col_inject";

   if( A->type != 'm') nrerror("type error in A in rowinject");else;
   if( V->type != 'v') nrerror("type error in V in rowinject");else;
   if( A->inrow != V->in ) nrerror("index error in rowinject");else;
   if( A->finrow != V->fin ) nrerror("index error 2 in rowinject");else;

   for(j=A->inrow;j<=A->finrow;j++){
    dum = Evec(V,j,rt); 
    Amat(A,j,n,dum,rt);
   }
}

void row_inject(struct matrix *A,int n,struct vector *V)
{
 int j;double dum;
 char rt[] = "row_inject";

   if( A->type != 'm') nrerror("type error in A in rowinject");else;
   if( V->type != 'v') nrerror("type error in V in rowinject");else;
   if( A->incol != V->in ) nrerror("index error in rowinject");else;
   if( A->fincol != V->fin ) nrerror("index error 2 in rowinject");else;

   for(j=A->incol;j<=A->fincol;j++){
    dum = Evec(V,j,rt); 
    Amat(A,n,j,dum,rt);
   }
}

void irow_inject(struct imatrix *A,int n,struct ivector *V)
{
 int j;double dum;
 char rt[] = "irow_inject";

   if( A->type != 'M') nrerror("type error in A in rowinject");else;
   if( V->type != 'V') nrerror("type error in V in rowinject");else;
   if( A->incol != V->in ) nrerror("index error in rowinject");else;
   if( A->fincol != V->fin ) nrerror("index error 2 in rowinject");else;

   for(j=A->incol;j<=A->fincol;j++){
    dum = Eivec(V,j,rt); 
    Aimat(A,n,j,dum,rt);
   }
}

void row_extract(struct matrix *A,int n,struct vector *V)
{
 int j;double dum;
 char rt[] = "row_extract";

   if( A->type != 'm') nrerror("type error in A in rowextract");else;
   if( V->type != 'v') nrerror("type error in V in rowextract");else;
   if( A->incol != V->in ) nrerror("index error in rowextract");else;
   if( A->fincol != V->fin ) nrerror("index error 2 in rowextract");else;

   for(j=A->incol;j<=A->fincol;j++){
    dum = Emat(A,n,j,rt);
    Avec(V,j,dum,rt); 
   }
}

double dotprod(struct vector *v,struct vector *w)
{ 
 int k;double sum;
 char rt[] = "dotprod";
    if(v->type != 'v') nrerror("type error in dotprod; first vector");else;
    if(w->type != 'v') nrerror("type error in dotprod; second vector");else;
    if( (v->in != w->in)||(v->fin != w->fin) )
    nrerror("index error in dotprod");else;
    sum = 0.0;
    for(k = v->in;k <= v->fin;k++) sum = sum + Evec(v,k,rt)*Evec(w,k,rt);
    return sum;
}

double semidotprod(struct vector *v,struct vector *w,int st1,int st2,int num)
{ int k;double sum = 0.0;
    if(v->type != 'v') nrerror("type error in dotprod 1");else;
    if(w->type != 'v') nrerror("type error in dotprod 2");else;
    if( (st1+1 < v->in )||(st1 + num > v->fin) )
    nrerror("index error in semidotprod 1");else;
    if( (st2+1 < w->in )||(st2 + num > w->fin) )
    nrerror("index error in semidotprod 2");else;
    if( (st1+1 > v->fin )||(st1 + num < v->in) )
    nrerror("index error in semidotprod 3");else;
    if( (st2+1 > w->fin )||(st2 + num < w->in) )
    nrerror("index error in semidotprod 2");else;
    for(k=1;k <= num;k++) sum = sum + (v->vec[st1 + k])*(w->vec[st2 + k]);
    return sum;
}

int matvecprod(struct matrix *A,struct vector *V,struct vector *Out)
{
 int i;
 char rt[] = "matvecprod";
 
  if ( V == Out ){
   struct vector dumOut;
   invec(&dumOut,Out->in,Out->fin,"dumOut");
   for(i=A->inrow;i<=A->finrow;i++){
   int k;double sum=0.0;
   for(k=A->incol;k<=A->fincol;k++){
    sum = sum + Evec(V,k,rt)*Emat(A,i,k,rt);
    }
    Avec(&dumOut,i,sum,rt);
   }
   vecassign(Out,&dumOut);
   freevec(&dumOut);
  }

  else{
   for(i=A->inrow;i<=A->finrow;i++){
   int k;double sum=0.0;
   for(k=A->incol;k<=A->fincol;k++){
    sum = sum + Evec(V,k,rt)*Emat(A,i,k,rt);
    }
    Avec(Out,i,sum,rt);
   }
  }
   return 1;
 }

int vecmatprod(struct vector *V,struct matrix *A,struct vector *Out)
{
 int i;
 char rt[] = "vecmatprod";
 
  if ( V == Out ){
   struct vector dumOut;
   invec(&dumOut,Out->in,Out->fin,"dumOut");
  
   for(i=A->incol;i<=A->fincol;i++){
   int k;double sum=0.0;
   for(k=A->inrow;k<=A->finrow;k++){
    sum = sum + Evec(V,k,rt)*Emat(A,k,i,rt);
    }
    Avec(&dumOut,i,sum,rt);
   }
   vecassign(Out,&dumOut);
   freevec(&dumOut);
  }

  else{
    for(i=A->incol;i<=A->fincol;i++){
     int k;double sum=0.0;
      for(k=A->inrow;k<=A->finrow;k++){
      sum = sum + Evec(V,k,rt)*Emat(A,k,i,rt);
      }
     Avec(Out,i,sum,rt);
    }
   }
   return 1;
 }


int veccombo(double a,struct vector *U,double b,struct vector *V,
struct vector *out)
{int i;
  if( (U->type != 'v')||(V->type != 'v')||(out->type != 'v') )
  nrerror("Type error in veccombo");
  if( (U->in != V->in)||(U->fin != V->fin) )
   nrerror("Index error 1 in veccombo");
   if( (U->in != out->in)||(U->fin != out->fin) )
   nrerror("Index error 2 in veccombo");

  for(i=U->in;i<=U->fin;i++){
      if( (a == 1.0)&&(b == 1.0) ) out->vec[i] = U->vec[i] + V->vec[i];
      else{
       if( a == 0.0 ) out->vec[i] = b*V->vec[i];
       else{
        if( b == 0.0 )  out->vec[i] = a*U->vec[i];
        else  out->vec[i] = a*U->vec[i] + b*V->vec[i];
       }
      }
  }
  return 1;
}


int matcombo(double u,struct matrix *A,double v,struct matrix *B,
struct matrix *out)
{int i,j;
  if( (A->type != 'm')||(B->type != 'm')||(out->type != 'm') )
  nrerror("Type error in matcombo");
  if( (A->inrow != B->inrow)||(A->finrow != B->finrow) )
   nrerror("Index error 1 in matcombo");
  if( (A->incol != B->incol)||(A->fincol != B->fincol) )
   nrerror("Index error 2 in matcombo");
   if( (B->inrow != out->inrow)||(B->finrow != out->finrow) )
   nrerror("Index error 3 in matcombo");
  if( (B->incol != out->incol)||(B->fincol != out->fincol) )
   nrerror("Index error 4 in matcombo");

  for(i=A->inrow;i<=A->finrow;i++){
  for(j=A->incol;j<=A->fincol;j++){
      if( ( u == 1.0)&&(v == 1.0) ) out->mat[i][j] = A->mat[i][j] + B->mat[i][j];
      else out->mat[i][j] = u*A->mat[i][j] + v*B->mat[i][j];
  }
  }
  return 1;
}


int matprod(int u,struct matrix *A,int v,struct matrix *B,
struct matrix *out) /* u=1 or v=1 means take tranpose of A or B (resp.) */
{

 int i,j,k;char rt[] = "matprod";

 if( (A->type != 'm')||(B->type != 'm')||(out->type != 'm') )
 nrerror("Error in matprod 1");

if( (u==0)&&(v==0) ){
  if( (A->incol != B->inrow)||(A->fincol != B->finrow) )
   nrerror("Error in matprod 2");
  if( (A->inrow != out->inrow)||(A->finrow != out->finrow) )
   nrerror("Error in matprod 3");
  if( (B->incol != out->incol)||(B->fincol != out->fincol) )
   nrerror("Error in matprod 4");

  if( (A == out)||(B == out) ){
   struct matrix dumout;
   inmat(&dumout,out->inrow,out->finrow,out->incol,out->fincol,"dumout");
  for(i=A->inrow;i<=A->finrow;i++){
  for(j=B->incol;j<=B->fincol;j++){
  double sum=0.0;
   for(k=A->incol;k<=A->fincol;k++){
    /* sum = sum + A->mat[i][k]*B->mat[k][j]; */
    sum = sum + Emat(A,i,k,rt)*Emat(B,k,j,rt);
   }
    /* dumout.mat[i][j] = sum; */
    Amat(&dumout,i,j,sum,rt);
  }
  }

  /* for(i=A->inrow;i<=A->finrow;i++){
  for(j=B->incol;j<=B->fincol;j++){
    out->mat[i][j] = dumout.mat[i][j];
  }
  } */
  matassign(out,&dumout);
  freemat(&dumout);
 }

 else{

  for(i=A->inrow;i<=A->finrow;i++){
  for(j=B->incol;j<=B->fincol;j++){
  int k;double sum=0.0;
  for(k=A->incol;k<=A->fincol;k++){
    /* sum = sum + (A->mat[i][k])*(B->mat[k][j]); */
    sum = sum + Emat(A,i,k,rt)*Emat(B,k,j,rt);
   }
  /*out->mat[i][j] = sum;*/
  Amat(out,i,j,sum,rt);
  }
  }

 }

}else;

if( (u==0)&&(v==1) ){

  if( (A->incol != B->incol)||(A->fincol != B->fincol) )
   nrerror("Error in matprod 2");
  if( (A->inrow != out->inrow)||(A->finrow != out->finrow) )
   nrerror("Error in matprod 3");
  if( (B->inrow != out->incol)||(B->finrow != out->fincol) )
   nrerror("Error in matprod 4");

  if( (A == out)||(B == out) ){
   struct matrix dumout1;
   inmat(&dumout1,out->inrow,out->finrow,out->incol,out->fincol,"dumout1");
  for(i=A->inrow;i<=A->finrow;i++){
  for(j=B->inrow;j<=B->finrow;j++){
  double sum=0.0;
  for(k=A->incol;k<=A->fincol;k++){
    /* sum = sum + A->mat[i][k]*B->mat[j][k]; */
    sum = sum + Emat(A,i,k,rt)*Emat(B,j,k,rt);
   }
  /* dumout1.mat[i][j] = sum; */
  Amat(&dumout1,i,j,sum,rt);
  }
  }

  /* for(i=out->inrow;i<=out->finrow;i++){
  for(j=out->incol;j<=out->fincol;j++){
    out->mat[i][j] = dumout1.mat[i][j];
  }
  } */
  matassign(out,&dumout1);
  freemat(&dumout1);
 }

 else{

  for(i=A->inrow;i<=A->finrow;i++){
  for(j=B->inrow;j<=B->finrow;j++){
  int k;double sum=0.0;
  for(k=A->incol;k<=A->fincol;k++){
    /* sum = sum + A->mat[i][k]*B->mat[j][k]; */
    sum = sum + Emat(A,i,k,rt)*Emat(B,j,k,rt);
   }
  /* out->mat[i][j] = sum; */
  Amat(out,i,j,sum,rt);
  }
  }

 }

} else;

if( (u==1)&&(v==0) ){
     if( (A->inrow != B->inrow)||(A->finrow != B->finrow) )
   nrerror("Error in matprod 2");
  if( (A->incol != out->inrow)||(A->fincol != out->finrow) )
   nrerror("Error in matprod 3");
  if( (B->incol != out->incol)||(B->fincol != out->fincol) )
   nrerror("Error in matprod 4");

  if( (A == out)||(B == out) ){
   struct matrix dumout2;
   inmat(&dumout2,out->inrow,out->finrow,out->incol,out->fincol,"dumout2");
  for(i=A->incol;i<=A->fincol;i++){
  for(j=B->incol;j<=B->fincol;j++){
  double sum=0.0;
   for(k=A->inrow;k<=A->finrow;k++){
   /* sum = sum + A->mat[k][i]*B->mat[k][j]; */
    sum = sum + Emat(A,k,i,rt)*Emat(B,k,j,rt);
   }
    /* dumout2.mat[i][j] = sum; */
    Amat(&dumout2,i,j,sum,rt);
  }
  }

  /* for(i=out->inrow;i<=out->finrow;i++){
  for(j=out->incol;j<=out->fincol;j++){
    out->mat[i][j] = dumout2.mat[i][j];
  }
  } */
  matassign(out,&dumout2);
  freemat(&dumout2);
 }

 else{

  for(i=A->incol;i<=A->fincol;i++){
  for(j=B->incol;j<=B->fincol;j++){
  int k;double sum=0.0;
  for(k=A->inrow;k<=A->finrow;k++){
    /* sum = sum + (A->mat[k][i])*(B->mat[k][j]); */
    sum = sum + Emat(A,k,i,rt)*Emat(B,k,j,rt);
   }
   /* out->mat[i][j] = sum; */
   Amat(out,i,j,sum,rt);
  }
  }

 }

}else;

 if( u*v > 0) nrerror("Error in matprod: both are transposed??");else;
 return 1;

}



/* ludcmp: derived from NUMC, p. 43; replaces a square matrix
"a" by an LU decomposition; indx is an ivector recording
permutation used for pivoting; "d" = +1.0 or -1.0 records parity of
interchanges. NOTE: THIS ROUTINE DESTROYS THE ORIGINAL MATRIX! */ 

int ludcmp(struct matrix *A,struct ivector *indx,double *d)
{ int n,i,imax,j,k,halt;
  double big,dum,sum,temp,TINY;
  struct vector vv;

   if( A->type != 'm') nrerror("type error in ludcmp");
   if( (A->inrow != 1)||(A->incol != 1) ) nrerror("index error in ludcmp");
   if( A->finrow != A->fincol ) nrerror("index error 2 in ludcmp");
    n = A->finrow; 
    TINY = 1.0e-20;
   if( (indx->in != 1)||(indx->fin != n )) nrerror("index error 4 in ludcmp");else;
   invec(&vv,1,n,"vv");
   *d =1.0;halt=0;
   for(i=1;i<=n;i++){
     big=0.0;
     for(j=1;j<=n;j++) if( (temp=fabs(A->mat[i][j])) > big) big=temp;else;
     if(big > 0.0) 
     vv.vec[i] = 1.0/big;
     else{
       halt++; /* nrerror("Singular matrix in ludcmp!"); */
     }
   }

 if(halt==0){

   for(j=1;j<=n;j++){
       for(i=1;i<j;i++){
         sum = A->mat[i][j]; 
	  for(k=1;k<i;k++) sum -= A->mat[i][k]*A->mat[k][j];
	 A->mat[i][j] = sum;
	}
	big=0.0;
	for(i=j;i<=n;i++){
	  sum = A->mat[i][j];
	  for(k=1;k<j;k++) sum -= A->mat[i][k]*A->mat[k][j];
	  A->mat[i][j] = sum;
	  if( (dum=vv.vec[i]*fabs(sum)) >= big){
	   big=dum;
	   imax=i;
	  } else;
        }
	if(j != imax){
	  for(k=1;k<=n;k++){
	   dum = A->mat[imax][k];
	   A->mat[imax][k] = A->mat[j][k];
	   A->mat[j][k] = dum;
          }
	  *d = -(*d);
	  vv.vec[imax] = vv.vec[j];
	} else;
	indx->vec[j] = imax;
        if(A->mat[j][j] == 0.0) A->mat[j][j] = TINY;
	if(j != n){
	  dum = 1.0/(A->mat[j][j]);
	  for(i=j+1;i<=n;i++) A->mat[i][j] *= dum;
	} else;
     }


  } else; /* end of halt conditional */

     freevec(&vv);
     return halt;

 }

/* lubksb: derived from NUMC, p. 43; solves Ax=b
given an LU decomp. of A. Returns x in vector "b".
Run ludcmp first. To solve many equations with same "A" but different "b",
run ludcmp once, then lubskb as desired.*/

void lubskb(struct matrix *A,struct ivector *indx,struct vector *b)
{int n,i,ii=0,ip,j;double sum;
  n = A->finrow;
  if( (b->in != 1)||(b->fin != n ) ) nrerror("index error in lubskb");else;
 for(i=1;i<=n;i++){
    ip = indx->vec[i];
    sum = b->vec[ip];
    b->vec[ip] = b->vec[i];
    if(ii) for(j=ii;j<=i-1;j++) sum -= A->mat[i][j]*b->vec[j];
    else if(sum) ii = i;
    b->vec[i] = sum;
 }
 for(i=n;i>=1;i--){
  sum = b->vec[i];
   for(j=i+1;j<=n;j++) sum -= A->mat[i][j]*b->vec[j];
   b->vec[i] = sum/(A->mat[i][i]);
   }
}


/* invert: derived from NUMC, p. 45; computes an inverse matrix.
 Computes det(A) and returns it; halts if it is 0, otherwise
 replaces A by inverse of A.
 NOTE: THIS ROUTINE DESTROYS THE ORIGINAL MATRIX A! */

double invert(struct matrix *A,char rt[])
{ int n,i,j;
  struct matrix Ainv;
  struct ivector indx;struct vector col;double d;
  if(A->type != 'm') nrerror("type error in invert");
  if( (A->inrow != 1)||(A->incol != 1) ) nrerror("index error 1 in invert");
  if( A->finrow != A->fincol)  nrerror("index error 2 in invert");
   n = A->finrow;
   inivec(&indx,1,n,"indx");invec(&col,1,n,"col");inmat(&Ainv,1,n,1,n,"Ainv");
  if( ludcmp(A,&indx,&d) == 0); 
  else{
    printf("Error in invert in routine %s: singular matrix: ...\n",rt);
    printmat(A);
    printf("\n");
    nrerror("error in invert: singular matrix");
  }
   for(j=1;j<=n;j++) d *= A->mat[j][j];
   if( d != 0.0 ){
    for(j=1;j<=n;j++){
     for(i=1;i<=n;i++) col.vec[i] = 0.0;
     col.vec[j] = 1.0;
     lubskb(A,&indx,&col);
     for(i=1;i<=n;i++) Ainv.mat[i][j] = col.vec[i];
    }
    matassign(A,&Ainv);
   } else;
   freeivec(&indx);freevec(&col);freemat(&Ainv);
   return d;
}

/* invert_alt below doesn't halt the program if the 
matrix is singular, but rather returns 0.0. */
 
double invert_alt(struct matrix *A)
{ int n,i,j;
  struct matrix Ainv;
  struct ivector indx;struct vector col;double d;
  if(A->type != 'm') nrerror("type error in invert");
  if( (A->inrow != 1)||(A->incol != 1) ) nrerror("index error 1 in invert");
  if( A->finrow != A->fincol)  nrerror("index error 2 in invert");
   n = A->finrow;
   inivec(&indx,1,n,"indx");invec(&col,1,n,"col");inmat(&Ainv,1,n,1,n,"Ainv");
  if( ludcmp(A,&indx,&d) == 0){ 
   for(j=1;j<=n;j++) d *= A->mat[j][j];
   if( d != 0.0 ){
    for(j=1;j<=n;j++){
     for(i=1;i<=n;i++) col.vec[i] = 0.0;
     col.vec[j] = 1.0;
     lubskb(A,&indx,&col);
     for(i=1;i<=n;i++) Ainv.mat[i][j] = col.vec[i];
    }
    matassign(A,&Ainv);
   } else;
  }
  else d = 0.0; /*nrerror("error in invert: singular matrix");*/
   freeivec(&indx);freevec(&col);freemat(&Ainv);
   return d;
}

double checkinverse(struct matrix *A,struct matrix *Ainv)
{
  struct matrix out; double error; int i,j,n;

  if( (A->type != 'm')||(Ainv->type != 'm') )
  nrerror("error 1 in checkinverse");
  if( (A->inrow != Ainv->inrow)||(A->incol != Ainv->incol) )
   nrerror("error 2 in checkinverse");
  if( (A->finrow != Ainv->finrow)||(A->fincol != Ainv->fincol) )
   nrerror("error 3 in checkinverse");
   n = A->finrow;
   if(n<1) nrerror("error 4 in checkinverse");
  inmat(&out,1,n,1,n,"out");

  matprod(0,A,0,Ainv,&out);

  error = 0.0;

  for(i=1;i<=n;i++){
  for(j=1;j<=n;j++){
  if(i==j) error = error + fabs(1.0 - out.mat[i][j]);
  else error = error + fabs(out.mat[i][j]);
  }
  }
   error = error;

  freemat(&out);
  return error;
}

/* det: computes derminant of a matrix. Can also use
"invert" which returns det. NOTE: THIS ROUTINE
DESTROYS THE MATRIX A!, replacing it by an LU decomposition. */

double det(struct matrix *A)
{ int n,j,halt;struct ivector indx;double d;
  if(A->type != 'm')  nrerror("type error 1 in invert");
  if( A->inrow != 1) nrerror("type error 2 in invert");
   n = A->finrow;
   inivec(&indx,1,n,"indx");
   if( (halt = ludcmp(A,&indx,&d)) == 0){
   for(j=1;j<=n;j++) d *= A->mat[j][j];
   } 
   else d = 0.0;
   freeivec(&indx);
   return d;
}

void transpose(struct matrix *A,struct matrix *A_tr)
{ 
  int j,k;double dum;
  char rt[] = "transpose";

  if( (A->type != 'm')||(A_tr->type != 'm') ) nrerror("Type error in transpose");else;
  if( (A->inrow != A_tr->incol)||(A->finrow != A_tr->fincol) ) nrerror("Index error 1 in transpose");else;
  if( (A->incol != A_tr->inrow)||(A->fincol != A_tr->finrow) ) nrerror("Index error 2 in transpose");else;

   if(A == A_tr){
   for(j=A->inrow;j<=A->finrow;j++){
    for(k=A->incol;k<j;k++){
    matswap(A,j,k,A,k,j); 
    }
   }
  }

 else{
   for(j=A->inrow;j<=A->finrow;j++){
    for(k=A->incol;k<=A->fincol;k++){
    dum = Emat(A,j,k,rt);
    Amat(A_tr,k,j,dum,rt);
    }
   }
  }

}


int diag(struct vector *V,struct matrix *A)
{int j,k;

   if( A->type != 'm') nrerror("error 1 in diag");else;
   if( V->type != 'v') nrerror("error 2 in diag");else;
   if( A->inrow != A->incol ) nrerror("error 3 in diag");else;
   if( A->finrow != A->fincol ) nrerror("error 4 in diag");else;
   if( A->inrow != V->in ) nrerror("error 5 in diag");else;
   if( A->finrow != V->fin ) nrerror("error 6 in diag");else;

   for(j=A->inrow;j<=A->finrow;j++){
    for(k=A->inrow;k<=A->finrow;k++){
      if( j == k) A->mat[j][j] = V->vec[j];
      else A->mat[j][k] = 0.0;
    }
   }

   return 1;
}

int printvec(int one_line,int print_name,struct vector *V)
{int q;
  char rt[] = "printvec";
 
   if( V->type != 'v' ) nrerror("Type error in printvec.");else;
   if(one_line == 1){
     if(print_name == 1) for(q=V->in;q<=V->fin;q++) printf("vector[%i] = %g ",q,Evec(V,q,rt));
     else for(q=V->in;q<=V->fin;q++) printf("%g ",Evec(V,q,rt));
     printf("\n");
   }
   else{
    if(print_name == 1) for(q=V->in;q<=V->fin;q++) printf("vector[%i] = %f\n",q,Evec(V,q,rt));
    else for(q=V->in;q<=V->fin;q++) printf("%f\n",Evec(V,q,rt));
    printf("\n");
   }
   return 1;
}

int printivec(int one_line,int print_name,struct ivector *V)
{int q;
  char rt[] = "printivec";
 
   if( V->type != 'V' ) nrerror("Type error in printvec.");else;
   if(one_line == 1){
     if(print_name == 1) for(q=V->in;q<=V->fin;q++) printf("vector[%i] = %i ",q,Eivec(V,q,rt));
     else for(q=V->in;q<=V->fin;q++) printf("%i ",Eivec(V,q,rt));
     printf("\n"); 
   }
   else{
    if(print_name == 1) for(q=V->in;q<=V->fin;q++) printf("vector[%i] = %i\n",q,Eivec(V,q,rt));
    else for(q=V->in;q<=V->fin;q++) printf("%i\n",Eivec(V,q,rt));
     printf("\n"); 
   }
   return 1;
}

int fileprintvec(FILE *ptr,int one_line,int print_name,struct vector *V)
{int q;
  char rt[] = "fileprintvec";
   
   if( V->type != 'v' ) nrerror("Type error in fileprintvec.");else;
   if(one_line == 1){
     if(print_name == 1) for(q=V->in;q<=V->fin;q++) fprintf(ptr,"vector[%i] = %g ",q,Evec(V,q,rt));
     else for(q=V->in;q<=V->fin;q++) fprintf(ptr,"%g ",Evec(V,q,rt));
     fprintf(ptr,"\n\n");
   }
   else{
    if(print_name == 1) for(q=V->in;q<=V->fin;q++) fprintf(ptr,"vector[%i] = %g\n",q,Evec(V,q,rt));
    else for(q=V->in;q<=V->fin;q++) fprintf(ptr,"%g\n",Evec(V,q,rt));
    fprintf(ptr,"\n");
   }
   return 1;
}

int fileprintivec(FILE *ptr,int one_line,int print_name,struct ivector *V)
{int q;
  char rt[] = "fileprintivec";
   
   if( V->type != 'V' ) nrerror("Type error in fileprintivec.");else;
   if(one_line == 1){
     if(print_name == 1) for(q=V->in;q<=V->fin;q++) fprintf(ptr,"vector[%i] = %i ",q,Eivec(V,q,rt));
     else for(q=V->in;q<=V->fin;q++) fprintf(ptr,"%i ",Eivec(V,q,rt));
     fprintf(ptr,"\n\n");
   }
   else{
    if(print_name == 1) for(q=V->in;q<=V->fin;q++) fprintf(ptr,"vector[%i] = %i\n",q,Eivec(V,q,rt));
    else for(q=V->in;q<=V->fin;q++) fprintf(ptr,"%i\n",Eivec(V,q,rt));
    fprintf(ptr,"\n");
   }
   return 1;
}

int printrow(struct matrix *M,int row)
{int p;
  char rt[] = "printrow";
 if( M->type != 'm' ) nrerror("type error in printrow");else;
 if( (row < M->inrow)||(row > M->finrow) ) nrerror("error in printrow: row outside index range");else;
  for(p=M->incol;p<=M->fincol;p++){
     printf("%f  ",Emat(M,row,p,rt));
    }
   printf("\n");
    return 1;
}

int printmat_a(struct matrix *M,int first,int last)
{int q,p;
  char rt[] = "printmat_a";
 if( M->type != 'm' ) nrerror("type error in printmat");else;
 if( (first < M->inrow)||(last > M->finrow) ) nrerror("Index error in printmat_a");else;
  for(q=first;q<=last;q++){
   for(p=M->incol;p<=M->fincol;p++){
     printf("%f  ",Emat(M,q,p,rt));
    }
   printf("\n");
  }
   printf("\n");
    return 1;
}

int printmat(struct matrix *M)
{int q,p;
  char rt[] = "printmat";
 if( M->type != 'm' ) nrerror("type error in printmat");else;
for(q=M->inrow;q<=M->finrow;q++){
 for(p=M->incol;p<=M->fincol;p++){
   printf("%f  ",Emat(M,q,p,rt));
  }
  printf("\n");
}
 printf("\n");
 return 1;
}

int printimat(struct imatrix *M)
{int q,p;
  char rt[] = "printimat";
 if( M->type != 'M' ) nrerror("type error in printmat");else;
for(q=M->inrow;q<=M->finrow;q++){
 for(p=M->incol;p<=M->fincol;p++){
   printf("%i  ",Eimat(M,q,p,rt));
  }
  printf("\n");
}
 printf("\n");
 return 1;
}

int fileprintmat(FILE *out,int printtype,struct matrix *M)
{
 int q,p;
 char rt[] = "fileprintmat";

 if( M->type != 'm' ) nrerror("type error in printmat");else;
 for(q=M->inrow;q<=M->finrow;q++){
  for(p=M->incol;p<=M->fincol;p++){
   if(printtype == 1) fprintf(out,"%f  ",Emat(M,q,p,rt));
   else fprintf(out,"%g  ",Emat(M,q,p,rt));
  }
  fprintf(out,"\n");
 }
 fprintf(out,"\n");
 return 1;
}

int fileprintmat_a(FILE *ptr,struct matrix *M,int first,int last)
{int q,p;
 char rt[] = "fileprintmat_a";
 if( M->type != 'm' ) nrerror("type error in printmat");else;
 if( (first < M->inrow)||(last > M->finrow) ) nrerror("Index error in printmat_a");else;
  for(q=first;q<=last;q++){
   for(p=M->incol;p<=M->fincol;p++){
    fprintf(ptr,"%f  ",Emat(M,q,p,rt));
    }
   fprintf(ptr,"\n");
  }
   fprintf(ptr,"\n");
    return 1;
}

