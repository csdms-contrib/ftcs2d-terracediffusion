#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

#define PI 3.141592653589793

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
      int  i,**m;

       /*allocate pointers to rows */
        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
      m -= nrl;

       /*allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
      }
       /* return pointer to array of pointers to rows */
        return m;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

main()
{    FILE *fp1,*fp2;
     int length_x,length_y,delta,lattice_size_x,lattice_size_y,**rills,**efficiencymask,i,j;
     float deltah,simulationlength,elapsedtime,timestep,dissectionrelief,kappa,depositslope;
     float depositaspect,northwestdatum,**topo,**topoold,**initialtopo; 

     fp1=fopen("rillmask","r");
     fp2=fopen("terracediffusion2d","w"); 
     length_x=300;             /* m */
     length_y=300;
     delta=1.0;                /* m */
     dissectionrelief=10.0;
     kappa=1.0;                /* m^2/kyr */
     depositslope=0.06;        /* m/m */
     depositaspect=300*PI/180; /* surface drains SE; degrees converted to radians */
     northwestdatum=0.0;       /* elevations are relative to NE corner = 0 m */
     lattice_size_x=length_x/delta;
     lattice_size_y=length_y/delta;
     topo=matrix(1,lattice_size_x,1,lattice_size_y);
     topoold=matrix(1,lattice_size_x,1,lattice_size_y);
     initialtopo=matrix(1,lattice_size_x,1,lattice_size_y);
     rills=imatrix(1,lattice_size_x,1,lattice_size_y);
     efficiencymask=imatrix(1,lattice_size_x,1,lattice_size_y);
     timestep=0.5*delta*delta/(2*kappa);
     simulationlength=50;      /* kyr */ 
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {fscanf(fp1,"%d",&rills[i][j]);
        if (rills[i][j]>0) rills[i][j]=0; else rills[i][j]=1; 
        efficiencymask[i][j]=1;
        if (i<lattice_size_x) if (rills[i+1][j]==1) efficiencymask[i][j]=0;
        if (i>1) if (rills[i-1][j]==1) efficiencymask[i][j]=0; 
        if (j<lattice_size_y) if (rills[i][j+1]==1) efficiencymask[i][j]=0; 
        if (j>1) if (rills[i][j-1]==1) efficiencymask[i][j]=0;
        topo[i][j]=northwestdatum-(i-1)*delta*depositslope*cos(depositaspect)+
         (j-1)*delta*depositslope*sin(depositaspect);
        if (rills[i][j]==1) topo[i][j]-=dissectionrelief;
        topoold[i][j]=topo[i][j];
        initialtopo[i][j]=topo[i][j];}
     elapsedtime=0.0;
     while (elapsedtime<=simulationlength)
      {printf("%f %f\n",elapsedtime,timestep);
       for (j=2;j<=lattice_size_y-1;j++)
        for (i=2;i<=lattice_size_x-1;i++)    
         if (efficiencymask[i][j]==0)
          {if (rills[i][j]==0)
            deltah=timestep*kappa/(delta*delta)*(topoold[i+1][j]+topoold[i-1][j]+
             topoold[i][j+1]+topoold[i][j-1]-4*topoold[i][j]);
           else
            deltah=0; 
           topo[i][j]+=deltah;
           if (fabs(topo[i][j]-initialtopo[i][j])>0.01*dissectionrelief)
            {efficiencymask[i+1][j]=0;
             efficiencymask[i-1][j]=0;
             efficiencymask[i][j+1]=0;
             efficiencymask[i][j-1]=0;}}
       elapsedtime+=timestep;
       i=1;
       for (j=1;j<=lattice_size_y;j++) 
        topo[i][j]=topo[i+1][j]+depositslope*cos(depositaspect)*delta;
       i=lattice_size_x;
       for (j=1;j<=lattice_size_y;j++)
        topo[i][j]=topo[i-1][j]-depositslope*cos(depositaspect)*delta;        
       j=1;
       for (i=1;i<=lattice_size_x;i++)
        topo[i][j]=topo[i][j+1]-depositslope*sin(depositaspect)*delta;    
       j=lattice_size_y;
       for (i=1;i<=lattice_size_x;i++)
        topo[i][j]=topo[i][j-1]+depositslope*sin(depositaspect)*delta;           
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         topoold[i][j]=topo[i][j];}
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       fprintf(fp2,"%f\n",topo[i][j]);
     fclose(fp1);
     fclose(fp2);
}  
