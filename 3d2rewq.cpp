#include <stdio.h>
//#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#define PIE 3.1415926

int main(int argc, char **argv)
{
    int i,j,k,kk,kkk,l;
    int nx,ny,nz,lt,nedge;
    int nleft,nright,nfront,nback,ntop,nbottom;
    float frequency;
    float velmax;
    float dt;
    int ncx_shot1,ncy_shot1,ncz_shot;
    int ishot,ncy_shot,ncx_shot;
    float unit;
    int nxshot,nyshot,dxshot,dyshot;
    FILE  *fin, *fout, *flog;
    struct timeval start,end;
    float all_time;

    float *u, *v, *w, *up, *up1, *up2, 
          *vp, *vp1, *vp2, *wp, *wp1, *wp2, 
          *us, *us1, *us2, *vs, *vs1, *vs2,
          *ws, *ws1, *ws2, *initial;
    float c[5][7];
    float *wave;
    float nshot,t0,tt,c0;
    float dtx,dtz,dr1,dr2;
    float xmax;
    float vvp2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
          tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;
    if(argc<4)
    {
        printf("please add 3 parameter: inpurfile, outfile, logfile\n");
        exit(0);
    }

    fin = fopen(argv[1], "r");
    if(fin == NULL)
    {
        printf("file %s is  not exist\n", argv[1]);
        exit(0);
    }
    fscanf(fin,"nx=%d\n",&nx);
    fscanf(fin,"ny=%d\n",&ny);
    fscanf(fin,"nz=%d\n",&nz);
    fscanf(fin,"lt=%d\n",&lt);
    fscanf(fin,"nedge=%d\n",&nedge);
    fscanf(fin,"ncx_shot1=%d\n",&ncx_shot1);
    fscanf(fin,"ncy_shot1=%d\n",&ncy_shot1);
    fscanf(fin,"ncz_shot=%d\n",&ncz_shot);
    fscanf(fin,"nxshot=%d\n",&nxshot);
    fscanf(fin,"nyshot=%d\n",&nyshot);
    fscanf(fin,"frequency=%f\n",&frequency);
    fscanf(fin,"velmax=%f\n",&velmax);
    fscanf(fin,"dt=%f\n",&dt);
    fscanf(fin,"unit=%f\n",&unit);
    fscanf(fin,"dxshot=%d\n",&dxshot);
    fscanf(fin,"dyshot=%d\n",&dyshot);
    fclose(fin);

    u       = (float*)malloc(sizeof(float)*nz*ny*nx);
    v       = (float*)malloc(sizeof(float)*nz*ny*nx);
    w       = (float*)malloc(sizeof(float)*nz*ny*nx);
    up      = (float*)malloc(sizeof(float)*nz*ny*nx);
    up1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    up2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vp      = (float*)malloc(sizeof(float)*nz*ny*nx);
    vp1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vp2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    wp      = (float*)malloc(sizeof(float)*nz*ny*nx);
    wp1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    wp2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    us      = (float*)malloc(sizeof(float)*nz*ny*nx);
    us1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    us2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vs      = (float*)malloc(sizeof(float)*nz*ny*nx);
    vs1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    vs2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    ws      = (float*)malloc(sizeof(float)*nz*ny*nx);
    ws1     = (float*)malloc(sizeof(float)*nz*ny*nx);
    ws2     = (float*)malloc(sizeof(float)*nz*ny*nx);
    initial = (float*)malloc(sizeof(float)*nz*ny*nx);
    wave = (float*)malloc(sizeof(float)*lt);

    nshot=nxshot*nyshot;
    t0=1.0/frequency;

    for(l=0;l<lt;l++)
    {
        tt=l*dt;
        tt=tt-t0;
        float sp=PIE*frequency*tt;
        float fx=100000.*exp(-sp*sp)*(1.-2.*sp*sp);
        wave[l]=fx;
    }

    c0=-2.927222164;
    c[0][0]=1.66666665;
    c[1][0]=-0.23809525;
    c[2][0]=0.03968254;
    c[3][0]=-0.004960318;
    c[4][0]=0.0003174603;

    c[0][1]=0.83333;
    c[1][1]=-0.2381;
    c[2][1]=0.0595;
    c[3][1]=-0.0099;
    c[4][1]=0.0008;

    for(i=0;i<5;i++)
        for(j=0;j<5;j++)
            c[j][2+i]=c[i][1]*c[j][1];

    dtx=dt/unit;
    dtz=dt/unit;

    dr1=dtx*dtx/2.;
    dr2=dtz*dtz/2.;

    fout=fopen(argv[2],"wb");

    for(i = 0; i < nz; i++)
        for(j = 0; j < ny; j++)
            for(k = 0; k < nx; k++)
                initial[i * ny * nx + j * nx + k] = 100.0f;

    int size = sizeof(float) * nx * ny * nz;

    for(ishot=1;ishot<=nshot;ishot++)
    {
        
        printf("shot=%d\n",ishot);
        ncy_shot=ncy_shot1+(ishot/nxshot)*dyshot;
        ncx_shot=ncx_shot1+(ishot%nxshot)*dxshot;
/*

此处是初始化参数，以前是0，后来改为100，也可改为其他或者随机赋值 ,up up1 up2可不修改
*/
        memcpy(u, initial, size);
        memcpy(v, initial, size);
        memcpy(w, initial, size);
        memset(up, 0, size);
        memset(up1, 0, size);
        memset(up2, 0, size);
        memcpy(vp, initial, size);
        memcpy(vp1, initial, size);
        memcpy(vp2, initial, size);
        memcpy(wp, initial, size);
        memcpy(wp1, initial, size);
        memcpy(wp2, initial, size);
        memcpy(us, initial, size);
        memcpy(us1, initial, size);
        memcpy(us2, initial, size);
        memcpy(vs, initial, size);
        memcpy(vs1, initial, size);
        memcpy(vs2, initial, size);
        memcpy(ws, initial, size);
        memcpy(ws1, initial, size);
        memcpy(ws2, initial, size);

        /*for(i=0;i<nz;i++)
            for(j=0;j<ny;j++)
                for(k=0;k<nx;k++)
                {
                    u[i*ny*nx+j*nx+k]=100.0f;
                    v[i*ny*nx+j*nx+k]=100.0f;
                    w[i*ny*nx+j*nx+k]=100.0f;
                    up[i*ny*nx+j*nx+k]=0.0f;
                    up1[i*ny*nx+j*nx+k]=0.0f;
                    up2[i*ny*nx+j*nx+k]=0.0f;
                    vp[i*ny*nx+j*nx+k]=100.0f;
                    vp1[i*ny*nx+j*nx+k]=100.0f;
                    vp2[i*ny*nx+j*nx+k]=100.0f;
                    wp[i*ny*nx+j*nx+k]=100.0f;
                    wp1[i*ny*nx+j*nx+k]=100.0f;
                    wp2[i*ny*nx+j*nx+k]=100.0f;
                    us[i*ny*nx+j*nx+k]=100.0f;
                    us1[i*ny*nx+j*nx+k]=100.0f;
                    us2[i*ny*nx+j*nx+k]=100.0f;
                    vs[i*ny*nx+j*nx+k]=100.0f;
                    vs1[i*ny*nx+j*nx+k]=100.0f;
                    vs2[i*ny*nx+j*nx+k]=100.0f;
                    ws[i*ny*nx+j*nx+k]=100.0f;
                    ws1[i*ny*nx+j*nx+k]=100.0f;
                    ws2[i*ny*nx+j*nx+k]=100.0f;
                }//for(k=0;k<nx;k++) end*/
        for(l=1;l<=lt;l++)
        {
           
            xmax=l*dt*velmax;
            nleft=ncx_shot-xmax/unit-10;
            nright=ncx_shot+xmax/unit+10;
            nfront=ncy_shot-xmax/unit-10;
            nback=ncy_shot+xmax/unit+10;
            ntop=ncz_shot-xmax/unit-10;
            nbottom=ncz_shot+xmax/unit+10;
            
            if(nleft<5) nleft=5;
            if(nright>nx-5) nright=nx-5;
            if(nfront<5) nfront=5;
           if(nback>ny-5) nback=ny-5;
            if(ntop<5) ntop=5;
           if(nbottom>nz-5) nbottom=nz-5;

            ntop = ntop-1;
            nfront = nfront-1;
            nleft = nleft-1;
          
            for(k=ntop;k<nbottom;k++)
                for(j=nfront;j<nback;j++)
                    for(i=nleft;i<nright;i++)
                    {
                        if(k < 210) {
                            vvp2 = 2300 * 2300;
                            vvs2 = 1232 * 1232;
                        } else if(k < 260) {
                            vvp2 = 2800 * 2800;
                            vvs2 = 1509 * 1509;
                        } else {
                            vvp2 = 3500 * 3500;
                            vvs2 = 1909 * 1909;
                        }

                        tempux2=0.0f;
                        tempuy2=0.0f;
                        tempuz2=0.0f;
                        tempvx2=0.0f;
                        tempvy2=0.0f;
                        tempvz2=0.0f;
                        tempwx2=0.0f;
                        tempwy2=0.0f;
                        tempwz2=0.0f;
                        tempuxz=0.0f;
                        tempuxy=0.0f;
                        tempvyz=0.0f;
                        tempvxy=0.0f;
                        tempwxz=0.0f;
                        tempwyz=0.0f;
                        for(kk=1;kk<=5;kk++)
                        {
                            tempux2=tempux2+c[kk-1][0]*(u[k*ny*nx+j*nx+(i+kk)]+u[k*ny*nx+j*nx+(i-kk)]);
                            tempuy2=tempuy2+c[kk-1][0]*(u[k*ny*nx+(j+kk)*nx+i]+u[k*ny*nx+(j-kk)*nx+i]);
                            tempuz2=tempuz2+c[kk-1][0]*(u[(k+kk)*ny*nx+j*nx+i]+u[(k -kk)*ny*nx+j*nx+i]);
                            
                           
                           
                            tempvx2=tempvx2+c[kk-1][0]*(v[k*ny*nx+j*nx+(i+kk)]+v[k*ny*nx+j*nx+(i-kk)]);
                            tempvy2=tempvy2+c[kk-1][0]*(v[k*ny*nx+(j+kk)*nx+i]+v[k*ny*nx+(j-kk)*nx+i]);
                            tempvz2=tempvz2+c[kk-1][0]*(v[(k+kk)*ny*nx+j*nx+i]+v[(k-kk)*ny*nx+j*nx+i]);

                            tempwx2=tempwx2+c[kk-1][0]*(w[k*ny*nx+j*nx+(i+kk)]+w[k*ny*nx+j*nx+(i-kk)]);
                            tempwy2=tempwy2+c[kk-1][0]*(w[k*ny*nx+(j+kk)*nx+i]+w[k*ny*nx+(j-kk)*nx+i]);
                            tempwz2=tempwz2+c[kk-1][0]*(w[(k+kk)*ny*nx+j*nx+i]+w[(k-kk)*ny*nx+j*nx+i]);

                        } //for(kk=1;kk<=5;kk++) end

                        tempux2=(tempux2+c0*u[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
                        tempuy2=(tempuy2+c0*u[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
                        tempuz2=(tempuz2+c0*u[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;

                        tempvx2=(tempvx2+c0*v[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
                        tempvy2=(tempvy2+c0*v[k*ny*nx+j*nx+i])*vvp2*dtx*dtx;
                        tempvz2=(tempvz2+c0*v[k*ny*nx+j*nx+i])*vvs2*dtz*dtz;

                        tempwx2=(tempwx2+c0*w[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
                        tempwy2=(tempwy2+c0*w[k*ny*nx+j*nx+i])*vvs2*dtx*dtx;
                        tempwz2=(tempwz2+c0*w[k*ny*nx+j*nx+i])*vvp2*dtz*dtz;

                        for(kk=1;kk<=5;kk++)
                        {
                            for(kkk=1;kkk<=5;kkk++)
                            {
                                tempuxz=tempuxz+c[kkk-1][1+kk]*(u[(k+kkk)*ny*nx+j*nx+(i+kk)]
                                                   -u[(k-kkk)*ny*nx+j*nx+(i+kk)]
                                                   +u[(k-kkk)*ny*nx+j*nx+(i-kk)]
                                                   -u[(k+kkk)*ny*nx+j*nx+(i-kk)]);
                                tempuxy=tempuxy+c[kkk-1][1+kk]*(u[k*ny*nx+(j+kkk)*nx+(i+kk)]
                                                   -u[k*ny*nx+(j-kkk)*nx+(i+kk)]
                                                   +u[k*ny*nx+(j-kkk)*nx+(i-kk)]
                                                   -u[k*ny*nx+(j+kkk)*nx+(i-kk)]);

                                tempvyz=tempvyz+c[kkk-1][1+kk]*(v[(k+kkk)*ny*nx+(j+kk)*nx+i]
                                                   -v[(k-kkk)*ny*nx+(j+kk)*nx+i]
                                                   +v[(k-kkk)*ny*nx+(j-kk)*nx+i]
                                                   -v[(k+kkk)*ny*nx+(j-kk)*nx+i]);
                                tempvxy=tempvxy+c[kkk-1][1+kk]*(v[k*ny*nx+(j+kkk)*nx+(i+kk)]
                                                   -v[k*ny*nx+(j-kkk)*nx+(i+kk)]
                                                   +v[k*ny*nx+(j-kkk)*nx+(i-kk)]
                                                   -v[k*ny*nx+(j+kkk)*nx+(i-kk)]);

                                tempwyz=tempwyz+c[kkk-1][1+kk]*(w[(k+kkk)*ny*nx+(j+kk)*nx+i]
                                                   -w[(k-kkk)*ny*nx+(j+kk)*nx+i]
                                                   +w[(k-kkk)*ny*nx+(j-kk)*nx+i]
                                                   -w[(k+kkk)*ny*nx+(j-kk)*nx+i]);
                                tempwxz=tempwxz+c[kkk-1][1+kk]*(w[(k+kkk)*ny*nx+j*nx+(i+kk)]
                                                   -w[(k-kkk)*ny*nx+j*nx+(i+kk)]
                                                   +w[(k-kkk)*ny*nx+j*nx+(i-kk)]
                                                   -w[(k+kkk)*ny*nx+j*nx+(i-kk)]);
                                   
                            } // for(kkk=1;kkk<=5;kkk++) end
                        } //for(kk=1;kk<=5;kk++) end
                         
                        up[k*ny*nx+j*nx+i]=2.*up1[k*ny*nx+j*nx+i]-up2[k*ny*nx+j*nx+i]
                                          +tempux2+tempwxz*vvp2*dtz*dtx
                                          +tempvxy*vvp2*dtz*dtx;
 
                        vp[k*ny*nx+j*nx+i]=2.*vp1[k*ny*nx+j*nx+i]-vp2[k*ny*nx+j*nx+i]
                                          +tempvy2+tempuxy*vvp2*dtz*dtx
                                          +tempwyz*vvp2*dtz*dtx;
                        wp[k*ny*nx+j*nx+i]=2.*wp1[k*ny*nx+j*nx+i]-wp2[k*ny*nx+j*nx+i]
                                          +tempwz2+tempuxz*vvp2*dtz*dtx
                                          +tempvyz*vvp2*dtz*dtx
                                          +((i == ncx_shot - 1 && j == ncy_shot - 1 && k == ncz_shot - 1) ? wave[l-1] : 0);
                        us[k*ny*nx+j*nx+i]=2.*us1[k*ny*nx+j*nx+i]-us2[k*ny*nx+j*nx+i]+tempuy2+tempuz2
                                          -tempvxy*vvs2*dtz*dtx-tempwxz*vvs2*dtz*dtx;
                        vs[k*ny*nx+j*nx+i]=2.*vs1[k*ny*nx+j*nx+i]-vs2[k*ny*nx+j*nx+i]+tempvx2+tempvz2
                                          -tempuxy*vvs2*dtz*dtx-tempwyz*vvs2*dtz*dtx;
                        ws[k*ny*nx+j*nx+i]=2.*ws1[k*ny*nx+j*nx+i]-ws2[k*ny*nx+j*nx+i]+tempwx2+tempwy2
                                          -tempuxz*vvs2*dtz*dtx-tempvyz*vvs2*dtz*dtx;
                    }//for(i=nleft;i<nright;i++) end

     
            for(k=ntop;k<nbottom;k++)
                for(j=nfront;j<nback;j++)
                    for(i=nleft;i<nright;i++)
                    {
                        u[k*ny*nx+j*nx+i]=up[k*ny*nx+j*nx+i]+us[k*ny*nx+j*nx+i];
                        v[k*ny*nx+j*nx+i]=vp[k*ny*nx+j*nx+i]+vs[k*ny*nx+j*nx+i];
                        w[k*ny*nx+j*nx+i]=wp[k*ny*nx+j*nx+i]+ws[k*ny*nx+j*nx+i];

                    }//for(i=nleft;i<nright;i++) end
        }//for(l=1;l<=lt;l++) end
        
        memcpy(up2, up1, size);
        memcpy(up1, up, size);
        memcpy(us2, us1, size);
        memcpy(us1, us, size);
        memcpy(vp2, vp1, size);
        memcpy(vp1, vp, size);
        memcpy(vs2, vs1, size);
        memcpy(vs1, vs, size);
        memcpy(wp2, wp1, size);
        memcpy(wp1, wp, size);
        memcpy(ws2, ws1, size);
        memcpy(ws1, ws, size);
   /*
  写入文件的数据
*/
        for(k=ntop;k<nbottom;k++)
                for(j=nfront;j<nback;j++)
                    for(i=nleft;i<nright;i++){
              
       fwrite(&up[k*ny*nx+j*nx+i], sizeof(float), 1, fout);
}
         
    }//for(ishot=1;ishot<=nshot;ishot++) end
    fclose(fout);

    free(u);
    free(v);
    free(w);
    free(up);
    free(up1);
    free(up2);
    free(vp);
    free(vp1);
    free(vp2);
    free(wp);
    free(wp1);
    free(wp2);
    free(us);
    free(us1);
    free(us2);
    free(vs);
    free(vs1);
    free(vs2);
    free(ws);
    free(ws1);
    free(ws2);
    free(initial);
    free(wave);

    return 0;
}
