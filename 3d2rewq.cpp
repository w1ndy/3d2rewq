#include <stdio.h>
//#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"

#define PIE 3.1415926
typedef unsigned char BYTE;

void memdup(BYTE *dst, BYTE *src, size_t dst_size, size_t src_size)
{
    if(dst_size < src_size) memcpy(dst, src, dst_size);
    int len = src_size, remain = dst_size - src_size;
    memcpy(dst, src, len);
    while(remain > len) {
        memcpy(dst + len, dst, len);
        remain -= len;
        len <<= 1;
    }
    memcpy(dst + len, dst, remain);
}

int main(int argc, char **argv)
{
    int nx,ny,nz,lt,nedge;
    float frequency;
    float velmax;
    float dt;
    int ncx_shot1,ncy_shot1,ncz_shot;
    int ishot, nshot;
    float unit;
    int nxshot,nyshot,dxshot,dyshot;
    FILE  *fin, *fout;

    float *initial, *wave;
    float c[5][7];
    float t0,tt,c0;
    float dtx,dtz,dr1,dr2;

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

    const size_t initial_count = 64;
    const size_t initial_size = initial_count * sizeof(float);
    initial = (float *)malloc(initial_size);
    wave = (float*)malloc(sizeof(float)*lt);

    nshot=nxshot*nyshot;
    t0=1.0/frequency;

    for(int l=0;l<lt;l++)
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

    for(int i=0;i<5;i++)
        for(int j=0;j<5;j++)
            c[j][2+i]=c[i][1]*c[j][1];

    dtx=dt/unit;
    dtz=dt/unit;

    dr1=dtx*dtx/2.;
    dr2=dtz*dtz/2.;

    fout=fopen(argv[2],"wb");

    for(int i = 0; i < initial_count; i++)
        initial[i] = 100.0f;

#pragma omp parallel for schedule(static,1) ordered num_threads(4)
    for(ishot=1;ishot<=nshot;ishot++)
    {
        int nleft,nright,nfront,nback,ntop,nbottom;
        float *u, *v, *w, *up, *up1, *up2,
              *vp, *vp1, *vp2, *wp, *wp1, *wp2, 
              *us, *us1, *us2, *vs, *vs1, *vs2,
              *ws, *ws1, *ws2, xmax;
        float vvp2,vvs2,tempux2,tempuy2,tempuz2,tempvx2,tempvy2,tempvz2,
              tempwx2,tempwy2,tempwz2,tempuxz,tempuxy,tempvyz,tempvxy,tempwxz,tempwyz;
        
        printf("shot=%d\n",ishot);
        int ncy_shot=ncy_shot1+(ishot/nxshot)*dyshot;
        int ncx_shot=ncx_shot1+(ishot%nxshot)*dxshot;
/*
此处是初始化参数，以前是0，后来改为100，也可改为其他或者随机赋值 ,up up1 up2可不修改
*/

        // alloc mem:
        xmax = lt * dt * velmax;
        int left_max    = ncx_shot - xmax / unit - 10;
        int right_max   = ncx_shot + xmax / unit + 10;
        int front_max   = ncy_shot - xmax / unit - 10;
        int back_max    = ncy_shot + xmax / unit + 10;
        int top_max     = ncz_shot - xmax / unit - 10;
        int bottom_max  = ncz_shot + xmax / unit + 10;

        int xwidth = right_max - left_max + 20;
        int ywidth = back_max - front_max + 20;
        int zwidth = bottom_max - top_max + 20;
        int ecount = xwidth * ywidth * zwidth;
        int block_size = sizeof(float) * ecount;

        float *all = (float*)malloc(21 * block_size);
        memdup((BYTE *)all, (BYTE *)initial, 21 * block_size, initial_size);

        u       = all;
        v       = all + ecount;
        w       = all + ecount * 2;
        up      = all + ecount * 3;
        up1     = all + ecount * 4;
        up2     = all + ecount * 5;
        vp      = all + ecount * 6;
        vp1     = all + ecount * 7;
        vp2     = all + ecount * 8;
        wp      = all + ecount * 9;
        wp1     = all + ecount * 10;
        wp2     = all + ecount * 11;
        us      = all + ecount * 12;
        us1     = all + ecount * 13;
        us2     = all + ecount * 14;
        vs      = all + ecount * 15;
        vs1     = all + ecount * 16;
        vs2     = all + ecount * 17;
        ws      = all + ecount * 18;
        ws1     = all + ecount * 19;
        ws2     = all + ecount * 20;

        memset(up, 0, 3 * block_size);

        for(int l=1;l<=lt;l++)
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
          
            ntop -= top_max - 10;
            nbottom -= top_max - 10;
            nfront -= front_max - 10;
            nback -= front_max - 10;
            nleft -= left_max - 10;
            nright -= left_max - 10;

            for(int k=ntop;k<nbottom;k++)
                for(int j=nfront;j<nback;j++)
                    for(int i=nleft;i<nright;i++)
                    {
                        if(k < 220 - top_max) {
                            vvp2 = 2300 * 2300;
                            vvs2 = 1232 * 1232;
                        } else if(k < 270 - top_max) {
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
                        for(int kk=1;kk<=5;kk++)
                        {
                            tempux2=tempux2+c[kk-1][0]*(u[k*ywidth*xwidth+j*xwidth+(i+kk)]+u[k*ywidth*xwidth+j*xwidth+(i-kk)]);
                            tempuy2=tempuy2+c[kk-1][0]*(u[k*ywidth*xwidth+(j+kk)*xwidth+i]+u[k*ywidth*xwidth+(j-kk)*xwidth+i]);
                            tempuz2=tempuz2+c[kk-1][0]*(u[(k+kk)*ywidth*xwidth+j*xwidth+i]+u[(k -kk)*ywidth*xwidth+j*xwidth+i]);
                           
                            tempvx2=tempvx2+c[kk-1][0]*(v[k*ywidth*xwidth+j*xwidth+(i+kk)]+v[k*ywidth*xwidth+j*xwidth+(i-kk)]);
                            tempvy2=tempvy2+c[kk-1][0]*(v[k*ywidth*xwidth+(j+kk)*xwidth+i]+v[k*ywidth*xwidth+(j-kk)*xwidth+i]);
                            tempvz2=tempvz2+c[kk-1][0]*(v[(k+kk)*ywidth*xwidth+j*xwidth+i]+v[(k-kk)*ywidth*xwidth+j*xwidth+i]);

                            tempwx2=tempwx2+c[kk-1][0]*(w[k*ywidth*xwidth+j*xwidth+(i+kk)]+w[k*ywidth*xwidth+j*xwidth+(i-kk)]);
                            tempwy2=tempwy2+c[kk-1][0]*(w[k*ywidth*xwidth+(j+kk)*xwidth+i]+w[k*ywidth*xwidth+(j-kk)*xwidth+i]);
                            tempwz2=tempwz2+c[kk-1][0]*(w[(k+kk)*ywidth*xwidth+j*xwidth+i]+w[(k-kk)*ywidth*xwidth+j*xwidth+i]);

                        } //for(kk=1;kk<=5;kk++) end

                        tempux2=(tempux2+c0*u[k*ywidth*xwidth+j*xwidth+i])*vvp2*dtx*dtx;
                        tempuy2=(tempuy2+c0*u[k*ywidth*xwidth+j*xwidth+i])*vvs2*dtx*dtx;
                        tempuz2=(tempuz2+c0*u[k*ywidth*xwidth+j*xwidth+i])*vvs2*dtz*dtz;

                        tempvx2=(tempvx2+c0*v[k*ywidth*xwidth+j*xwidth+i])*vvs2*dtx*dtx;
                        tempvy2=(tempvy2+c0*v[k*ywidth*xwidth+j*xwidth+i])*vvp2*dtx*dtx;
                        tempvz2=(tempvz2+c0*v[k*ywidth*xwidth+j*xwidth+i])*vvs2*dtz*dtz;

                        tempwx2=(tempwx2+c0*w[k*ywidth*xwidth+j*xwidth+i])*vvs2*dtx*dtx;
                        tempwy2=(tempwy2+c0*w[k*ywidth*xwidth+j*xwidth+i])*vvs2*dtx*dtx;
                        tempwz2=(tempwz2+c0*w[k*ywidth*xwidth+j*xwidth+i])*vvp2*dtz*dtz;

                        for(int kk=1;kk<=5;kk++)
                        {
                            for(int kkk=1;kkk<=5;kkk++)
                            {
                                tempuxz += c[kkk-1][1+kk]*(u[(k+kkk)*ywidth*xwidth+j*xwidth+(i+kk)]
                                                   -u[(k-kkk)*ywidth*xwidth+j*xwidth+(i+kk)]
                                                   +u[(k-kkk)*ywidth*xwidth+j*xwidth+(i-kk)]
                                                   -u[(k+kkk)*ywidth*xwidth+j*xwidth+(i-kk)]);
                                tempuxy += c[kkk-1][1+kk]*(u[k*ywidth*xwidth+(j+kkk)*xwidth+(i+kk)]
                                                   -u[k*ywidth*xwidth+(j-kkk)*xwidth+(i+kk)]
                                                   +u[k*ywidth*xwidth+(j-kkk)*xwidth+(i-kk)]
                                                   -u[k*ywidth*xwidth+(j+kkk)*xwidth+(i-kk)]);

                                tempvyz += c[kkk-1][1+kk]*(v[(k+kkk)*ywidth*xwidth+(j+kk)*xwidth+i]
                                                   -v[(k-kkk)*ywidth*xwidth+(j+kk)*xwidth+i]
                                                   +v[(k-kkk)*ywidth*xwidth+(j-kk)*xwidth+i]
                                                   -v[(k+kkk)*ywidth*xwidth+(j-kk)*xwidth+i]);
                                tempvxy += c[kkk-1][1+kk]*(v[k*ywidth*xwidth+(j+kkk)*xwidth+(i+kk)]
                                                   -v[k*ywidth*xwidth+(j-kkk)*xwidth+(i+kk)]
                                                   +v[k*ywidth*xwidth+(j-kkk)*xwidth+(i-kk)]
                                                   -v[k*ywidth*xwidth+(j+kkk)*xwidth+(i-kk)]);

                                tempwyz += c[kkk-1][1+kk]*(w[(k+kkk)*ywidth*xwidth+(j+kk)*xwidth+i]
                                                   -w[(k-kkk)*ywidth*xwidth+(j+kk)*xwidth+i]
                                                   +w[(k-kkk)*ywidth*xwidth+(j-kk)*xwidth+i]
                                                   -w[(k+kkk)*ywidth*xwidth+(j-kk)*xwidth+i]);
                                tempwxz += c[kkk-1][1+kk]*(w[(k+kkk)*ywidth*xwidth+j*xwidth+(i+kk)]
                                                   -w[(k-kkk)*ywidth*xwidth+j*xwidth+(i+kk)]
                                                   +w[(k-kkk)*ywidth*xwidth+j*xwidth+(i-kk)]
                                                   -w[(k+kkk)*ywidth*xwidth+j*xwidth+(i-kk)]);
                                   
                            } // for(kkk=1;kkk<=5;kkk++) end
                        } //for(kk=1;kk<=5;kk++) end
                         
                        up[k*ywidth*xwidth+j*xwidth+i]=2.*up1[k*ywidth*xwidth+j*xwidth+i]-up2[k*ywidth*xwidth+j*xwidth+i]
                                          +tempux2+tempwxz*vvp2*dtz*dtx
                                          +tempvxy*vvp2*dtz*dtx;
 
                        vp[k*ywidth*xwidth+j*xwidth+i]=2.*vp1[k*ywidth*xwidth+j*xwidth+i]-vp2[k*ywidth*xwidth+j*xwidth+i]
                                          +tempvy2+tempuxy*vvp2*dtz*dtx
                                          +tempwyz*vvp2*dtz*dtx;
                        wp[k*ywidth*xwidth+j*xwidth+i]=2.*wp1[k*ywidth*xwidth+j*xwidth+i]-wp2[k*ywidth*xwidth+j*xwidth+i]
                                          +tempwz2+tempuxz*vvp2*dtz*dtx
                                          +tempvyz*vvp2*dtz*dtx
                                          +((i + left_max - 10 == ncx_shot - 1 && j + front_max - 10 == ncy_shot - 1 && k + top_max - 10 == ncz_shot - 1) ? wave[l-1] : 0);
                        us[k*ywidth*xwidth+j*xwidth+i]=2.*us1[k*ywidth*xwidth+j*xwidth+i]-us2[k*ywidth*xwidth+j*xwidth+i]+tempuy2+tempuz2
                                          -tempvxy*vvs2*dtz*dtx-tempwxz*vvs2*dtz*dtx;
                        vs[k*ywidth*xwidth+j*xwidth+i]=2.*vs1[k*ywidth*xwidth+j*xwidth+i]-vs2[k*ywidth*xwidth+j*xwidth+i]+tempvx2+tempvz2
                                          -tempuxy*vvs2*dtz*dtx-tempwyz*vvs2*dtz*dtx;
                        ws[k*ywidth*xwidth+j*xwidth+i]=2.*ws1[k*ywidth*xwidth+j*xwidth+i]-ws2[k*ywidth*xwidth+j*xwidth+i]+tempwx2+tempwy2
                                          -tempuxz*vvs2*dtz*dtx-tempvyz*vvs2*dtz*dtx;
                       
                    }//for(i=nleft;i<nright;i++) end

     
            for(int k=ntop;k<nbottom;k++)
                for(int j=nfront;j<nback;j++)
                    for(int i=nleft;i<nright;i++)
                    {
                        u[k*ywidth*xwidth+j*xwidth+i]=up[k*ywidth*xwidth+j*xwidth+i]+us[k*ywidth*xwidth+j*xwidth+i];
                        v[k*ywidth*xwidth+j*xwidth+i]=vp[k*ywidth*xwidth+j*xwidth+i]+vs[k*ywidth*xwidth+j*xwidth+i];
                        w[k*ywidth*xwidth+j*xwidth+i]=wp[k*ywidth*xwidth+j*xwidth+i]+ws[k*ywidth*xwidth+j*xwidth+i];

                        up2[k*ywidth*xwidth+j*xwidth+i]=up1[k*ywidth*xwidth+j*xwidth+i];
                        up1[k*ywidth*xwidth+j*xwidth+i]=up[k*ywidth*xwidth+j*xwidth+i];
                        us2[k*ywidth*xwidth+j*xwidth+i]=us1[k*ywidth*xwidth+j*xwidth+i];
                        us1[k*ywidth*xwidth+j*xwidth+i]=us[k*ywidth*xwidth+j*xwidth+i];
                        vp2[k*ywidth*xwidth+j*xwidth+i]=vp1[k*ywidth*xwidth+j*xwidth+i];
                        vp1[k*ywidth*xwidth+j*xwidth+i]=vp[k*ywidth*xwidth+j*xwidth+i];
                        vs2[k*ywidth*xwidth+j*xwidth+i]=vs1[k*ywidth*xwidth+j*xwidth+i];
                        vs1[k*ywidth*xwidth+j*xwidth+i]=vs[k*ywidth*xwidth+j*xwidth+i];
                        wp2[k*ywidth*xwidth+j*xwidth+i]=wp1[k*ywidth*xwidth+j*xwidth+i];
                        wp1[k*ywidth*xwidth+j*xwidth+i]=wp[k*ywidth*xwidth+j*xwidth+i];
                        ws2[k*ywidth*xwidth+j*xwidth+i]=ws1[k*ywidth*xwidth+j*xwidth+i];
                        ws1[k*ywidth*xwidth+j*xwidth+i]=ws[k*ywidth*xwidth+j*xwidth+i];
                    }//for(i=nleft;i<nright;i++) end

        }//for(l=1;l<=lt;l++) end
#pragma omp ordered
        {
            for(int k=ntop;k<nbottom;k++)
                for(int j=nfront;j<nback;j++)
                    fwrite(&up[k*ywidth*xwidth+j*xwidth+nleft], sizeof(float), nright - nleft, fout);
        }
        free(all);
    }//for(ishot=1;ishot<=nshot;ishot++) end
    fclose(fout);
    free(initial);
    free(wave);

    return 0;
}
