
#include <stdio.h> /*include standard c header (cabeça~ho padrão da linguagem C)*/
#include <process.h>
#include <ctype.h>
#include <math.h>
#define NCLE 0
#define NCOU 8
#define NPT 1500
#define NDEL 0
#define NPAR 5
#define NTOT NCOU *NPT

float pi = 3.14159;
float pid = 6.28318;
float hor = 0.000000, te = 0.000000;
float tmax = 0.000000, tmin = 0.000000;
float gt[3], xs[10], x[8], ks[10][2], isi[3], is[3], dei[4], vf[3];
float vsa = 1., vsb = -0.5, vsc = -0.5, angf, iex = 1.13, vfi = 0.;
float isa = 0., isb = 0., isd = 0., isq = 0.;
float isdi = 0.0, isqi = 0.0, isai = 0.0, isbi = 0.0, fsai = 0.0, fsbi = 0.0;
float frdi = 0.0, frqi = 0.0, fsdi = 0.0, fsqi = 0.0;
float fsa, ce, wr, isab, isbb, ceb, rrv = 0;

/* Dados de placa da maquina */
/* Maquina 1 */
float rr = 3.421, rs = 5.793, p = 2.0, ls = 0.386, lr = 0.386, lm = 0.363;
float jj = 0.0267, ff = 0.0297, cm = 0.0, ceref = 0.0;

/* */
/* Maquina 2 */
/*float rr=1.410, rs=0.390, p= 2.0, ls=0.094, lr=0.094, lm=0.091;
float jj=0.0400, ff=0.0100, cm=0.0, ceref=0.0;*/
/* */
float fsdf = 0.0, fsqf = 0.0, frdf = 0.0, frqf = 0.0;
float fsaf = 0.0, fsbf = 0.0, fraf = 0.0, frbf = 0.0;
float isaf = 0.0, isbf = 0.0, isdf = 0.0, isqf = 0.0;
float tetf = 0.0, wrf = 0., cef = 0.;
float vslf = 0.0, vs2f = 0.0, vs3f = 0.0;
float cmf = 0., wd = 314.159, psi = 0.;
double fldbl, fldb2;
float fls = 0.8, flr = 0.8, ikd, ikq, flsb, flrb, qrs;
float fsb, fra, frb, frd, frq, fed, fsq, flrr = 0.8, flsr = 0.88, rlmd;
float fsdb, fsqb, vsd, vsq, war, wa, wbr, wbi, angbi = 0.0, angai = 0.0;
float tm1 = 10.00, tm1a = 15., tm1b = 0.00, tm2 = 0., tm3 = 0.;
float rq2 = 1.414214, rq3 = 1.732051, rq23 = 0.816496;
/* Dados dos reguladores*/
float kptot = 0., xerf = 0.0, xerw = 0.0, erf = 0.0, erw = 0.0;
float fkp = 0., fki = 0., fkp2 = 0., fki2 = 0.;
float ekid = 0., ekiq = 0.;
float anga = 0., ctt, stt;
float hinv = 0.0;
float dkr = 0.0, dka = 0.0;
/* */
float cgbr = 0.0, sgbr = 0.0, frdr = 0.0, frqr = 0.0, cgbi = 0.0;
float sigma = 0.0, istl = 0.0;
float vsd = 0., vsq = 0., cel = 0., cd = 0.0;
float vsar = 0.0, vsbr = 0.0;
float aux[3], gtf[3],xf[7],xd[5][7],dv[5][7],xa[5][73; 
float
 esd1=0.0,esq1=0.0,esq2=0.0,esq3=0.0,ts=0.0,tr=0.0; 
float Is, Vsl,Vs2,Vs3,Vs; 
int kmod= 1, klim=0, fonte=0, fcorr=0, regul=0; 
FILE *stream; 
int label[NPAR]; 
float abscisse[NPT]; 
float coordinate[NCOU][NPT]; 

int ia; 
long offset; 
int numread; 
int numwritten; 
float t = 0.0, tm10 = 0.00, hm10 = 0.0, np2 = 0.0; 
int j,id,ip,np,i2,g,g2; 

main() 
{
    /* */
    id = -1;
    ip = 0;
    np = 0;
    offset = 4 * NPT;
    label[0] = NCLE;
    label[1] = NCOU;
    label[2] = NPT;
    label[3] = NDEL;
    label[4] = NDEL;
    /**/
    printf("\n Entre periodo amostragem (te):");
    scanf("%f", &te);
    printf("\n Entre tempo maximo simulacao'(tmax): ");
    scanf("%f", &tmax);
    printf("\n Entre passo de integracao (hor) :");
    scanf("%f", &hor);
    /*******************************************************/
    /*Não tem a opção kmod==2 Fluxo Rotorico escorregamento estator*/
    /********************************************/
    /* */
    sigma = 1. - (lm * lm / ls / lr);
    tr = lr / rr;
    ts = ls / rs;
    hinv = tr / lm / te;
    cd = 1. / ls / lr / sigma;
    esd1 = (lm / sigma / ts / lr);
    esq1 = lm / lr;
    esq2 = lm / tr / lr;
    esq3 = sigma * ls;

    /* */
    /*Inicialização de parâmetros*/
    /* FILE *fp_out: */
    /* Calculo do passo de escritura hm10 */
    hm10 = (tmax - tm10) / NPT;
    if (hm10 < hor)
        hm10 = hor;

    /* / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / / */

    /* INICIO DE PROCESSAMENTO */

    while (t < tmax)
    {
        if (t >= tm2)
        {
            bloco2();
        }
        if (t >= tm1b)
        {
            bloc1b();
        }
        if ((t >= tm10) && (id < (NPT - 1)))
        {
            id = id + 1;
            tm10 = tm10 + hm10;
            ip = ip + 1;
            if (ip >= 50)
            {
                ip = 0;
                printf("%f\n", t);
            }
            abscisse[id] = tm10 - hm10;
            coordinate[0][id] = (float)(sqrt(isa * isa + isb * isb));
            coordinate[1][id] = (float)(sqrt(fsa * fsa + fsb * fsb));
            coordinate[2][id] = (float)(sqrt(fra * fra + frb * frb));
            coordinate[3][id] = ce;
            coordinate[4][id] = ceref;
            coordinate[5][id] = wr;
            coordinate[6][id] = vsar;
            coordinate[7][id] = vsa;
        }
        t = t + hor;
    }
    printf("%d\n", id);
    printf("%f %f\n", "war = ", war);

    /*File opened in binary (b) mode */
    if ((stream = fopen("invtd.des", "w+b")) != NULL)
    {
        /*Writing labels*/
        numwritten = fwrite((int *)label, sizeof(int), NPAR, stream);
        fseek(stream, offset, SEEK_SET);
        printf("Wrote %d items\n", numwritten);
        /* Writing independent variable */
        numwritten = fwrite((float *)abscisse, sizeof(float), NPT, stream);
        printf("Wrote %d items\n", numwritten);
        /* Writing dependent variables */
        numwritten = fwrite((float *)coordinate, sizeof(float), NTOT, stream);
        printf("Wrote %d items\n", numwritten);
    }
    else
    printf("Problem opening the file");
    fclose(stream);
    /* Attempt to read in NPT long integers */
    if ((stream = fopen("invtd.des", "r+b")) != NULL)
    {
        numread = fread((int *)label, sizeof(int), NPT, stream);
        printf("Number of items read = %d\n", numread);
        fclose(stream);
    }
    else
    {
        printf("Was not able to open the file");
    }

    bloc1b()
    {
        /* modulo controle de conjugado e fluxo estatorico */
        /* utilizando fonte PWM para conj. e fluxo */
        float hm1b;
        static float tck = 0., tk = 0, tk1 = 0, to = 0.;
        static float ec = 900.;
        static float cont = 0.;
        static float vonf, vtf = 0.;
        static int se, kst2, kor = 0.;
        static float ist = 0.;
        int i;
        static float vsd, vsq;
        static float cga, sga, cgb, sgb, angb = 0.;
        /**/
        hm1b = hor;
        tm1b = tm1b + hm1b;
        /**/
        vsa = velf;
        vsb = ve2f;
        vsc = vs3f;
        gt[0] = gtf[0];
        gt[1] = gtf[1];
        gt[2] = gtf[2];
        /**/
        if (tm1b > 0.03)
            ceref = 7.52;
        if (tm1b > 0.15)
            ceref = -7.52;
        ist = ceref * lr / lm / flrr / 2.;
        war = ceref * ls * ls * rr / p / flsr / flsr / lm / lm;
        wbr = lm * rr * ist / flrr / lr;
        /**/
        cont = cont + hm1b;

        if (cont >= te)
        {
            /* amostragem: medicao e calculo dos reguladores e PWM */
            /* calculo modulos e fases dos fluxos rotoricos e estatoricos */
            fldb1 = (double)(x[0] * x[0]);
            fldb2 = (double)(x[1] * x[1]);
            fldb1 = sqrt(fldb1 + fldb2);
            fls = (float)(fldb1);
            fldb1 = (double)(x[2] * x[2]);
            fldb2 = (double)(x[3] * x[3]);
            fldb1 = sqrt(fldb1 + fldb2);
            flr = (float)(fldb1);
            cga = fsa / fls;
            sga = fsb / fls;
            cgb = fra / flr;
            sgb = frb / flr;
            /* */

            /***************************************************************/
            /* CONTROLE DE FLUXO ROTÓRICO (escorregamento estator kmod=2) */
            /***************************************************************/

            wbi = wbr + wr;
            angbi = angbi + wbi * te;
            if (angbi >= pid)
                angbi = angbi - pid;
            /* */
            frdi = flrr * cos(angbi);
            frqi = flrr * sin(angbi);
            /* */
            frd = fra;
            frq = frb;
            /* */
            if (regul == 0)
            {
                wbi = wbr + wr;
                angbi = angbi + wbi * te;
                if (angbi >= pid)
                    angbi = angbi - pid;
                fsdi = flrr * cos(angbi);
                fsqi = flrr * sin(angbi);
                vsd = -(wbi)*fsqi;
                vsq = (wbi)*fsdi;
                vsar = vsd;
                vsbr = vsq;
            }
            /* */
            if (regul == 1)
            {
                dei[2] = frdi - frd;
                dei[3] = frqi - frq;
                /* */
                isdi = erf + (fkp + fki) * dei[2];
                erf = erf + fki * dei[2];
                if (erf >= 15.)
                    erf = 15;
                if (erf <= -15.)
                    erf = -15.;
                isdi = isdi + (wr * frd * tr / lm);
                /* */
                isqi - erw + (fkp + fki) * dei[3];
                erw = erw + fki * dei[3];
                if (erw >= 15.)
                    erw = 15.;
                if (erw <= -15.)
                    erw = -15.;
                isqi = isqi - (wr * frq * tr / lm);
            }
            if (regul == 2)
            {
                /* */
                /* Regulador preditivo */
                /* */
                isdi = hinv * frdi - hinv * (1. - te / tr) * frd + hinv * wr * te * frq;
                isqi = hinv * frqi - hinv * wr * te * frd - hinv * (1. - te / tr) * frq;
                /*printf("%s %f %s %f\n","hinv =",hinv,"frdi =",frdi);
                printf("%s %f %s %f\n","te = ",te,"tr =",tr);
                printf("%s %f %s %f\n","frd =",frd,"frq =",frq);*/
            }
            /* */
            if (isdi >= 8.)
                isdi = 8.;
            if (isdi <= -8.)
                isdi = -8.;
            if (isqi >= 8.)
                isqi = 8.;

            if (isqi <= -8.)
                isqi = -8.;
            /* */
            isai = isdi;
            isbi = isqi;
            /* */
            /* Fonte de corrente no estator */
            /* */
            dei[0] = isai - isa;
            dei[1] = isbi - isb;
            vsd = ekid + (dkr + dka) * dei[0] - (wr + wbr) * lm * flr * sgb / lr;
            vsq = ekiq + (dkr + dka) * dei[1] + (wr + wbr) * lm * flr * cgb / lr;
            ekid = ekid + dka * dei[0];
            ekiq = ekiq + dka * dei[1];
            vsar = vsd;
            vsbr = vsq;
            /**/
            if (vsar >= 311.)
                vsar = 311.;
            if (vsar <= -311.)
                vsar = -311.;
            if (vsbr >= 311.)
                vsbr = 311.;
            if (vsbr <= -311.)
                vsbr = -311;
            // break;
            /*
            angb = angb + (wr+wbr)*te;
            if(angb >= pid) angb = angb - pid;
            cgb=cos(angb); sgb=sin(angb);
            frd = fra*cgb + frb*sgb;
            frq = frb*cgb - fra*sgb;
            dei[2] = f l r r - frd;
            vsd = erf + (fkp+fki)*dei[2];
            erf = erf + fki*dei[2];
            isdi = vsd/(lm*rr);*/
            /* */
            /*
             isqi=ist;
            derw = erw + (fkp+fki)*(-frq);
            erw = erw + fki*(-frq);
            wbr = lm*rr*isqi/( frd*lr) - derw/frd;*/
            /* mudança de coordenadas */
            /*isai = isdi*cgb - isqi*sgb;
            isbi = isqi*cgb + isdi*sgb;
            */

            /********************************************/
            /* Fonte PWM */
            /********************************************/
            /* calculo do modulador PWM */
            /********************************************/

            cont = 0.;
            klim = 0;
            kor = kor + 1.;
            if (kor >= 2)
                kor = 0;
            tck = 2. * te / rq3;

            /********************************************/
            /* Determinacao do setor vetor tensao de operacao */
            vsar = vsar / rq23;
            vsbr = vsbr / rq23;
            aux[0] = rq3 * vsar + vsbr;
            aux[1] = -rq3 * vsar + vsbr;
            aux[2] = -2. * vsbr;

            /**/
            if ((aux[0] >= 0.) && (aux[1] <= 0.) && (aux[2] <= 0.))
                se = 1;
            else
            {
                if ((aux[0] >= 0) && (aux[1] >= 0) && (aux[2] <= 0))
                    se = 2;
                else
                {
                    if ((aux[0] <= 0) && (aux[1] >= 0) && (aux[2] <= 0))
                        se = 3;
                    else
                    {
                        if ((aux[0] <= 0) && (aux[1] >= 0) && (aux[2] >= 0))
                            se = 4;
                        else
                        {
                            if ((aux[0] <= 0) && (aux[1] <= 0) && (aux[2] >= 0))
                                se = 5;
                            else
                            {
                                if ((aux[0] >= 0) && (aux[1] <= 0) && (aux[2] >= 0))
                                    se = 6;
                            }
                        }}}
            }

            /* calculo dos tempos de aplicacao dos vetores tensao */

            tk = tck * (vsar * sin(se * pi / 3.) - vsbr * cos(se * pi / 3.)) / ec;
            tk1 = -tck * (vsar * sin((se - 1) * pi / 3.) - vsbr * cos((se - 1) * pi / 3)) / ec;
            if (tk <= 0.)
                tk = 0.;
            if ((tk >= te) && (tk >= tk1))
            {
                kst2 = se;
                klim = 1;
            }
            if ((tk1 >= te) && (tk1 >= tk))
            {
                kst2 = se + 1; klim = 1; if (kst2 > 6) kst2 = 1;  }

            if (tk1 <= 0.)
                tk1 = 0.;
            to = te - tk - tk1;
            if (to >= te)
                to = te;
            if (to <= 0.)
                to = 0.;

            /* */
            /* Determinacao dos vetores tensao instantaneos (funcao timer) */
            /* */

            if (klim == 0)
            {
                if (cont <= to / 2.)
                    kst2 = 7;
                else
                {
                    if (kor == 0)
                    {
                        if (cont < (to / 2. + tk))
                            kst2 = se;
                        else
                        {
                            if (cont < (to / 2. + tk + tk1))
                            {
                                kst2 = se + 1;
                                if (kst2 > 6)
                                    kst2 = 1;
                            }
                            else
                            {
                                kst2 = 7;
                            }
                        }
                    }
                    else
                    {
                        if (cont < (to / 2. + tk1))
                        {
                            kst2 = se + 1;
                            if (kst2 > 6)
                                kst2 = 1;
                        }
                        else
                        {
                            if (cont < (to / 2. + tk + tk1))
                                kst2 = se;
                                else { kst2 = 7; }
                        }
                    }
                }
            }
switch(kst2) 
{ 
case 1: gtf[0] =1; gtf[1] =0; gtf[2] =0; vtf= 1 .;
break; 
case 2: gtf[0] =1; gtf[1] = 1; gtf[2] =0; vtf= 2.;
break; 
case 3: gtf[0] =0; gtf[1] = 1; gtf[2] =0; vtf= 3.;
break; 
case 4: gtf[0] =0; gtf[1] =1; gtf[2] =1; vtf= 4.;
break; 
case 5: gtf[0] =0; gtf[1] =0; gtf[2] =1; vtf= 5.;
break; 
case 6: gtf[0] =1; gtf[1] =0; gtf[2] = 1; vtf= 6.;
break; 
case 7: vtf=0.; 
if(gt[0] == gt[1]) 
gtf[2] = gt[0]; 
else {
if(gt[0] == gt[2]) 
gtf[1]=gt[2]; 
else { gtf[0] = gt[1]; }
}
}

/*calculo das variáveis de saida **/ 
/* */ 
vf[0]=gtf[0]*(ec/2) + (gtf[0]-1)*(ec/2) ;
vf[1]=gtf[1]*(ec/2) + (gtf[1]-1)*(ec/2) ;
vf[2]=gtf[2]*(ec/2) + (gtf[2]-1)*(ec/2) ;
/* */ 
vonf=(1./3.)*(vf[0]+vf[l]+vf[2]); 
vslf=vf[0]-vonf; 
vs2f=vf[1]-vonf; 
vs3f=vf[2]-vonf; 
/* */ 
} 
bloco2() 
{ 
float hm2; 
static float vrd=0., vrq=0.; 
/* */
static int ki3=0, ne=5; 
int 1,j; 
hm2 = hor; 
tm2 = tm2 + hm2; 
rlmd = rs*cd*lm; 
/* */
isa=isaf; 
isb=isbf; 
isd = isdf;
isq = isqf; 
wr=wrf; 
ce=cef; 
fsd = fsdf; fsq = fsqf; 
frd = frdf; frq = frqf; 
fsa = fsaf; fsb = fsbf; 
fra = fraf; frb = frbf; 
if ( t >= tmin) 
{ 
printf("%s %f %s %f\n"," wr = ",wr," ce = ",cef); 
tmin = tmin + 0.01; 
}

/* */ 
psi = psi + wd*hm2; 
if(psi >= pid) psi=psi-pid; 
/* */
ved=rq23*(vsa*cos(psi)+vsb*cos(psi-pid/3.)+vsc*cos(psi+pid/3.)); 
vsq=-rq23*(vsa*sin(psi)+vsb*sin(psi-pid/3.)+vsc*sin(psi+pid/3.)); 
/* */ 
/* Integração: runge kutta 4a ordem */ 
/* */ 
x[0]=fsdf; x[1]=fsqf; x[2] frdf; x[3]=frqf; x[4]=tetf; x[5]=wrf; 
for(i=0; i<=ne; ++i)
{
xd[0][i]= x[i];
}
for(j=0; j<=3; ++j)
{
dv[j][0]=vsd-(cd*rs)*(xd[j][0]*lr-xd[j][2]*lm)+wd*xd[j][1]; 
dv[j][1]=vsq-(cd*rs)*(xd[j][1]*lr-xd[j][3]*lm)-wd*xd[j][0];
dv[j][2]=vrd-(cd*rr)*(xd [j][2]*ls-xd[j][0]*lm)+(wd-xd[j][5])
*xd[j][3];
dv[j][3]=vrq-(cd*rr)*(xd[j][3]*ls-xd[j][1]*lm)-(wd-xd[j][5])
*xd[j][2];
dv[j][4]=xd[j][5];
cel = p*cd*lm* (xd [j][1]*xd[j][2]-xd[j][0]*xd[j][3]); 
dv[j][5] = (p*cel-ff*xd [j][5]-p*cm)/jj;
for(i=0; i<=ne; ++i)
    {
     xa[j][i] hm2*dv[j][i]; 
     switch(j)
        {
            case 0: xd[1][i] = x[i] + xa[0][i]/2.;
            break;
            case 1: xd[2][i] = x[i] + xa[1][i]/2.;
            break;
            case 2: xd [3][i] = x[i] + xa[2][i];
            }
        }
    }
for(i=0; i<=ne; ++i)
{
xf[i] = x[i]+xa[0][i]/6.+xa [1] [i]/3.+xa [2][i]/3.+xa[3][i]/6.;
}
fsdf=xf[0]; fsqf=xf[1]; frdf=xf[2]; frqf=xf[3];

/* */
tetf=xf[4]; wrf=xf[5];
if(tetf>= pid) tetf=tetf-pid;
isdf = cd*(xf[0]*lr-xf[2]*lm);
isqf = cd*(xf[1]*lr-xf[3]*lm);
/* */
isaf = (isdf*cos(psi)-isqf*sin(psi));
isbf (isqf*cos(psi)+isdf*sin(psi));
fsaf = (fsdf*cos(psi)-fsqf*sin(psi));
fsbf = (feqf*cos(psi)+fsdf*sin(psi));
fraf = (frdf*cos(psi)-frqf*sin(psi));
frbf (frqf*cos(psi)+frdf*sin(psi));
cef = p*cd*lm* (xf[1]*xf[2]-xf[0]*xf[3]);
ki3 = ki3 + 1;
if(ki3>=1000)
{
ki3 0;
}
}