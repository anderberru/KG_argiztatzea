//	Program developed by
//	
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc *.c -lGL -lGLU -lglut -lm
//
// 
//


#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cargar-triangulo.h"

#define PI 3.14159265358979323846

typedef struct mlist
    {
    double m[16];
    struct mlist *hptr;
    } mlist;

typedef struct koordlist
    {
    double koord[3];
    struct koordlist *hptr;
    } koordlist;
    
typedef struct triobj
    {
    hiruki *triptr;
    int num_triangles;
    unsigned char *rgbptr;
    mlist *mptr;
    struct triobj *hptr;
    } triobj;

typedef struct argia
    {
    int piztu;
    double I[3];    // Intentsitatea (Ir, Ig, Ib)
    double F[3];    // Fokuaren bektorea (Fx, Fy, Fz)
    koordlist *dir; // Eguzkiaren direkzioaren bektorea (dirX, dirY, dirZ)
    koordlist *pos; // Bonbillaren posizioa (x, y, z)
    int fokua;      // fokua bada, 1, bestela 0
    double kd[3];   // isla barreiatua (difuse)
    double ks[3];   // ispilu isla (specular)
    double foku_irekiera;
    } argia;

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy,dimentsioa;

int indexx;
hiruki *triangulosptr;
triobj *foptr;
triobj *sel_ptr;
triobj *fCamptr;
triobj *selCam_ptr;
argia *fArgi_ptr;
argia *argiak_ptr;
argia *selArgi_ptr;
int denak;
int lineak;
int objektuak;
int kontrola;
int kameraOBJ;
int analisi;
char aldaketa;
int ald_lokala;
int perspektiba;
int back_culling;
int gorria;
double mesa[16];
double desplazamendua;
int argia_index;
double ikuste_bol_aldaketa;
double IaR, IaG, IaB;
double KaR, KaG, KaB;
double triangeluaren_I[3];
double f_d;

char fitxiz[100];


void ordenatu(punto *p1ptr, punto *p2ptr, punto *p3ptr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr);
void triangelua_bete(triobj *optr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr, float tartea);
void triangelua_bete_kasu_berezia(triobj *optr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr, float tartea);
void ordenatu_x(punto *p1ptr, punto *p2ptr, punto *p3ptr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr);
void mxm(double mBerria[16], double m1[16], double m2[16]);
void modelView_kalkulatu(triobj *optr, double modelView[16]);
void mP_paraleloa_kalkulatu(triobj *optr, double mP[16]);
void mP_perspektiba_kalkulatu(triobj *optr, double mP[16]);
int mxp_zati_w(punto *pptr, double m[16], punto p);
void mxN(double m[16], hiruki *tptr, double vBerria[3]);
void mxv(double m[16], double koordenatuak[3], double vBerria[3]);
void kamerak_hasieratu();
void analisi_bektoreak();
void biderketa_bektoriala(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double *x, double *y, double *z);
double biderketa_eskalarra(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z);
void N_kalkulatu(triobj *optr);
void print_egoerak();
void mESA_kalkulatu(double m[16]);
void analisi_biraketa(double vx, double vy, double vz, double alpha);
void mESA_eguneratu();
int talka();
void argiak_hasieratu();
void fokuak_kalkulatu(int hasieratu);
void piztu_itzali(int index);
void print_piztutako_argiak();
void foko_angelu_aldaketa(int handitu);
void undo_argiak();
double f_dist(double d);
void intentsitatea_kalkulatu(double N[3], double x, double y, double z, double *Ir, double *Ig, double *Ib);
double puntuen_arteko_distantzia(double x1, double y1, double z1, double x2, double y2, double z2);
void puntuen_arteko_bektorea(double x1, double y1, double z1, double x2, double y2, double z2, double *vx, double *vy, double *vz);
void H_halkulatu(double Vx, double Vy, double Vz, double Lx, double Ly, double Lz, double *Hx, double *Hy, double *Hz);
double I_formula(double I, double kd, double ks, double NxL, double NxH);
int foko_barruan_dago(double F[3], double Lx, double Ly, double Lz, double irekiera);
void mxv_2(double m[16], double vx, double vy, double vz, double *x, double *y, double *z);

void argiak_hasieratu() {
    int i;
    argia *argi_ptr;
    mlist *mlag;

    argiak_ptr = (argia *)malloc(4*sizeof(argia));

    // eguzkia
    argiak_ptr[0].fokua = 0;
    argiak_ptr[0].piztu = 1;
    
    argiak_ptr[0].dir = (koordlist *)malloc(sizeof(koordlist));
    argiak_ptr[0].dir->koord[0] = 0;
    argiak_ptr[0].dir->koord[1] = -1;
    argiak_ptr[0].dir->koord[2] = 0;
    argiak_ptr[0].dir->hptr = 0;

    argiak_ptr[0].I[0] = 1;
    argiak_ptr[0].I[1] = 1;
    argiak_ptr[0].I[2] = 1;

    argiak_ptr[0].ks[0] = 1;
    argiak_ptr[0].ks[1] = 0;
    argiak_ptr[0].ks[2] = 0;

    argiak_ptr[0].kd[0] = 0.25;
    argiak_ptr[0].kd[1] = 0.25;
    argiak_ptr[0].kd[2] = 0.25;


    // bonbila
    argiak_ptr[1].fokua = 0;
    argiak_ptr[1].piztu = 0;
    
    argiak_ptr[1].pos = (koordlist *)malloc(sizeof(koordlist));
    argiak_ptr[1].pos->koord[0] = -0.5;
    argiak_ptr[1].pos->koord[1] = 0.3;
    argiak_ptr[1].pos->koord[2] = 0.5;
    argiak_ptr[1].pos->hptr = 0;

    argiak_ptr[1].I[0] = 1;
    argiak_ptr[1].I[1] = 1;
    argiak_ptr[1].I[2] = 1;

    argiak_ptr[1].ks[0] = 0;
    argiak_ptr[1].ks[1] = 1;
    argiak_ptr[1].ks[2] = 0;

    argiak_ptr[1].kd[0] = 0.2;
    argiak_ptr[1].kd[1] = 0.2;
    argiak_ptr[1].kd[2] = 0.2;

    fokuak_kalkulatu(1);
    
}

void fokuak_kalkulatu(int hasieratu) {
    argia *argi_ptr;
    triobj *lag_ptr;
    int index;

    if (hasieratu == 1) {

        // objektuaren fokua
        argiak_ptr[2].fokua = 1;
        argiak_ptr[2].piztu = 0;

        (argiak_ptr[2]).F[0] = -(sel_ptr->mptr->m[2]);
        (argiak_ptr[2]).F[1] = -(sel_ptr->mptr->m[6]);
        (argiak_ptr[2]).F[2] = -(sel_ptr->mptr->m[10]);
        argiak_ptr[2].foku_irekiera = 0.5;

        argiak_ptr[2].I[0] = 1;
        argiak_ptr[2].I[1] = 1;
        argiak_ptr[2].I[2] = 1;

        argiak_ptr[2].ks[0] = 1;
        argiak_ptr[2].ks[1] = 1;
        argiak_ptr[2].ks[2] = 0;

        argiak_ptr[2].kd[0] = 0.2;
        argiak_ptr[2].kd[1] = 0.2;
        argiak_ptr[2].kd[2] = 0.2;

        // kameraren fokua
        argiak_ptr[3].fokua = 1;
        argiak_ptr[3].piztu = 0;

        (argiak_ptr[3]).F[0] = -(selCam_ptr->mptr->m[2]);
        (argiak_ptr[3]).F[1] = -(selCam_ptr->mptr->m[6]);
        (argiak_ptr[3]).F[2] = -(selCam_ptr->mptr->m[10]);
        (argiak_ptr[3]).foku_irekiera = 0.5;

        argiak_ptr[3].I[0] = 1;
        argiak_ptr[3].I[1] = 1;
        argiak_ptr[3].I[2] = 1;

        argiak_ptr[3].kd[0] = 0.25;
        argiak_ptr[3].kd[1] = 0.25;
        argiak_ptr[3].kd[2] = 0.25;

    } else if (hasieratu == 0) { // kameraren edo objektuaren fokuak eguneratu
        if (kontrola == 0) {
            lag_ptr = sel_ptr;
            index = 2;
        } else if (kontrola == 1) {
            lag_ptr = selCam_ptr;
            index = 3;
        }
        argiak_ptr[index].F[0] = -(lag_ptr->mptr->m[2]);
        argiak_ptr[index].F[1] = -(lag_ptr->mptr->m[6]);
        argiak_ptr[index].F[2] = -(lag_ptr->mptr->m[10]);
        //printf("FOKUA: x: %f, y: %f, z:%f\n", argiak_ptr[index].F[0], argiak_ptr[index].F[1], argiak_ptr[index].F[2]);
    } 

}

void foko_angelu_aldaketa(int handitu) {
    double alpha, irekieraBerria;
    if (handitu == 1) alpha = 0.1;
    else alpha = -0.1;

    irekieraBerria = selArgi_ptr->foku_irekiera + alpha;
    if (irekieraBerria < 0 || irekieraBerria > PI) return;

    selArgi_ptr->foku_irekiera = irekieraBerria;
    //printf("Foku irekiera: %f\n", selArgi_ptr->foku_irekiera);
}

double f_dist(double d)
{
    double a, b, c, izendatzailea;
    a = 0.05; b = 0.4; c = 0.5;

    izendatzailea = (a*pow(d, 2.0) + b*d + c);
    if (izendatzailea < 1 || izendatzailea == 0) return 1;
    else return (1 / izendatzailea);
}

void intentsitatea_kalkulatu(double N[3], double x, double y, double z, double *IR, double *IG, double *IB) {
    double Lx, Ly, Lz, Hx, Hy, Hz, Vx, Vy, Vz, NxL, NxH, intenR, intenG, intenB, dist, foko_kam[3], L_norm;
    triobj *aux_ptr;
    int foko_barruan;

    intenR = IaR * KaR;
    intenG = IaG * KaG;
    intenB = IaB * KaB;

    if (kameraOBJ == 1) {
        aux_ptr = sel_ptr;
    } else if (kameraOBJ == 0) {
        aux_ptr = selCam_ptr;
    }

    // V
    //puntuen_arteko_bektorea(x, y, z, aux_ptr->mptr->m[3], aux_ptr->mptr->m[7], aux_ptr->mptr->m[11], &Vx, &Vy, &Vz);
    //mxv_2(mesa, Vx, Vy, Vz, &Vx, &Vy, &Vz);
    Vx = -x;
    Vy = -y;
    Vz = -z;

    // EGUZKIA
    f_d = 1;
    // L
    Lx = -(argiak_ptr[0].dir->koord[0]);
    Ly = -(argiak_ptr[0].dir->koord[1]);
    Lz = -(argiak_ptr[0].dir->koord[2]);
    mxv_2(mesa, Lx, Ly, Lz, &Lx, &Ly, &Lz);
    L_norm = sqrt(pow(Lx, 2) + pow(Ly, 2) + pow(Lz, 2));
    Lx = Lx / L_norm;
    Ly = Ly / L_norm;
    Lz = Lz / L_norm;
    // H
    H_halkulatu(Vx, Vy, Vz, Lx, Ly, Lz, &Hx, &Hy, &Hz);
    //mxv_2(mesa, Hx, Hy, Hz, &Hx, &Hy, &Hz);

    NxL = biderketa_eskalarra(N[0], N[1], N[2], Lx, Ly, Lz);
    if (NxL < 0) NxL = 0;
    NxH = biderketa_eskalarra(N[0], N[1], N[2], Hx, Hy, Hz);
    if (NxH < 0) NxH = 0;

    if (NxL >= 1 || NxH >= 1) printf("NxL: %f, NxH: %f\n", NxL, NxH);

    intenR += (argiak_ptr[0].piztu) * I_formula(argiak_ptr[0].I[0], argiak_ptr[0].kd[0], argiak_ptr[0].ks[0], NxL, NxH);
    intenG += (argiak_ptr[0].piztu) * I_formula(argiak_ptr[0].I[1], argiak_ptr[0].kd[1], argiak_ptr[0].ks[1], NxL, NxH);
    intenB += (argiak_ptr[0].piztu) * I_formula(argiak_ptr[0].I[2], argiak_ptr[0].kd[2], argiak_ptr[0].ks[2], NxL, NxH);

    // BONBILA
    dist = puntuen_arteko_distantzia(x, y, z, argiak_ptr[1].pos->koord[0], argiak_ptr[1].pos->koord[1], argiak_ptr[1].pos->koord[2]);
    f_d = f_dist(dist);
    //f_d = 1;
    // L
    puntuen_arteko_bektorea(x, y, z, argiak_ptr[1].pos->koord[0], argiak_ptr[1].pos->koord[1], argiak_ptr[1].pos->koord[2], &Lx, &Ly, &Lz);
    mxv_2(mesa, Lx, Ly, Lz, &Lx, &Ly, &Lz);
    L_norm = sqrt(pow(Lx, 2) + pow(Ly, 2) + pow(Lz, 2));
    Lx = Lx / L_norm;
    Ly = Ly / L_norm;
    Lz = Lz / L_norm;

    // H
    H_halkulatu(Vx, Vy, Vz, Lx, Ly, Lz, &Hx, &Hy, &Hz);
    //mxv_2(mesa, Hx, Hy, Hz, &Hx, &Hy, &Hz);

    NxL = biderketa_eskalarra(N[0], N[1], N[2], Lx, Ly, Lz);
    if (NxL < 0) NxL = 0;
    NxH = biderketa_eskalarra(N[0], N[1], N[2], Hx, Hy, Hz);
    if (NxH < 0) NxH = 0;
    

    intenR += (argiak_ptr[1].piztu) * I_formula(argiak_ptr[1].I[0], argiak_ptr[1].kd[0], argiak_ptr[1].ks[0], NxL, NxH);
    intenG += (argiak_ptr[1].piztu) * I_formula(argiak_ptr[1].I[1], argiak_ptr[1].kd[1], argiak_ptr[1].ks[1], NxL, NxH);
    intenB += (argiak_ptr[1].piztu) * I_formula(argiak_ptr[1].I[2], argiak_ptr[1].kd[2], argiak_ptr[1].ks[2], NxL, NxH);

    // OBJEKTUA
    dist = puntuen_arteko_distantzia(x, y, z, sel_ptr->mptr->m[3], sel_ptr->mptr->m[7], sel_ptr->mptr->m[11]);
    f_d = f_dist(dist);
    //f_d = 1;
    
    // L
    if (kameraOBJ == 0) {
        puntuen_arteko_bektorea(x, y, z, sel_ptr->mptr->m[3], sel_ptr->mptr->m[7], sel_ptr->mptr->m[11], &Lx, &Ly, &Lz);
        mxv_2(mesa, Lx, Ly, Lz, &Lx, &Ly, &Lz);
    } else {
        Lx = Vx;
        Ly = Vy;
        Lz = Vz;
    }
    L_norm = sqrt(pow(Lx, 2) + pow(Ly, 2) + pow(Lz, 2));
    Lx = Lx / L_norm;
    Ly = Ly / L_norm;
    Lz = Lz / L_norm;
    // H
    H_halkulatu(Vx, Vy, Vz, Lx, Ly, Lz, &Hx, &Hy, &Hz);
    //mxv_2(mesa, Hx, Hy, Hz, &Hx, &Hy, &Hz);

    NxL = biderketa_eskalarra(N[0], N[1], N[2], Lx, Ly, Lz);
    if (NxL < 0) NxL = 0;
    NxH = biderketa_eskalarra(N[0], N[1], N[2], Hx, Hy, Hz);
    if (NxH < 0) NxH = 0;

    foko_barruan = foko_barruan_dago(argiak_ptr[2].F, Lx, Ly, Lz, argiak_ptr[2].foku_irekiera);

    intenR += (argiak_ptr[2].piztu) * foko_barruan * I_formula(argiak_ptr[2].I[0], argiak_ptr[2].kd[0], argiak_ptr[2].ks[0], NxL, NxH);
    intenG += (argiak_ptr[2].piztu) * foko_barruan * I_formula(argiak_ptr[2].I[1], argiak_ptr[2].kd[1], argiak_ptr[2].ks[1], NxL, NxH);
    intenB += (argiak_ptr[2].piztu) * foko_barruan * I_formula(argiak_ptr[2].I[2], argiak_ptr[2].kd[2], argiak_ptr[2].ks[2], NxL, NxH);

    // KAMERA
    dist = puntuen_arteko_distantzia(x, y, z, selCam_ptr->mptr->m[3], selCam_ptr->mptr->m[7], selCam_ptr->mptr->m[11]);
    f_d = f_dist(dist);
    //f_d = 1;
    // L
    //puntuen_arteko_bektorea(x, y, z, selCam_ptr->mptr->m[3], selCam_ptr->mptr->m[7], selCam_ptr->mptr->m[11], &Lx, &Ly, &Lz);
    //mxv_2(mesa, Lx, Ly, Lz, &Lx, &Ly, &Lz);
    Lx = Vx;
    Ly = Vy;
    Lz = Vz;
    L_norm = sqrt(pow(Lx, 2) + pow(Ly, 2) + pow(Lz, 2));
    Lx = Lx / L_norm;
    Ly = Ly / L_norm;
    Lz = Lz / L_norm;
    // H
    H_halkulatu(Vx, Vy, Vz, Lx, Ly, Lz, &Hx, &Hy, &Hz);
    //mxv_2(mesa, Hx, Hy, Hz, &Hx, &Hy, &Hz);

    NxL = biderketa_eskalarra(N[0], N[1], N[2], Lx, Ly, Lz);
    if (NxL < 0) NxL = 0;
    NxH = biderketa_eskalarra(N[0], N[1], N[2], Hx, Hy, Hz);
    if (NxH < 0) NxH = 0;

    foko_kam[0] = 0; foko_kam[1] = 0; foko_kam[2] = -1;
    foko_barruan = foko_barruan_dago(foko_kam, Lx, Ly, Lz, argiak_ptr[3].foku_irekiera);
    //foko_barruan = foko_barruan_dago(argiak_ptr[3].F, Lx, Ly, Lz, argiak_ptr[3].foku_irekiera);

    intenR += (argiak_ptr[3].piztu) * foko_barruan * I_formula(argiak_ptr[3].I[0], argiak_ptr[3].kd[0], argiak_ptr[3].ks[0], NxL, NxH);
    intenG += (argiak_ptr[3].piztu) * foko_barruan * I_formula(argiak_ptr[3].I[1], argiak_ptr[3].kd[1], argiak_ptr[3].ks[1], NxL, NxH);
    intenB += (argiak_ptr[3].piztu) * foko_barruan * I_formula(argiak_ptr[3].I[2], argiak_ptr[3].kd[2], argiak_ptr[3].ks[2], NxL, NxH);

    //printf("Ir: %f, Ig: %f, Ib: %f\n", intenR, intenG, intenB);

    triangeluaren_I[0] = intenR;
    triangeluaren_I[1] = intenG;
    triangeluaren_I[2] = intenB;
    if (triangeluaren_I[0] > 1) triangeluaren_I[0] = 1;
    if (triangeluaren_I[1] > 1) triangeluaren_I[1] = 1;
    if (triangeluaren_I[2] > 1) triangeluaren_I[2] = 1;

}

double I_formula(double I, double kd, double ks, double NxL, double NxH) {
    return (I*(kd*NxL + ks*pow(NxH, 4*50)));
}

int foko_barruan_dago(double F[3], double Lx, double Ly, double Lz, double irekiera) {
    double F_x_minusL, fx, fy, fz;

    if (F[0] != 0 || F[1] != 0 || F[2] != -1) mxv_2(mesa, F[0], F[1], F[2], &fx, &fy, &fz);
    else {fx = F[0]; fy = F[1]; fz = F[2];}

    F_x_minusL = biderketa_eskalarra(fx, fy, fz, -Lx, -Ly, -Lz);
    if (F_x_minusL > cos(irekiera)) {
        return 1;
    } else {
        return 0;
    }
}

void H_halkulatu(double Vx, double Vy, double Vz, double Lx, double Ly, double Lz, double *Hx, double *Hy, double *Hz) {
    double VgehiL[3], VgehiL_norm;

    VgehiL[0] = Vx + Lx;    VgehiL[1] = Vy + Ly;    VgehiL[2] = Vz + Lz;
    VgehiL_norm = sqrt(pow(VgehiL[0], 2.0) + pow(VgehiL[1], 2.0) + pow(VgehiL[2], 2.0));
    *Hx = VgehiL[0] / VgehiL_norm;
    *Hy = VgehiL[1] / VgehiL_norm;
    *Hz = VgehiL[2] / VgehiL_norm;
}

void puntuen_arteko_bektorea(double x1, double y1, double z1, double x2, double y2, double z2, double *vx, double *vy, double *vz) {
    *vx = x2 - x1;
    *vy = y2 - y1;
    *vz = z2 - z1;
}

double puntuen_arteko_distantzia(double x1, double y1, double z1, double x2, double y2, double z2) {
    double vx, vy, vz, vnorm;

    puntuen_arteko_bektorea(x1, y1, z1, x2, y2, z2, &vx, &vy, &vz);
    vnorm = sqrt(pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0));
    return vnorm;
}


void objektuari_aldaketa_sartu_ezk(double m[16])
{
    int j;
    mlist *mlag;
    triobj *aux_ptr;
    koordlist *koordlag;
    double mBerria[16];
    double koordBerriak[3];

    if (kontrola == 1) {
        aux_ptr = selCam_ptr;
    } else if(kontrola == 0) {
        aux_ptr = sel_ptr;
    } else if (kontrola == 2) {

        koordlag = (koordlist *)malloc(sizeof(koordlist));

        if (argia_index == 0) {
            mxv(m, selArgi_ptr->dir->koord, koordBerriak);
            koordlag->hptr = selArgi_ptr->dir;
            selArgi_ptr->dir = koordlag;
        } else if (argia_index == 1) {
            mxv(m, selArgi_ptr->pos->koord, koordBerriak);
            koordlag->hptr = selArgi_ptr->pos;
            selArgi_ptr->pos = koordlag;
        }
        for (j=0; j<3; j++){
            koordlag->koord[j]=koordBerriak[j];
        }
        //printf("x: %f, y: %f, z: %f\n", koordlag->koord[0], koordlag->koord[1], koordlag->koord[2]);
        return;
        
    }

    mlag = (mlist *)malloc(sizeof(mlist));

    mxm(mBerria, m, aux_ptr->mptr->m);

    for (j=0; j<16; j++){
        mlag->m[j]=mBerria[j];
    }
    mlag->hptr=aux_ptr->mptr;
    aux_ptr->mptr=mlag;
    fokuak_kalkulatu(0);

}

void objektuari_aldaketa_sartu_esk(double m[16])
{
    int j;
    mlist *mlag;
    triobj *aux_ptr;
    double mBerria[16];

    if (kontrola == 1) {
        aux_ptr = selCam_ptr;
    } else if(kontrola == 0) {
        aux_ptr = sel_ptr;
    }
    
    mlag = (mlist *)malloc(sizeof(mlist));
    
    mxm(mBerria, aux_ptr->mptr->m, m);

    for (j=0; j<16; j++){
        mlag->m[j]=mBerria[j];
    }
    mlag->hptr=aux_ptr->mptr;
    aux_ptr->mptr=mlag;
    fokuak_kalkulatu(0);
    
}



// TODO
// funtzio honek u eta v koordenatuei dagokien pointerra itzuli behar du.
// debe devolver el pointer correspondiente a las coordenadas u y v
unsigned char * color_textura(float u, float v)
{
int dx, dy;
char * lag;
//printf("texturan...%x\n",bufferra);

dx = u * (dimx-1);

dy = (1-v)*(dimy-1);

int desp = dx + dimx*dy;

lag = (unsigned char *)bufferra;

if (u < 0 || u > 1 || v < 0 || v > 1) {
    return(lag);
} else {
    return(lag+3*desp);
}


}


// TODO
// lerroa marrazten du, baina testuraren kodea egokitu behar da
// dimentsioa --> marrazgunearen pixel kopurua
// dibuja una linea pero hay que codificar la textura
// dimentsioa --> número de pixels de la ventana de dibujo
void  dibujar_linea_z(triobj *optr, float linea,float c1x,float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v)
{
float xkoord,zkoord;
float cx, cz, cu, cv;
float u,v;
int i,puntukop;
unsigned char r,g,b;
unsigned char *colorv;

// TODO x balioak -1 eta 1 artekoak direla ziurtatu eta ondorioz z, u eta v egokitu

 // x balioak -1 eta 1 artekoa izan behar du, ezkerretik eskuinera 2 unitateko aldaketa dauka x balioak, eta 
 // dimentsioa adina pixel behar dira ezkerretik eskuinera joateko, beraz, 2-ko aldaketari "dimentsioa" pixel dagozkio.
 // ondorioz c1x-tik c2x-ra doan tartean behar ditudan pixelak = (c2x-c1x)dimentsioa/2 pixel behar ditut
 // ondorioz, Z-ren u-ren eta v-ren aldaketan zenbaki hori erabili behar dut.
 // para un cambio de -1 a 1 en x hay que dibujar "dimentsioa" pixels, para un cambio que va de c1x a c2x, cuantos?
 puntukop = (c2x-c1x)*(float)dimentsioa;

if (puntukop > 0) 
    {
    // pixel batetik bestera dagoen distantzia x koordenatuan: cx = (c2x-c1x)/(float)puntukop;
    // diustancia en x de un pixel al siguiente: cx
    cx = (c2x-c1x)/(float)puntukop;
    cz = (c2z-c1z)/(float)puntukop;
    cu = (c2u-c1u)/(float)puntukop;
    cv = (c2v-c1v)/(float)puntukop;
    }
  else 
    {
    cx = 0;
    cz = 0;
    cu = 0;
    cv = 0;
    }
 

glBegin( GL_POINTS );
for (i = -1, xkoord=c1x, zkoord=c1z, u=c1u, v=c1v; 
     i <=puntukop /*xkoord <= c2x*/; 
     i++,xkoord += cx /* 500 puntu -1 eta 1 artean */)
    {
    if (gorria == 1) {
        r = 255.0;
        g = 0.0;
        b = 0.0;
    } else if (gorria == 0) {
        if (optr->rgbptr == 0){
        colorv=  color_textura(u, v); // si esta función es correcta se ve la foto en la ventana
        r= colorv[0];
        g=colorv[1];
        b=colorv[2];    
        } else {
            r = optr->rgbptr[0];
            g = optr->rgbptr[1];
            b = optr->rgbptr[2];
        }
        r += triangeluaren_I[0]*255;
        g += triangeluaren_I[1]*255;
        b += triangeluaren_I[2]*255;
        if (r > 255) r = 255;
        if (g > 255) g = 255;
        if (b > 255) b = 255;
    }
    
    glColor3ub(r,g,b);
    glVertex3f(xkoord, linea, zkoord );
    zkoord+=cz;
    u+=cu;
    v+=cv;
    }
glEnd();
}


void print_matrizea(char *str)
{
int i;

printf("%s\n",str);
for (i = 0;i<4;i++)
   printf("%lf, %lf, %lf, %lf\n",sel_ptr->mptr->m[i*4],sel_ptr->mptr->m[i*4+1],sel_ptr->mptr->m[i*4+2],
                                 sel_ptr->mptr->m[i*4+3]);
}


// TODO (transform)
// supone que la cuarta coordenada es un 1
void mxp(punto *pptr, double m[16], punto p)
{
//print_matrizea("objektuaren matrizea\n");
//printf("puntua = %lf,%lf,%lf\n ",p.x,p.y,p.z);
pptr->x = m[0]*p.x + m[1]*p.y + m[2]*p.z + m[3];
pptr->y = m[4]*p.x + m[5]*p.y + m[6]*p.z + m[7];
pptr->z = m[8]*p.x + m[9]*p.y + m[10]*p.z + m[11];
pptr->u = p.u;
pptr->v = p.v;

}

int mxp_zati_w(punto *pptr, double m[16], punto p)
{
double w;

w = m[12]*p.x + m[13]*p.y + m[14]*p.z + m[15];
if (w==0.0) return 0;
pptr->x = (m[0]*p.x + m[1]*p.y + m[2]*p.z + m[3])/w;
pptr->y = (m[4]*p.x + m[5]*p.y + m[6]*p.z + m[7])/w;
pptr->z = (m[8]*p.x + m[9]*p.y + m[10]*p.z + m[11])/w;
pptr->u = p.u;
pptr->v = p.v;
return 1;

}

void mxN(double m[16], hiruki *tptr, double vBerria[3])
{
double vx, vy , vz, vx2, vy2, vz2;

vx = tptr->N[0];
vy = tptr->N[1];
vz = tptr->N[2];

// N-ren azken koordenatua 0: (vx, vy, vz, 0)
vx2 = m[0]*vx + m[1]*vy + m[2]*vz;
vy2 = m[4]*vx + m[5]*vy + m[6]*vz;
vz2 = m[8]*vx + m[9]*vy + m[10]*vz;

vBerria[0] = vx2;
vBerria[1] = vy2;
vBerria[2] = vz2;
}

void mxv(double m[16], double koordenatuak[3], double vBerria[3]) {
    double vx, vy , vz, vx2, vy2, vz2;

    vx = koordenatuak[0];
    vy = koordenatuak[1];
    vz = koordenatuak[2];

    vx2 = m[0]*vx + m[1]*vy + m[2]*vz + m[3];
    vy2 = m[4]*vx + m[5]*vy + m[6]*vz + m[7];
    vz2 = m[8]*vx + m[9]*vy + m[10]*vz + m[11];

    vBerria[0] = vx2;
    vBerria[1] = vy2;
    vBerria[2] = vz2;

}

void mxv_2(double m[16], double vx, double vy, double vz, double *x, double *y, double *z) {
    
    *x = m[0]*vx + m[1]*vy + m[2]*vz + m[3];
    *y = m[4]*vx + m[5]*vy + m[6]*vz + m[7];
    *z = m[8]*vx + m[9]*vy + m[10]*vz + m[11];

}

// bi matrize biderkatzen ditu eta mBerria matrize berrian gordetzen du emaitza
void mxm(double mBerria[16], double m1[16], double m2[16]) {
    int i, pos;

    for (i=0; i<16; i++) {
        mBerria[i] = 0;
        for (pos = 0; pos<4; pos++) {
            mBerria[i] += (m1[((i/4)*4)+pos] * m2[(i % 4)+(4*pos)]);
        }   
    }
}

void modelView_kalkulatu(triobj *optr, double modelView[16]) {
    
    mxm(modelView, mesa, optr->mptr->m);
}

void analisi_bektoreak() {
    double zx, zy, zz, xx, xy, xz, yx, yy, yz, znorm, xnorm;
    
    // zk = (E - at) / ||E - at||
    zx = selCam_ptr->mptr->m[3] - sel_ptr->mptr->m[3];
    zy = selCam_ptr->mptr->m[7] - sel_ptr->mptr->m[7];
    zz = selCam_ptr->mptr->m[11] - sel_ptr->mptr->m[11];
    znorm = sqrt(pow(zx, 2.0)+pow(zy, 2.0)+pow(zz, 2.0));
    zx = zx/znorm;
    zy = zy/znorm;
    zz = zz/znorm;
   
    // xk = (Vup ^ zk) / ||Vup ^ zk||, non Vup=kameraren hasierako y bektorea
    biderketa_bektoriala(selCam_ptr->mptr->m[1], selCam_ptr->mptr->m[5], selCam_ptr->mptr->m[9], zx, zy, zz, &xx, &xy, &xz);
    xnorm = sqrt(pow(xx, 2.0)+pow(xy, 2.0)+pow(xz, 2.0));
    xx = xx/xnorm;
    xy = xy/xnorm;
    xz = xz/xnorm;

    // yk = zk ^ xk
    biderketa_bektoriala(zx, zy, zz, xx, xy, xz, &yx, &yy, &yz);
    
    selCam_ptr->mptr->m[0]=xx; selCam_ptr->mptr->m[4]=xy; selCam_ptr->mptr->m[8]=xz;
    selCam_ptr->mptr->m[1]=yx; selCam_ptr->mptr->m[5]=yy; selCam_ptr->mptr->m[9]=yz;
    selCam_ptr->mptr->m[2]=zx; selCam_ptr->mptr->m[6]=zy; selCam_ptr->mptr->m[10]=zz;

}

double biderketa_eskalarra(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z) {
    return (v1x*v2x + v1y*v2y + v1z*v2z);
}

void biderketa_bektoriala(double v1x, double v1y, double v1z, double v2x, double v2y, double v2z, double *x, double *y, double *z) {
    *x = v1y*v2z - v1z*v2y;
    *y = -(v1x*v2z - v1z*v2x);
    *z = v1x*v2y - v1y*v2x;
    if (*x == -0.0) *x=0.0;
    if (*y == -0.0) *y=0.0;
    if (*z == -0.0) *z=0.0;
}

void mP_paraleloa_kalkulatu(triobj *optr, double mP[16]) {
    int i;
    double /*m1[16], m2[16],*/ right, left, top, bottom, near, far;
    mlist *mObj;

    if (kameraOBJ == 1 && foptr!=0) {
        mObj = sel_ptr->mptr;
    } else {
        mObj = selCam_ptr->mptr;
    }   
    
    right = 1 + ikuste_bol_aldaketa*10;
    left = -right;
    top = 1 + ikuste_bol_aldaketa*10;
    bottom = -top;
    near = 0;
    far = -100; 

    for (i=0; i<16; i++) {
        mP[i] = 0;
    }
   
    mP[0] = 2.0 / (right-left);
    mP[5] = 2.0 / (top-bottom);
    mP[10] = 2.0 / (near-far);
    mP[15] = 1;
    mP[3] = -((right+left) / (right-left));
    mP[7] = -((top+bottom) / (top-bottom));
    mP[11] = -((near+far) / (near-far));

}

void mP_perspektiba_kalkulatu(triobj *optr, double mP[16]) {
    int i;
    double right, left, top, bottom, near, far;
    mlist *mObj;

    if (kameraOBJ == 1 && foptr!=0) {
        mObj = sel_ptr->mptr;
    } else {
        mObj = selCam_ptr->mptr;
    }   
    right = 0.1 + ikuste_bol_aldaketa;
    left = -right;
    top = 0.1 + ikuste_bol_aldaketa;
    bottom = -top;
    near = 0.1;
    far = 100.0;

    for (i=0; i<16; i++) {
        mP[i] = 0;
    }

    mP[0] = (2.0*near) / (right-left);
    mP[5] = (2.0*near) / (top-bottom);
    mP[2] = (right+left) / (right-left);
    mP[6] = (top+bottom) / (top-bottom);
    mP[10] = (-(near+far) / (near-far));
    mP[11] = ((-2.0*near*far) / (near-far));
    mP[14] = -1.0;

}


// TODO (transform)
// objektua munduan kokatzen duen matrizetik abiatuta objktuaren erreferentzi sistemara pasatzen duen matrizea lortu
void obtener_CSR_partiendo_de_M(GLdouble* M, GLdouble* MCSR) {
    MCSR[0] = M[0]; MCSR[4] = M[1]; MCSR[8] = M[2]; MCSR [12] = 0;
    MCSR[1] = M[4]; MCSR[5] = M[5]; MCSR[9] = M[6]; MCSR [13] = 0;
    MCSR[2] = M[8]; MCSR[6] = M[9]; MCSR[10] = M[10]; MCSR [14] = 0;
    MCSR[3] = 0;    MCSR[7] = 0;    MCSR[11] = 0;     MCSR [15] = 1;
}

void mESA_kalkulatu(double m[16]) {
    int i, j;

    for (i=0; i<4; i++) { // objektuaren matrizearen iraulia xk, yk, eta zk lortzeko
        for (j=0; j<4; j++) {
            mesa[4*j + i] = m[4*i + j];
        }
    }

    // -E.x, -E.y, -E.z
    mesa[3] = -(biderketa_eskalarra(m[3], m[7], m[11], m[0], m[4], m[8]));
    mesa[7] = -(biderketa_eskalarra(m[3], m[7], m[11], m[1], m[5], m[9]));
    mesa[11] = -(biderketa_eskalarra(m[3], m[7], m[11], m[2], m[6], m[10]));

    // azken errenkada
    mesa[12]=0; mesa[13]=0; mesa[14]=0;
    mesa[15] = 1;
}

void mESA_eguneratu() {
    if (foptr!=0 && kameraOBJ==1) {
        mESA_kalkulatu(sel_ptr->mptr->m);
    } else {      
        mESA_kalkulatu(selCam_ptr->mptr->m);
    }
}

void dibujar_triangulo(triobj *optr, int ti)
{
hiruki *tptr;

punto *pgoiptr, *pbeheptr, *perdiptr;
float x1,y1,z1,u1,v1,x2,y2,z2,u2,v2,x3,y3,z3,u3,v3, p1x, p1y, p1z;
float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v;
float y; // \in  (-1,1)
float lerrotartea,cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
int lerrokop, marraztu, retval1, retval2, retval3;
punto p1,p2,p3;
double modelView[16], mP[16], mPxmodelView[16], nBerria[3], vExnBerria, Ir, Ig, Ib;

gorria = 0;
marraztu = 1;
//printf("hirukian\n");
if (ti >= optr->num_triangles) return;
tptr = optr->triptr+ti;

modelView_kalkulatu(optr, modelView);

mxN(modelView, tptr, nBerria); // N kameraren erreferentzia sisteman jarri
mxp(&p1,modelView,tptr->p1);
mxp(&p2,modelView,tptr->p2);
mxp(&p3,modelView,tptr->p3);

intentsitatea_kalkulatu(nBerria, p1.x, p1.y, p1.z, &Ir, &Ig, &Ib);

if (perspektiba == 0) {
    mP_paraleloa_kalkulatu(optr, mP);
} else {
    if (p1.z>-0.1 || p2.z>-0.1 || p3.z>-0.1) return;
    mP_perspektiba_kalkulatu(optr, mP);
}

mxm(mPxmodelView, mP, modelView);

retval1=mxp_zati_w(&p1,mPxmodelView,tptr->p1);
retval2=mxp_zati_w(&p2,mPxmodelView,tptr->p2);
retval3=mxp_zati_w(&p3,mPxmodelView,tptr->p3);
if (retval1==0 || retval2==0 || retval3==0) return;


// triangeluaren bektore normala marrazteko
/*
glBegin(GL_LINES);
    glVertex3d(p1.x, p1.y, p1.z);
    glVertex3d(p1.x+nBerria[0], p1.y + nBerria[1], p1.z + nBerria[2]);
glEnd();
*/

if (perspektiba == 0 && nBerria[2] < 0) {
    marraztu = 0;
} else if (perspektiba == 1) {
    p1x = -(p1.x);
    p1y = -(p1.y);
    p1z = -(p1.z);
    vExnBerria = biderketa_eskalarra(p1x, p1y, p1z, nBerria[0], nBerria[1], nBerria[2]);
    if (vExnBerria < 0) marraztu = 0;
}


if (marraztu == 0) {
    gorria = 1;
    if (lineak == 1 && back_culling == 0){
        glBegin(GL_POLYGON);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3d(p1.x, p1.y, p1.z);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3d(p2.x, p2.y, p2.z);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    } else if (back_culling == 1) return;
    
}


if (lineak == 1 && marraztu == 1)
        {
        glBegin(GL_POLYGON);
        glColor3f(1.0f, 1.0f, 1.0f);
        glVertex3d(p1.x, p1.y, p1.z);
        glColor3f(1.0f, 1.0f, 1.0f);
        glVertex3d(p2.x, p2.y, p2.z);
        glColor3f(1.0f, 1.0f, 1.0f);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
        }
        //  else 



 //   lehenengo hiru erpinak ordenatu behar ditut
 //   ordenar los puntos
 
 ordenatu(&p1, &p2, &p3, &pgoiptr, &pbeheptr, &perdiptr);

// marrazgunearen neurria da "dimentsioa", alegia zenbat pixel dauden zabalean edota garaieran (txikiena)
// lerrotik lerrorako distantzia:
// distancia entre dos lineas horizontales:
lerrotartea = 2.0/(float)dimentsioa;
 
 
 //     kasu bereziak
 //     casos especiales
 if (pgoiptr->x == pbeheptr->x && pgoiptr->x == perdiptr->x) {
    triangelua_bete_kasu_berezia(optr, &pgoiptr, &pbeheptr, &perdiptr, lerrotartea); // x guztiak berdinak direnez, lerro bertikal bat marraztu, goiko puntutuk behekora arte
 } else if (pgoiptr->y == pbeheptr->y && pgoiptr->y == perdiptr->y) {
    ordenatu_x(&p1, &p2, &p3, &pgoiptr, &pbeheptr, &perdiptr); // x koordenatuak txikienetik handienera ordenatu
    dibujar_linea_z(optr, pgoiptr->y,pbeheptr->x,pbeheptr->z,pbeheptr->u,pbeheptr->v,pgoiptr->x,pgoiptr->z,pgoiptr->u,pgoiptr->v); // lerro horizontal bakar bat marraztu
} else {
    //     kasu orokorra: bi triangelu
    //     caso general: dos triangulos
    triangelua_bete(optr, &pgoiptr, &pbeheptr, &perdiptr, lerrotartea);
}

}

void ordenatu(punto *p1ptr, punto *p2ptr, punto *p3ptr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr) {
    
    if (p1ptr->y > p2ptr->y) {
        *pgoiptrptr=p1ptr;
        *pbeheptrptr=p2ptr;
    } else {
        *pgoiptrptr=p2ptr;
        *pbeheptrptr=p1ptr;
    }
    
    if (p3ptr->y > (*pgoiptrptr)->y ){
        *perdiptrptr=*pgoiptrptr;
        *pgoiptrptr=p3ptr;
    } else {
        if ((*pbeheptrptr)->y > p3ptr->y ) {
            *perdiptrptr=*pbeheptrptr;
            *pbeheptrptr=p3ptr;
        } else {
            *perdiptrptr = p3ptr;
        }
    }
   
    
}

void triangelua_bete(triobj *optr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr, float tartea) {
   
    float goi, behe, erdi;
    goi = (*pgoiptrptr)->y;
    erdi = (*perdiptrptr)->y;
    behe = (*pbeheptrptr)->y;
    float p1x, p1z, p1u, p1v, p2x, p2z, p2u, p2v;
    float yFinal1, yFinal2;
    float y;
    unsigned char r,g,b;
    unsigned char *colorv;

    for (y = goi; y>=behe; y-=tartea) {
        //printf("y: %f\n", y);

        // PUNTUEN ARTEKO ZUZENAREN FUNTZIOA: (x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
        yFinal1 = (y - (*pgoiptrptr)->y) / ((*pbeheptrptr)->y - (*pgoiptrptr)->y);  // (y-y1)/(y2-y1)
        p1x = yFinal1*((*pbeheptrptr)->x - (*pgoiptrptr)->x) + (*pgoiptrptr)->x;    // [(y-y1)/(y2-y1)]*(x2-x1) + x1
        p1z = yFinal1*((*pbeheptrptr)->z - (*pgoiptrptr)->z) + (*pgoiptrptr)->z;    // [(y-y1)/(y2-y1)]*(z2-z1) + z1
        p1u = yFinal1*((*pbeheptrptr)->u - (*pgoiptrptr)->u) + (*pgoiptrptr)->u;    // [(y-y1)/(y2-y1)]*(u2-u1) + u1
        p1v = yFinal1*((*pbeheptrptr)->v - (*pgoiptrptr)->v) + (*pgoiptrptr)->v;    // [(y-y1)/(y2-y1)]*(v2-v1) + v1

        if (y < erdi){
            //printf("y: %f\n", y);
            
            yFinal2 = (y - (*perdiptrptr)->y) / ((*pbeheptrptr)->y - (*perdiptrptr)->y);
            p2x = yFinal2*((*pbeheptrptr)->x - (*perdiptrptr)->x) + (*perdiptrptr)->x;
            p2z = yFinal2*((*pbeheptrptr)->z - (*perdiptrptr)->z) + (*perdiptrptr)->z;
            p2u = yFinal2*((*pbeheptrptr)->u - (*perdiptrptr)->u) + (*perdiptrptr)->u;
            p2v = yFinal2*((*pbeheptrptr)->v - (*perdiptrptr)->v) + (*perdiptrptr)->v;

            //printf("BEHE: x1: %f, z1: %f, x2: %f, z2:%f\n", p1x, p1z, p2x, p2z);
        } else {
            yFinal2 = (y - (*pgoiptrptr)->y) / ((*perdiptrptr)->y - (*pgoiptrptr)->y);
            p2x = yFinal2*((*perdiptrptr)->x - (*pgoiptrptr)->x) + (*pgoiptrptr)->x;
            p2z = yFinal2*((*perdiptrptr)->z - (*pgoiptrptr)->z) + (*pgoiptrptr)->z;
            p2u = yFinal2*((*perdiptrptr)->u - (*pgoiptrptr)->u) + (*pgoiptrptr)->u;
            p2v = yFinal2*((*perdiptrptr)->v - (*pgoiptrptr)->v) + (*pgoiptrptr)->v;
            //printf("ERDI: x1: %f, z1: %f, x2: %f, z2:%f\n", p1x, p1z, p2x, p2z);
        }
        

        if (p1x < p2x){
            dibujar_linea_z(optr, y,p1x,p1z,p1u,p1v,p2x,p2z,p2u,p2v);
        } else if (p1x > p2x) {
            dibujar_linea_z(optr, y,p2x,p2z,p2u,p2v,p1x,p1z,p1u,p1v);
        } else if (p1x == p2x){ // hiru puntuak zeharkako lerro bat osatzen dutenean puntuak zuzenean marrazteko
            glBegin( GL_POINTS );

            if (optr->rgbptr == 0){
                colorv=  color_textura(p1u, p1v); // si esta función es correcta se ve la foto en la ventana
                r= colorv[0];
                g=colorv[1];
                b=colorv[2];    
            } else {
                r = optr->rgbptr[0];
                g = optr->rgbptr[1];
                b = optr->rgbptr[2];
            }   
            r += triangeluaren_I[0]*255;
            g += triangeluaren_I[1]*255;
            b += triangeluaren_I[2]*255;
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;
            glColor3ub(r,g,b);
            glVertex3f(p1x, y, p1z );
        
            glEnd();
        }

        
        
    }
    

}

void triangelua_bete_kasu_berezia(triobj *optr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr, float tartea) {
    //printf("KASU BEREZIA\n");
    //printf("Goikoa: %f, %f, %f\n Erdikoa: %f, %f %f\n Behekoa: %f, %f, %f\n", (*pgoiptrptr)->x, (*pgoiptrptr)->y, (*pgoiptrptr)->z, (*perdiptrptr)->x, (*perdiptrptr)->y, (*perdiptrptr)->z, (*pbeheptrptr)->x, (*pbeheptrptr)->y, (*pbeheptrptr)->z);
    float goi, behe, erdi;
    goi = (*pgoiptrptr)->y;
    erdi = (*perdiptrptr)->y;
    behe = (*pbeheptrptr)->y;
    float p1x, p1z, p1u, p1v, p2x, p2z, p2u, p2v;
    float yFinal1, yFinal2;
    float y;
    unsigned char r,g,b;
    unsigned char *colorv;

    for (y = goi; y>=behe; y-=tartea) {
        //printf("y: %f\n", y);

        // PUNTUEN ARTEKO ZUZENAREN FUNTZIOA: (x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
        yFinal1 = (y - (*pgoiptrptr)->y) / ((*pbeheptrptr)->y - (*pgoiptrptr)->y);  // (y-y1)/(y2-y1)
        p1x = yFinal1*((*pbeheptrptr)->x - (*pgoiptrptr)->x) + (*pgoiptrptr)->x;    // [(y-y1)/(y2-y1)]*(x2-x1) + x1
        p1z = yFinal1*((*pbeheptrptr)->z - (*pgoiptrptr)->z) + (*pgoiptrptr)->z;    // [(y-y1)/(y2-y1)]*(z2-z1) + z1
        p1u = yFinal1*((*pbeheptrptr)->u - (*pgoiptrptr)->u) + (*pgoiptrptr)->u;    // [(y-y1)/(y2-y1)]*(u2-u1) + u1
        p1v = yFinal1*((*pbeheptrptr)->v - (*pgoiptrptr)->v) + (*pgoiptrptr)->v;    // [(y-y1)/(y2-y1)]*(v2-v1) + v1
        
       
        glBegin( GL_POINTS );

        if (optr->rgbptr == 0){
            colorv=  color_textura(p1u, p1v); // si esta función es correcta se ve la foto en la ventana
            r= colorv[0];
            g=colorv[1];
            b=colorv[2];    
        } else {
            r = optr->rgbptr[0];
            g = optr->rgbptr[1];
            b = optr->rgbptr[2];
        }
        r += triangeluaren_I[0]*255;
        g += triangeluaren_I[1]*255;
        b += triangeluaren_I[2]*255;
        if (r > 255) r = 255;
        if (g > 255) g = 255;
        if (b > 255) b = 255;
        glColor3ub(r,g,b);
        glVertex3f(p1x, y, p1z );
    
        glEnd();

        
    } 
 

}

void ordenatu_x(punto *p1ptr, punto *p2ptr, punto *p3ptr, punto **pgoiptrptr, punto **pbeheptrptr, punto **perdiptrptr) {
    
    if (p1ptr->x > p2ptr->x) {
        *pgoiptrptr=p1ptr;
        *pbeheptrptr=p2ptr;
    } else {
        *pgoiptrptr=p2ptr;
        *pbeheptrptr=p1ptr;
    }
    
    if (p3ptr->x > (*pgoiptrptr)->x ){
        *perdiptrptr=*pgoiptrptr;
        *pgoiptrptr=p3ptr;
    } else {
        if ((*pbeheptrptr)->x > p3ptr->x ) {
            *perdiptrptr=*pbeheptrptr;
            *pbeheptrptr=p3ptr;
        } else {
            *perdiptrptr = p3ptr;
        }
    }
       
}


static void marraztu(void)
{
float u,v;
int i,j;
triobj *auxptr;
double nBerria[3];
/*
unsigned char* colorv;
unsigned char r,g,b;
*/



// clear viewport...
if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
    else 
      {
      if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
      }


 // kamerak marraztu
 /*
for (auxptr =fCamptr; auxptr != 0; auxptr = auxptr->hptr)
    {
        if (auxptr!=selCam_ptr) {
            for (i =0; i < auxptr->num_triangles; i++)
            {
                dibujar_triangulo(auxptr,i);
            }
        }
    }
*/
//argien bektoreak
/*
for (i=0; i<4; i++) {
    if (argiak_ptr[i].fokua == 0) {
        glBegin(GL_LINES);
            glVertex3d(argiak_ptr[i].mptr->m[3], argiak_ptr[i].mptr->m[7], argiak_ptr[i].mptr->m[11]);
            glVertex3d(argiak_ptr[i].mptr->m[3] + argiak_ptr[i].dir[0], argiak_ptr[i].mptr->m[7] + argiak_ptr[i].dir[1], argiak_ptr[i].mptr->m[11] + argiak_ptr[i].dir[2]);
        glEnd();
    } else {
        if (i == 2) auxptr = sel_ptr;
        else if (i==3) auxptr = selCam_ptr;
        glBegin(GL_LINES);
            glVertex3d(auxptr->mptr->m[3], auxptr->mptr->m[7], auxptr->mptr->m[11]);
            glVertex3d(auxptr->mptr->m[3] + argiak_ptr[i].F[0], auxptr->mptr->m[7] + argiak_ptr[i].F[1], auxptr->mptr->m[11] + argiak_ptr[i].F[2]);
        glEnd();
    }
    
}
*/
 // marrazteko objektuak behar dira
  // no se puede dibujar sin objetos
if (foptr ==0) {glFlush(); return;}
/*
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
glOrtho(-500.0, 500.0, -500.0, 500.0,-500.0, 500.0);
*/
//glFrustum(-0.5, 0.5, -0.5, 0.5,0.5, 1000.0);
// por ahora dibujamos todos los pixels de la ventana de 500x500 con el color que devuelve la función color_textura
// pero hay que llamar a la función que dibuja un triangulo con la textura mapeada:
//printf("marraztera\n");
//printf("AKI\n");
triangulosptr = sel_ptr->triptr;
//printf("AKI2\n");
if (objektuak == 1)
    {
    if (denak == 1)
        {
        for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr)
            {   // objektu bat kamera bezala erabiltzean, objektu hori ez marraztu perspektibako proiekzioan
                if (!(kameraOBJ==1 && auxptr==sel_ptr)) {
                    for (i =0; i < auxptr->num_triangles; i++)
                        {
                            dibujar_triangulo(auxptr,i);
                        }
                }
            }
        }
      else
        {
        for (i =0; i < sel_ptr->num_triangles; i++)
            {
            dibujar_triangulo(sel_ptr,i);
            }
        }
    }
  else 
    {
     dibujar_triangulo(sel_ptr,indexx);
    }
/*
for (i=0;i<500;i++)
    for (j=0;j<500;j++)
        {
        u = i/500.0;
        v = j/500.0;
        colorv=  color_textura(u, v); // si esta función es correcta se ve la foto en la ventana
        r= colorv[0];
        g=colorv[1];
        b=colorv[2];     
	glBegin( GL_POINTS );
	glColor3ub(r,g,b);
	glVertex3f(i,j,0.);
	glEnd();
	}
*/
glFlush();
}



void read_from_file(char *fitx, int isCam)
{
int i,retval;
triobj *optr;

    printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    optr->rgbptr=0;
    retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &(optr->rgbptr));
    // TODO (transform...)
    // int cargar_triangulos_color(char *fitxiz, int *hkopptr, hiruki **hptrptr,unsigned char **rgbptr);
    // retval = cargar_triangulos_color(...)
    //printf("retval: %d", retval);
    if (retval !=15 && retval != 9) 
         {
         printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n",fitxiz);
         free(optr);
         }
       else
         {
         triangulosptr = optr->triptr;
         printf("objektuaren matrizea...\n");
         optr->mptr = (mlist *)malloc(sizeof(mlist));
         for (i=0; i<16; i++) optr->mptr->m[i] =0;
         
         N_kalkulatu(optr);
         printf("objektu zerrendara doa informazioa...\n");
         if (isCam == 0) {
            optr->mptr->m[0] = 1.0;
            optr->mptr->m[5] = 1.0;
            optr->mptr->m[10] = 1.0;
            optr->mptr->m[15] = 1.0;
            optr->mptr->hptr = 0;

            optr->hptr = foptr;
            foptr = optr;
            sel_ptr = optr;
         } else if (isCam == 1) {
            optr->mptr->hptr = 0;
            optr->hptr = fCamptr;
            fCamptr = optr;
            selCam_ptr = optr;
         }
         
         }
     printf("datuak irakurrita\n");
}


void N_kalkulatu(triobj *optr) {
    int ti;
    hiruki *tptr;
    double v1x, v1y, v1z, v2x, v2y, v2z, nx, ny, nz, nnorm;

    for (ti = 0; ti<optr->num_triangles; ti++) {
        tptr = optr->triptr+ti;

        v1x = tptr->p2.x - tptr->p1.x;
        v1y = tptr->p2.y - tptr->p1.y;
        v1z = tptr->p2.z - tptr->p1.z;

        v2x = tptr->p3.x - tptr->p1.x;
        v2y = tptr->p3.y - tptr->p1.y;
        v2z = tptr->p3.z - tptr->p1.z;

        biderketa_bektoriala(v1x, v1y, v1z, v2x, v2y, v2z, &nx, &ny, &nz);
        nnorm = sqrt(pow(nx, 2.0) + pow(ny, 2.0) + pow(nz, 2.0));
        tptr->N[0] = nx/nnorm;
        tptr->N[1] = ny/nnorm;
        tptr->N[2] = nz/nnorm;

    }

}

void x_aldaketa(int dir)
{
    double m[16];
    int i, aldaketa_eginda;
    float desp, alpha;
    aldaketa_eginda = 0;

    if ((kontrola == 0 && foptr != 0) || (kontrola == 1 && fCamptr != 0) || kontrola == 2) {
        if ((aldaketa == 't' && kontrola == 0) || (aldaketa == 't' && kontrola==2 && argia_index==1)) {
            aldaketa_eginda = 1;
            if (dir == 1) {
                desp = desplazamendua;
            } else {
                desp = -desplazamendua;
            }
            for (i = 0; i<16; i++) {
                if (i==0 || i==5 || i==10 || i==15) {
                    m[i] = 1;
                } else if (i==3) {
                    m[i] = desp;
                } else {
                    m[i] = 0;
                }
            }
        
        } else if ((aldaketa == 'r' || kontrola == 1) && kontrola != 2) {
            aldaketa_eginda = 1;
            if (dir == 1) {
                alpha = 0.1;
            } else {
                alpha = -0.1;
            }
            if (analisi == 1 && kontrola == 1 && foptr != 0) {
                analisi_biraketa(selCam_ptr->mptr->m[0], selCam_ptr->mptr->m[4], selCam_ptr->mptr->m[8], alpha);
                return;
            }
            for (i = 0; i<16; i++) {
                if (i==0 || i==15) {
                    m[i] = 1;
                } else if (i==5 || i == 10) {
                    m[i] = cos(alpha);
                } else if (i==6){
                    m[i] = -sin(alpha);
                } else if (i==9){
                    m[i] = sin(alpha);
                } else {
                    m[i] = 0;
                }
            }
            
        }

        if (aldaketa_eginda == 0) {
            return;
        }
        if ((ald_lokala == 1 && kontrola == 0) || kontrola == 1 ) {
            objektuari_aldaketa_sartu_esk(m);
        } else {
            objektuari_aldaketa_sartu_ezk(m);
        }
    }

}


void y_aldaketa(int dir)
{
    double m[16], jatorrira[16], lekura_itzuli[16], mXjatorrira[16], lekuraXmXjatorrira[16];
    int i, aldaketa_eginda;
    float desp, alpha;
    aldaketa_eginda = 0;

    if ((kontrola == 0 && foptr != 0) || (kontrola == 1 && fCamptr != 0) || kontrola == 2) {
        if ((aldaketa == 't' && kontrola == 0) || (aldaketa == 't' && kontrola==2 && argia_index==1)) {
            aldaketa_eginda = 1;
                if (dir == 1) {
                    desp = desplazamendua;
                } else {
                    desp = -desplazamendua;
                }
                for (i = 0; i<16; i++) {
                    if (i==0 || i==5 || i==10 || i==15) {
                        m[i] = 1;
                    } else if (i==7) {
                        m[i] = desp;
                    } else {
                        m[i] = 0;
                    }
                }
            
        } else if ((aldaketa == 'r' || kontrola == 1) && kontrola != 2) {
            aldaketa_eginda = 1;
            if (dir == 1) {
                alpha = 0.1;
            } else {
                alpha = -0.1;
            }
            if (analisi == 1 && kontrola == 1 && foptr != 0) {
                analisi_biraketa(selCam_ptr->mptr->m[1], selCam_ptr->mptr->m[5], selCam_ptr->mptr->m[9], alpha);
                return;
            }
            for (i = 0; i<16; i++) {
                if (i==5 || i==15) {
                    m[i] = 1;
                } else if (i==0 || i == 10) {
                    m[i] = cos(alpha);
                } else if (i==8){
                    m[i] = -sin(alpha);
                } else if (i==2){
                    m[i] = sin(alpha);
                } else {
                    m[i] = 0;
                }
            }
            
        }

        if (aldaketa_eginda == 0) {
            return;
        }


        if ((ald_lokala == 1 && kontrola == 0) || kontrola == 1 ) {
            objektuari_aldaketa_sartu_esk(m);
        } else {
            objektuari_aldaketa_sartu_ezk(m);
        }
    }
}



void z_aldaketa(int dir)
{
    double m[16];
    int i, aldaketa_eginda;
    float desp, alpha;
    aldaketa_eginda = 0;
    
    if ((kontrola == 0 && foptr != 0) || (kontrola == 1 && fCamptr != 0) || kontrola == 2) {
        if ((aldaketa == 't' || kontrola == 1) || (aldaketa == 't' && kontrola==2 && argia_index==1)) {
            aldaketa_eginda = 1;
            if (dir == 1) {
                desp = desplazamendua;
            } else {
                desp = -desplazamendua;
            }
            if (kontrola == 1) desp = -desp;

            for (i = 0; i<16; i++) {
                if (i==0 || i==5 || i==10 || i==15) {
                    m[i] = 1;
                } else if (i==11) {
                    m[i] = desp;
                } else {
                    m[i] = 0;
                }
            }
            
        } else if ((aldaketa == 'r' && kontrola == 0) || (kontrola==2 && argia_index==0)) {
            aldaketa_eginda = 1;
            if (dir == 1) {
                alpha = 0.1;
            } else {
                alpha = -0.1;
            }
            for (i = 0; i<16; i++) {
                if (i==10 || i==15) {
                    m[i] = 1;
                } else if (i==0 || i == 5) {
                    m[i] = cos(alpha);
                } else if (i==1){
                    m[i] = -sin(alpha);
                } else if (i==4){
                    m[i] = sin(alpha);
                } else {
                    m[i] = 0;
                }
            }
            
        }

        if (aldaketa_eginda == 0) {
            return;
        }
        if ((ald_lokala == 1 && kontrola == 0) || kontrola == 1 ) {
            objektuari_aldaketa_sartu_esk(m);
        } else {
            objektuari_aldaketa_sartu_ezk(m);
        }
    }
    
}

int talka() {
    double vx, vy, vz, vnorm;
    vx = sel_ptr->mptr->m[3] - selCam_ptr->mptr->m[3];
    vy = sel_ptr->mptr->m[7] - selCam_ptr->mptr->m[7];
    vz = sel_ptr->mptr->m[11] - selCam_ptr->mptr->m[11];
    vnorm = sqrt(pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0));
    //printf("DISTANTZIA: %f\n ", vnorm);

    if (vnorm < desplazamendua) { // 0.1 aldaketetan egiten den desplazamendua da
        return 1;
    }
    return 0;
}

void analisi_biraketa(double vx, double vy, double vz, double alpha) {
    int i;
    double m[16], jatorrira[16], lekura_itzuli[16], mXjatorrira[16], lekuraXmXjatorrira[16], cos_alpha, sin_alpha;

    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);

    for (i = 0; i<16; i++) {
        m[i] = 0;
        if (i==0 || i==5 || i==10 || i==15) {
            jatorrira[i] = 1;
            lekura_itzuli[i] = 1;
        } else {
            jatorrira[i] = 0;
            lekura_itzuli[i] = 0;
        }
    }

    // Rodrigues matrizea
    m[15] = 1;
    m[0] = cos_alpha + (1 - cos_alpha)*pow(vx, 2.0);
    m[5] = cos_alpha + (1 - cos_alpha)*pow(vy, 2.0);
    m[10] = cos_alpha + (1 - cos_alpha)*pow(vz, 2.0);

    m[1] = (1 - cos_alpha)*vx*vy - vz*sin_alpha;
    m[4] = (1 - cos_alpha)*vx*vy + vz*sin_alpha;
    m[2] = (1 - cos_alpha)*vx*vz + vy*sin_alpha;
    m[8] = (1 - cos_alpha)*vx*vz - vy*sin_alpha;
    m[6] = (1 - cos_alpha)*vy*vz - vx*sin_alpha;
    m[9] = (1 - cos_alpha)*vy*vz + vx*sin_alpha;

    // jatorri aldaketako matrizeak
    jatorrira[3] = -sel_ptr->mptr->m[3];
    lekura_itzuli[3] = sel_ptr->mptr->m[3];
    jatorrira[7] = -sel_ptr->mptr->m[7];
    lekura_itzuli[7] = sel_ptr->mptr->m[7];
    jatorrira[11] = -sel_ptr->mptr->m[11];
    lekura_itzuli[11] = sel_ptr->mptr->m[11];

    mxm(mXjatorrira, m, jatorrira);
    mxm(lekuraXmXjatorrira, lekura_itzuli, mXjatorrira);
    objektuari_aldaketa_sartu_ezk(lekuraXmXjatorrira);
  
}

void eskala_aldatu(int handitu)
{
    double m[16];
    int i;
    float tam;

    if (foptr != 0) {
        if (handitu == 1) {
            tam = 1.5;
        } else {
            tam = 1/1.5;
        }
        for (i = 0; i<16; i++) {
            if (i==0 || i==5 || i==10) {
                m[i] = tam;
            } else if (i==15) {
                m[i] = 1;
            } else {
                m[i] = 0;
            }
        }
        
        objektuari_aldaketa_sartu_esk(m);
        
    }
}


void undo()
{
    mlist * lag;
    triobj *auxfoptr, *selAux;
    if (kontrola == 1) {
        auxfoptr = fCamptr;
        selAux = selCam_ptr;
    } else if (kontrola == 0){
        auxfoptr = foptr;
        selAux = sel_ptr;
    } else if (kontrola == 2) {
        undo_argiak();
        return;
    }

    if (auxfoptr != 0) {
        if (selAux-> mptr -> hptr != 0){
            lag = (selAux->mptr->hptr);
            selAux->mptr = lag;
        }/* else {
            printf("CANT UNDO\n");
        } */
        fokuak_kalkulatu(0);
    }
}

void undo_argiak() {
    koordlist *lag;
    
    if (argia_index == 0) { // eguzkia
        if (selArgi_ptr->dir->hptr != 0){
            lag = (selArgi_ptr->dir->hptr);
            selArgi_ptr->dir = lag;
            //printf("x: %f, y: %f, z: %f\n", selArgi_ptr->dir->koord[0], selArgi_ptr->dir->koord[1], selArgi_ptr->dir->koord[2]);
        }
    } else if (argia_index == 1) { // bonbila
        if (selArgi_ptr->pos->hptr != 0){
            lag = (selArgi_ptr->pos->hptr);
            selArgi_ptr->pos = lag;
            //printf("x: %f, y: %f, z: %f\n", selArgi_ptr->pos->koord[0], selArgi_ptr->pos->koord[1], selArgi_ptr->pos->koord[2]);
        }
    }
    
}

void delete_obj() {
    hiruki *tptr;
    mlist *mlag, *mlag2;
    triobj *obj_auxptr, *objptr_iterate, *faux, *sel_aux;

    if (kontrola == 0){
        faux = foptr;
        sel_aux = sel_ptr;
    } else if (kontrola == 1) {
        faux = fCamptr;
        sel_aux = selCam_ptr;
    }

    if ((kontrola == 0 && foptr != 0) || (kontrola == 1 && fCamptr->hptr != 0)) {
        
        // triangelu zerrenda askatu
        tptr = sel_aux->triptr;
        free(tptr);
        
        // matrize zerrenda askatu
        mlag =sel_aux->mptr;
        while (mlag != 0){
            mlag2 = mlag->hptr;
            free(mlag);
            mlag = mlag2;
        }

        // objetuen zerrenda zeharkatu, aukeratutakoa aurkitu, aurreko objetua sel_ptr-ren(edo selCam_ptr) hurrengora apuntatu, eta sel_ptr(edo selCam_ptr) askatu
       if (sel_aux == faux) {
            obj_auxptr = sel_aux->hptr;
            faux = obj_auxptr;
            free(sel_aux->rgbptr);
            free(sel_aux);        
            sel_aux = obj_auxptr;
            if (sel_aux == 0) {
                faux = 0;
                indexx =0; 
                glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT ); // azken objektua ezabatu lehiotik
                glFlush();
            } 
        } else {    
            objptr_iterate = faux;
            while (objptr_iterate->hptr != 0){

                if (objptr_iterate->hptr == sel_aux){
                    
                    obj_auxptr = sel_aux->hptr;
                    free(sel_aux->rgbptr);
                    free(sel_aux);
                    if (obj_auxptr != 0){
                        sel_aux = obj_auxptr;
                    } else {
                        sel_aux = faux;
                    }
                    
                    objptr_iterate->hptr = obj_auxptr;
                    break;

                } else {
                    objptr_iterate = objptr_iterate->hptr; 
                }
            }
        }
        
        if (kontrola == 0){
            foptr = faux;
            sel_ptr = sel_aux;
        } else if (kontrola == 1){
            fCamptr = faux;
            selCam_ptr = sel_aux;
        }     
    }
    
}

void kamerak_hasieratu() {
    int i, kont;
    triobj *auxptr;

    for (i=0; i<3; i++) {
        read_from_file("kam-1+1.txt", 1);
    }

     
    kont = 0;
    for (auxptr =fCamptr; auxptr != 0; auxptr = auxptr->hptr) {
        
        switch (kont) {
            case 2:
                auxptr->mptr->m[2] = 1;
                auxptr->mptr->m[3] = 0.99;
                auxptr->mptr->m[5] = 1;
                auxptr->mptr->m[8] = -1;
                auxptr->mptr->m[15] = 1;
                break;
            case 1:
                auxptr->mptr->m[1] = -1;
                auxptr->mptr->m[6] = 1;
                auxptr->mptr->m[7] = 0.99;
                auxptr->mptr->m[8] = -1;
                auxptr->mptr->m[15] = 1;
                break;
            case 0:
                auxptr->mptr->m[0] = 1.0;
                auxptr->mptr->m[5] = 1.0;
                auxptr->mptr->m[10] = 1.0;
                auxptr->mptr->m[15] = 1.0;

                auxptr->mptr->m[3] = 0;
                auxptr->mptr->m[7] = 0;
                auxptr->mptr->m[11] = 2;
                break;
        }
        kont++;
        
    }

}

void piztu_itzali(int index) {
    if (argiak_ptr[index].piztu == 1) argiak_ptr[index].piztu = 0;
    else argiak_ptr[index].piztu = 1;

    print_egoerak();
}


// This function will be called whenever the user pushes one key
static void teklatua (unsigned char key, int x, int y)
{
int retval;
int i;
FILE *obj_file;
triobj *auxptr;

switch(key)
	{
	case 13: 
	        if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
	                         // si no hay objeto que no haga nada
	            {
	            indexx ++;  // azkena bada lehenengoa bihurtu
		                // pero si es el último? hay que controlarlo!
		    if (indexx == sel_ptr->num_triangles) 
		        {
		        indexx = 0;
		        if ((denak == 1) && (objektuak == 0))
		            {
		            glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
		            glFlush();
		            }
		        }
		    }
		break;
    case '1':
        piztu_itzali(0);
        break;
    case '2':
        piztu_itzali(1);
        break;
    case '3':
        piztu_itzali(2);
        break;
    case '4':
        piztu_itzali(3);
        break;
	case 'd':
		if (denak == 1) denak = 0;
		    else denak = 1;
        print_egoerak();
		break;
    case 'p':
        if (perspektiba == 1) perspektiba = 0;
		    else perspektiba = 1;
        print_egoerak();
		break;
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
        print_egoerak();
		break;
	case 'c':
        if (kameraOBJ == 0){
            if (kontrola == 1) {
                kontrola = 2; printf("Argiak kontrolatzen\n"); 
                if (argia_index<2) ald_lokala=0;
                if (argia_index == 0) aldaketa = 'r'; else if (argia_index == 1) aldaketa = 't';
                } else if (kontrola == 0){ kontrola = 1; printf("Kamera kontrolatzen\n");}
                else if (kontrola == 2) {kontrola = 0; printf("Objektuak kontrolatzen\n");}
        print_egoerak();
        }
		break;
    case 'C':
		if (kameraOBJ == 1){ kameraOBJ = 0; printf("KameraOBJ=0\n");}
		    else {kameraOBJ = 1; kontrola = 0; analisi = 0; ald_lokala = 1; printf("KameraOBJ=1\n");}
        print_egoerak();
        mESA_eguneratu();
		break;
    case 'b':
        if (back_culling == 1){ back_culling = 0; printf("Back culling desaktibatuta\n");}
		    else {back_culling = 1; printf("Back culling aktibatuta\n");}
        print_egoerak();
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
        print_egoerak();
		break;
	case 't':
            if (kontrola != 2 || argia_index!=0) aldaketa = 't';
            print_egoerak();
		break;
	case 'r':
		if (kontrola != 2 || argia_index!=1) aldaketa = 'r';
        print_egoerak();
		break;
	case 'g':
        if (kontrola!=2) {
            if (kontrola == 0 && kameraOBJ == 0) {
                if (ald_lokala == 1) ald_lokala = 0;
                else ald_lokala = 1;
            } else if (kontrola == 1) {
                if (analisi == 1) analisi = 0;
                else analisi = 1;
            }
            print_egoerak();
            if (analisi==1 && foptr!=0) analisi_bektoreak();
            if (kontrola == 1 || kameraOBJ == 1) mESA_eguneratu();
        }
		break;
        case 'x':
                if (kontrola == 2 && argia_index == 0) {z_aldaketa(1); break;}
                if (kontrola != 1) x_aldaketa(1); else y_aldaketa(1);
                if (kontrola == 0 && analisi==1 && foptr!=0) analisi_bektoreak();
                if (kontrola == 1 || kameraOBJ == 1 || (kontrola == 0 && analisi==1 && foptr!=0)) mESA_eguneratu();
                break;
        case 'y':
                if (kontrola == 2 && argia_index == 0) break;
                if (kontrola != 1) y_aldaketa(1); else x_aldaketa(1);
                if (kontrola == 0 && analisi==1 && foptr!=0) analisi_bektoreak();
                if (kontrola == 1 || kameraOBJ == 1 || (kontrola == 0 && analisi==1 && foptr!=0)) mESA_eguneratu();
                break;
        case 'z':
                if (kontrola == 2 && argia_index == 0) break;
                if (kontrola == 1 && analisi == 1 && talka() == 1) break;
                z_aldaketa(1);
                if (kontrola == 0 && analisi==1 && foptr!=0) analisi_bektoreak();
                if (kontrola == 1 || kameraOBJ == 1 || (kontrola == 0 && analisi==1 && foptr!=0)) mESA_eguneratu();
                break;
        case 'X':
                if (kontrola == 2 && argia_index == 0) {z_aldaketa(0); break;}
                if (kontrola != 1) x_aldaketa(0); else y_aldaketa(0);
                if (kontrola == 0 && analisi==1 && foptr!=0) analisi_bektoreak();
                if (kontrola == 1 || kameraOBJ == 1 || (kontrola == 0 && analisi==1 && foptr!=0)) mESA_eguneratu();
                break;
        case 'Y':
                if (kontrola == 2 && argia_index == 0) break;
                if (kontrola != 1) y_aldaketa(0); else x_aldaketa(0);
                if (kontrola == 0 && analisi==1 && foptr!=0) analisi_bektoreak();
                if (kontrola == 1 || kameraOBJ == 1 || (kontrola == 0 && analisi==1 && foptr!=0)) mESA_eguneratu();
                break;
        case 'Z':
                if (kontrola == 2 && argia_index == 0) break;
                z_aldaketa(0);
                if (kontrola == 0 && analisi==1 && foptr!=0) analisi_bektoreak();
                if (kontrola == 1 || kameraOBJ == 1 || (kontrola == 0 && analisi==1 && foptr!=0)) mESA_eguneratu();
                break;
        case '+':
                if (kontrola == 0 && kameraOBJ == 0){
                    eskala_aldatu(1);
                    if (kameraOBJ == 1) mESA_eguneratu();
                } else if (kontrola == 2 && selArgi_ptr->fokua == 1) {
                    foko_angelu_aldaketa(1);
                } else if (kontrola == 1 || kameraOBJ == 1) {
                    ikuste_bol_aldaketa += 0.01;
                } 
                break;
        case '-':
                if (kontrola == 0 && kameraOBJ == 0){
                    eskala_aldatu(0);
                    if (kameraOBJ == 1) mESA_eguneratu();
                } else if (kontrola == 2 && selArgi_ptr->fokua == 1) {
                    foko_angelu_aldaketa(0);
                } else if (kontrola == 1 || kameraOBJ == 1) {
                    ikuste_bol_aldaketa -= 0.01;
                } 
                break;
        case 'u':
                undo();
                if (kontrola == 0 && analisi==1 && foptr!=0) analisi_bektoreak();
                if (kontrola == 1 || kameraOBJ == 1 || (kontrola == 0 && analisi==1 && foptr!=0)) mESA_eguneratu();
                break;
	case 'f':
	        /*Ask for file*/
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
	        read_from_file(fitxiz, 0);
	        indexx = 0;
            if (analisi==1 && foptr!=0) {
                analisi_bektoreak();
            }
            mESA_eguneratu();
                break;
       /* case 'S':  // save to file
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
                if ((obj_file = fopen(fitxiz, "w")) == NULL)
                         {
                         printf("ezin fitxategia ireki\n");
                         }
                     else
                         {
                         for (i =0; i < sel_ptr->num_triangles; i++)
                            {
                            fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                                 sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z, 
                                 sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                                 sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z, 
                                 sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                                 sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z, 
                                 sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                            }
                         fclose(obj_file);
                         }
                break; */
        case 9: /* <TAB> */
        if (kontrola == 1){ // kamera kontrolatzen
            if (fCamptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
                {
                selCam_ptr = selCam_ptr->hptr;
                /*The selection is circular, thus if we move out of the list we go back to the first element*/
                if (selCam_ptr == 0) selCam_ptr = fCamptr;
                indexx =0; // the selected polygon is the first one
                }
        } else if (kontrola == 0){ // objektuak kontrolatzen
            if (foptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
                {
                sel_ptr = sel_ptr->hptr;
                /*The selection is circular, thus if we move out of the list we go back to the first element*/
                if (sel_ptr == 0) sel_ptr = foptr;
                indexx =0; // the selected polygon is the first one
                }
        } else if (kontrola == 2) { // argiak kontrolatzen
            argia_index++;
            if (argia_index >= 4) argia_index = 0;
            selArgi_ptr = &(argiak_ptr[argia_index]);
            if (argia_index<2) ald_lokala=0;
            if (argia_index == 0) aldaketa = 'r'; else if (argia_index == 1) aldaketa = 't';
            print_egoerak();
        }
        if (analisi==1 && foptr!=0) {
            analisi_bektoreak();
        }
        mESA_eguneratu();
            break;
        case 127: // <SUPR>
            delete_obj();
            if (analisi==1 && foptr!=0) {
                analisi_bektoreak();
            }
            mESA_eguneratu();
            break;
	case 27:  // <ESC>
		exit( 0 );
		break;
	default:
		printf("%d %c\n", key, key );
	}
    
    
    
// The screen must be drawn to show the new triangle
glutPostRedisplay();
}

void print_egoerak() {

    printf("\nEGOERAK:\n[d]Obj_guztiak: %d,                       [o]Obj_osoa: %d,                [l]Lineak: %d,\n", denak, objektuak, lineak);
    printf("[t]Aldaketak(t)/[r]biraketak(r): %c,      [g]Lokala(1)/globala(0): %d,    [c]Kontrolatzen: argia(2)/kamera(1)/objektua(0): %d\n", aldaketa, ald_lokala, kontrola);
    printf("[p]Perspektiba(1)/paraleloa(0): %d,       [C]Objektuaren ikuspegia: %d,   [g, kamera moduan]Analisia(1)/hegaldia(0): %d\n", perspektiba, kameraOBJ, analisi);
    printf("[b]Back_culling: %d\n", back_culling);

    if (kameraOBJ == 1) {
        printf("AUKERATUTAKO OBJEKTUAREN IKUSPEGIA (OBJEKTUA EZ DA ERAKUSTEN)\n");
    } else printf("AUKERATUTAKO KAMERAREN IKUSPEGIA\n");
    printf("AUKERATUTAKO ARGIA: ");
    switch (argia_index)
    {
    case 0:
        printf("EGUZKIA\n");
        break;
    case 1:
        printf("BONBILLA\n");
        break;
    case 2:
        printf("OBJEKTUAREN FOKOA\n");
        break;
    case 3:
        printf("KAMERAREN FOKOA\n");
        break;
    default:
        break;
    }
    print_piztutako_argiak();
    
}

void print_piztutako_argiak() {
    printf("PIZTUTAKO ARGIAK[1..4]: ");
    if (argiak_ptr[0].piztu == 1) printf("[1]EGUZKIA ");
    if (argiak_ptr[1].piztu == 1) printf("[2]BONBILLA ");
    if (argiak_ptr[2].piztu == 1) printf("[3]OBJEKTUA ");
    if (argiak_ptr[3].piztu == 1) printf("[4]KAMERA ");
    printf("\n");
}


void viewportberria (int zabal, int garai)
{
if (zabal < garai)  dimentsioa = zabal;
    else  dimentsioa = garai;
glViewport(0,0,dimentsioa,dimentsioa);
printf("linea kopuru berria = %d\n",dimentsioa);
}
       
int main(int argc, char** argv)
{
int retval;


	printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	dimentsioa = 500;
	glutInitWindowSize ( dimentsioa, dimentsioa );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow( "KBG/GO praktika" );

	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	glutReshapeFunc( viewportberria);
	/* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */ 
        retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
        if (retval !=1) 
            {
            printf("Ez dago testuraren fitxategia (testura.ppm)\n");
            exit(-1);
            }
        
	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        glDepthFunc(GL_GREATER); // Handiena marraztu dezan
        glClearDepth(0.0); // Buferra hasieratzeko balioa 0 izan dadila 1 izan beharrean.

        denak = 1;
        lineak = 0;
        objektuak = 1;
        kontrola = 0;
        kameraOBJ = 0;
        foptr = 0;
        sel_ptr = 0;
        fCamptr = 0;
        selCam_ptr = 0;
        aldaketa = 'r';
        ald_lokala = 1;
        perspektiba = 1;
        analisi = 0;
        back_culling = 1;
        gorria = 0;
        desplazamendua = 0.1;
        ikuste_bol_aldaketa = 0.0;
        IaR = 1;
        IaG = 1;
        IaB = 1;
        KaR = 0.1;
        KaG = 0.1;
        KaB = 0.1;
        triangeluaren_I[0] = IaR * KaR;
        triangeluaren_I[1] = IaG * KaG;
        triangeluaren_I[2] = IaB * KaB;

        kamerak_hasieratu();

        if (argc>1) read_from_file(argv[1], 0);
            else {
                /*
                read_from_file("r_falke-1+1.txt", 0);
                foptr->mptr->m[3] = -0.5;
                foptr->rgbptr[0] = 0;
                foptr->rgbptr[1] = 0;
                foptr->rgbptr[2] = 0;
                read_from_file("abioia-1+1.txt", 0);
                foptr->mptr->m[3] = 0.5;
                foptr->mptr->m[11] = -0.5;
                foptr->rgbptr[0] = 0;
                foptr->rgbptr[1] = 0;
                foptr->rgbptr[2] = 0;
                */
               read_from_file("abioia-1+1.txt", 0);
                if (sel_ptr != 0) 
                    { sel_ptr->mptr->m[3] = -1.0;
                    if (sel_ptr->rgbptr !=0) sel_ptr->rgbptr[0]=0;
                    }   
                read_from_file("abioia-1+1.txt", 0);
                if (sel_ptr != 0) 
                    { sel_ptr->mptr->m[3] = 1.0;
                    if (sel_ptr->rgbptr !=0) sel_ptr->rgbptr[1]=0;
                    }   
                read_from_file("abioia-1+1.txt", 0);
                if (sel_ptr != 0) 
                    { 
                    sel_ptr->mptr->m[7] = -0.4;
                    sel_ptr->mptr->m[11] = 0.5;
                    if (sel_ptr->rgbptr !=0) sel_ptr->rgbptr[2]=0;
                    }   
                read_from_file("abioia-1+1.txt", 0);
                if (sel_ptr != 0) 
                    { sel_ptr->mptr->m[11] = -0.7;
                    if (sel_ptr->rgbptr !=0) sel_ptr->rgbptr[2]=0;
                    }        
                read_from_file("abioia-1+1.txt", 0);
                if (sel_ptr != 0) 
                    { sel_ptr->mptr->m[7] = 0.7;
                    if (sel_ptr->rgbptr !=0) sel_ptr->rgbptr[2]=0;
                    }                        
                read_from_file("abioia-1+1.txt", 0);
                    if (sel_ptr->rgbptr !=0) sel_ptr->rgbptr[2]=0;
            }
        
        
        if (fCamptr != 0) mESA_kalkulatu(selCam_ptr->mptr->m);
        argiak_hasieratu();
        selArgi_ptr = &(argiak_ptr[0]);
        argia_index=0;
        print_egoerak();
	glutMainLoop();

	return 0;   
}
