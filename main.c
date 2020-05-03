#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

typedef struct pixel
{
    unsigned char B, G, R;
} pixel;

typedef struct detectie
{
    double corelatie;
    unsigned int coord_inaltime, coord_latime;
    int culoare;

}detectie;

void xorshift (unsigned int *R, unsigned int latime_img, unsigned int inaltime_img, unsigned int seed)
{
    unsigned int i;
    R[0] = seed;
    for (i = 1; i <= 2 * latime_img * inaltime_img - 1; i++)
    {
        R[i] = R[i - 1];
        R[i] = R[i] ^ R[i] << 13;
        R[i] = R[i] ^ R[i] >> 17;
        R[i] = R[i] ^ R[i] << 5;
    }
}

void incarcare_BMP (char* nume_img, pixel **P, unsigned int *latime_img, unsigned int *inaltime_img)
{
    FILE* f = fopen(nume_img, "rb");
    fseek(f, 18, SEEK_SET);
    fread(latime_img, sizeof(unsigned int), 1, f);
    fread(inaltime_img, sizeof(unsigned int), 1, f);
    fseek(f, 54, SEEK_SET);

    unsigned int padding = 0;
    if (*latime_img % 4 != 0)
        padding = 4 - (3 * (*latime_img)) % 4;

    int  i, j;
    *P = (pixel*)malloc(*latime_img * (*inaltime_img) * sizeof(pixel));
    for (i = *inaltime_img - 1; i >= 0; i--)
    {
        for (j = 0; j < *latime_img; j++)
            fread(&(*P)[i * (*latime_img) + j], sizeof(pixel), 1, f);
        fseek(f, padding, SEEK_CUR);
    }

    fclose(f);
}

void salvare_BMP (char* nume_img, pixel *P, unsigned int latime_img, unsigned int inaltime_img, unsigned char* header_img)
{
    FILE* f = fopen(nume_img, "wb");
    int i, j;
    for (i = 0; i < 54; i++)
        fwrite(&header_img[i], sizeof(unsigned char), 1, f);

    unsigned int padding = 0;
    if (latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;

    for (i = inaltime_img - 1; i >= 0; i--)
    {
        for (j = 0; j < latime_img; j++)
            fwrite(&P[i * latime_img + j], sizeof(pixel), 1, f);
        char x = 0;
        for (j = 0; j < padding; j++)
            fwrite(&x, 1, 1, f);
    }
    fclose(f);
}


void criptare_BMP (char* nume_img_initiala, char* nume_img_criptata, char* nume_cheie_secreta, pixel *P, unsigned int latime_img, unsigned int inaltime_img)
{
    FILE* f_initiala = fopen(nume_img_initiala, "rb");
    FILE* f_cheie = fopen(nume_cheie_secreta, "r");

    unsigned int i;
    unsigned char* header_img;
    header_img = (unsigned char*)malloc(54 * sizeof(unsigned char));
    if (header_img == NULL)
    {
        printf("Eroare la alocarea memoriei pentru header");
        return;
    }
    for (i = 0; i < 54; i++)
    {
        unsigned char c;
        fread(&c, sizeof(unsigned char), 1, f_initiala);
        header_img[i] = c;
    }

    unsigned int* R;
    R = (unsigned int*)malloc(2 * latime_img * inaltime_img * sizeof(unsigned int));
    if (R == NULL)
    {
        printf("Nu s-a putut aloca memorie pentru R");
        return;
    }

    unsigned int seed;
    fscanf(f_cheie, "%u", &seed);
    xorshift(R, latime_img, inaltime_img, seed);

    unsigned int* perm = (unsigned int*)malloc(latime_img * inaltime_img * sizeof(unsigned int));
    if (perm == NULL)
    {
        printf("Nu s-a putut aloca memorie pentru perm");
        return;
    }
    for (i = 0; i < latime_img * inaltime_img; i++)
        perm[i] = i;
    for (i = latime_img * inaltime_img - 1; i >= 1; i--)
    {
        unsigned int j = R[latime_img * inaltime_img - i] %(i + 1);
        unsigned int aux = perm[i];
        perm[i] = perm[j];
        perm[j] = aux;
    }

    pixel* auxP;
    auxP = (pixel*)malloc(latime_img * inaltime_img * sizeof(pixel));
    if (auxP == NULL)
    {
        printf("Nu s-a putut aloca memorie pentru auxP");
        return;
    }

    for (i = 0; i < latime_img * inaltime_img; i++)
        auxP[perm[i]] = P[i];
    free(perm);
    for (i = 0; i < latime_img * inaltime_img; i++)
        P[i] = auxP[i];

    unsigned int SV;
    fscanf(f_cheie, "%u", &SV);
    unsigned char *octSV = (unsigned char*)malloc(4 * sizeof(unsigned char));
    for (i = 0; i < 4; i++)
    {
        octSV[i] = (unsigned char)SV;
        SV = SV >> 8;
    }

    unsigned char *octR = (unsigned char*)malloc(4 * sizeof(unsigned char));
    unsigned int copieR = R[latime_img * inaltime_img];
    for (i = 0; i < 4; i++)
    {
        octR[i] = (unsigned char)copieR;
        copieR = copieR >> 8;
    }

    auxP[0].R = octSV[2] ^ P[0].R ^ octR[2];
    auxP[0].G = octSV[1] ^ P[0].G ^ octR[1];
    auxP[0].B = octSV[0] ^ P[0].B ^ octR[0];

    for (i = 1; i < latime_img * inaltime_img; i++)
    {
        int j;
        unsigned int copieR = R[latime_img * inaltime_img + i];
        for (j = 0; j < 4; j++)
        {
            octR[j] = (unsigned char)copieR;
            copieR = copieR >> 8;
        }
        auxP[i].R = auxP[i - 1].R ^ P[i].R ^ octR[2];
        auxP[i].G = auxP[i - 1].G ^ P[i].G ^ octR[1];
        auxP[i].B = auxP[i - 1].B ^ P[i].B ^ octR[0];
    }
    free(octR);
    free(octSV);
    for (i = 0; i < latime_img * inaltime_img; i++)
        P[i] = auxP[i];
    free(auxP);
    free(R);
    fclose(f_initiala);
    fclose(f_cheie);

    salvare_BMP(nume_img_criptata, P, latime_img, inaltime_img, header_img);
    free(header_img);

}

void decriptare_BMP (char* nume_img_criptata, char* nume_img_decriptata, char* nume_cheie_secreta, pixel *P, unsigned int latime_img, unsigned int inaltime_img)
{
    FILE* f_criptata = fopen(nume_img_criptata, "rb");
    FILE* f_cheie = fopen(nume_cheie_secreta, "r");

    unsigned int i;
    unsigned char* header_img;
    header_img = (unsigned char*)malloc(54 * sizeof(unsigned char));
    if (header_img == NULL)
    {
        printf("Eroare la alocarea memoriei pentru header");
        return;
    }
    for (i = 0; i < 54; i++)
    {
        unsigned char c;
        fread(&c, sizeof(unsigned char), 1, f_criptata);
        header_img[i] = c;
    }

    unsigned int* R;
    R = (unsigned int*)malloc(2 * latime_img * inaltime_img * sizeof(unsigned int));
    if (R == NULL)
    {
        printf("Nu s-a putut aloca memorie pentru R");
        return;
    }

    unsigned int seed;
    fscanf(f_cheie, "%u", &seed);
    xorshift(R, latime_img, inaltime_img, seed);

    unsigned int* perm = (unsigned int*)malloc(latime_img * inaltime_img * sizeof(unsigned int));
    if (perm == NULL)
    {
        printf("Nu s-a putut aloca memorie pentru perm");
        return;
    }
    for (i = 0; i < latime_img * inaltime_img; i++)
        perm[i] = i;
    for (i = latime_img * inaltime_img - 1; i >= 1; i--)
    {
        unsigned int j = R[latime_img * inaltime_img - i] %(i + 1);
        unsigned int aux = perm[i];
        perm[i] = perm[j];
        perm[j] = aux;
    }

    unsigned int* perm2 = (unsigned int*)malloc(latime_img * inaltime_img * sizeof(unsigned int));
    if (perm2 == NULL)
    {
        printf("Nu s-a putut aloca memorie pentru perm2");
        return;
    }
    for (i = 0; i < latime_img * inaltime_img; i++)
        perm2[perm[i]] = i;
    free(perm);

    pixel* auxP;
    auxP = (pixel*)malloc(latime_img * inaltime_img * sizeof(pixel));
    if (auxP == NULL)
    {
        printf("Nu s-a putut aloca memorie pentru auxP");
        return;
    }

    unsigned int SV;
    fscanf(f_cheie, "%u", &SV);
    unsigned char *octSV = (unsigned char*)malloc(4 * sizeof(unsigned char));
    for (i = 0; i < 4; i++)
    {
        octSV[i] = (unsigned char)SV;
        SV = SV >> 8;
    }

    unsigned char *octR = (unsigned char*)malloc(4 * sizeof(unsigned char));
    unsigned int copieR = R[latime_img * inaltime_img];
    for (i = 0; i < 4; i++)
    {
        octR[i] = (unsigned char)copieR;
        copieR = copieR >> 8;
    }

    auxP[0].R = octSV[2] ^ P[0].R ^ octR[2];
    auxP[0].G = octSV[1] ^ P[0].G ^ octR[1];
    auxP[0].B = octSV[0] ^ P[0].B ^ octR[0];

    for (i = 1; i < latime_img * inaltime_img; i++)
    {
        int j;
        unsigned int copieR = R[latime_img * inaltime_img + i];
        for (j = 0; j < 4; j++)
        {
            octR[j] = (unsigned char)copieR;
            copieR = copieR >> 8;
        }
        auxP[i].R = P[i - 1].R ^ P[i].R ^ octR[2];
        auxP[i].G = P[i - 1].G ^ P[i].G ^ octR[1];
        auxP[i].B = P[i - 1].B ^ P[i].B ^ octR[0];
    }
    free(octR);
    free(octSV);

    for (i = 0; i < latime_img * inaltime_img; i++)
        P[perm2[i]] = auxP[i];
    free(perm2);
    free(auxP);
    free(R);
    fclose(f_criptata);
    fclose(f_cheie);

    salvare_BMP(nume_img_decriptata, P, latime_img, inaltime_img, header_img);
    free(header_img);

}

void valori_chi_square (char* nume_img, int latime_img, int inaltime_img)
{
    FILE* f = fopen(nume_img, "rb");
    pixel p;
    double *frR = (double*)calloc(256, sizeof(double));
    double *frG = (double*)calloc(256, sizeof(double));
    double *frB = (double*)calloc(256, sizeof(double));
    fseek(f, 54, SEEK_SET);
    while (fread(&p, sizeof(pixel), 1, f) == 1)
    {
        frR[p.R]++;
        frG[p.G]++;
        frB[p.B]++;
    }

    double fr_estimata = (latime_img * inaltime_img) / 256;
    double chi_squareR = 0, chi_squareG = 0, chi_squareB = 0;
    int i;
    for (i = 0; i <= 255; i++)
        chi_squareR += ((frR[i] - fr_estimata) * (frR[i] - fr_estimata)) / fr_estimata;
    for (i = 0; i <= 255; i++)
        chi_squareG += ((frG[i] - fr_estimata) * (frG[i] - fr_estimata)) / fr_estimata;
    for (i = 0; i <= 255; i++)
        chi_squareB += ((frB[i] - fr_estimata) * (frB[i] - fr_estimata)) / fr_estimata;
    printf("\nR: %0.2f\nG: %0.2lf\nB: %0.2lf\n", chi_squareR, chi_squareG, chi_squareB);

    free(frR);
    free(frG);
    free(frB);
    fclose(f);
}

void grayscale_image(pixel* P, unsigned int latime_img, unsigned int inaltime_img)
{
    int i, j;
	for (i = 0; i < inaltime_img; i++)
        for (j = 0; j < latime_img; j++)
        {
            unsigned char aux = 0.299 * P[i * latime_img + j].R + 0.587 * P[i * latime_img + j].G + 0.114 * P[i * latime_img + j].B;
            P[i * latime_img + j].R = P[i * latime_img + j].G = P[i * latime_img + j].B = aux;
        }
}

void template_matching (pixel* P, unsigned int latime_img, unsigned int inaltime_img, pixel* S, unsigned int latime_sablon, unsigned int inaltime_sablon, float ps, detectie** fereastra, int* nr_ferestre, int culoare)
{
    *fereastra = NULL;
    *nr_ferestre = 0;

    float suma_intensitati_sablon = 0;
    int i, j;
    for (i = 0; i < inaltime_sablon; i++)
        for (j = 0; j < latime_sablon; j++)
            suma_intensitati_sablon += S[i * latime_sablon + j].R;

    float medie_intensitati_sablon;
    float nr_pixeli_sablon = latime_sablon * inaltime_sablon;
    medie_intensitati_sablon = suma_intensitati_sablon / nr_pixeli_sablon;

    float deviatie_sablon = 0;
    for (i = 0; i < inaltime_sablon; i++)
        for (j = 0; j < latime_sablon; j++)
            deviatie_sablon += (S[i * latime_sablon + j].R - medie_intensitati_sablon) * (S[i * latime_sablon + j].R - medie_intensitati_sablon);
    deviatie_sablon *= (1 / (nr_pixeli_sablon - 1));
    deviatie_sablon = sqrt(deviatie_sablon);

    int x, y;
    for (x = 0; x <= inaltime_img - inaltime_sablon; x++)
        for (y = 0; y <= latime_img - latime_sablon; y++)
        {
            float suma_intensitati_fereastra = 0;
            for (i = x; i < x + inaltime_sablon; i++)
                for (j = y; j < y + latime_sablon; j++)
                    suma_intensitati_fereastra += P[i * latime_img + j].R;

            float medie_intensitati_fereastra;
            medie_intensitati_fereastra = suma_intensitati_fereastra / nr_pixeli_sablon;

            float deviatie_fereastra = 0;
            for (i = x; i < x + inaltime_sablon; i++)
                for (j = y; j < y + latime_sablon; j++)
                    deviatie_fereastra += (P[i * latime_img + j].R - medie_intensitati_fereastra) * (P[i * latime_img + j].R - medie_intensitati_fereastra);
            deviatie_fereastra *= (1 / (nr_pixeli_sablon - 1));
            deviatie_fereastra = sqrt(deviatie_fereastra);

            float corelatie = 0;
            for (i = 0; i < inaltime_sablon; i++)
                for (j = 0; j < latime_sablon; j++)
                    corelatie += (1 / (deviatie_sablon * deviatie_fereastra)) * (P[(x + i) * latime_img + y + j].R - medie_intensitati_fereastra) * (S[i * latime_sablon + j].R - medie_intensitati_sablon);
            corelatie *= (1 / nr_pixeli_sablon) ;
            if (corelatie > ps)
            {
                (*fereastra) = (detectie*)realloc(*fereastra, (*nr_ferestre + 1) * sizeof(detectie));

                if((*fereastra) == NULL)
                {
                    printf("Eroare la alocarea memoriei pentru fereastra");
                    return;
                }

                (*fereastra)[*nr_ferestre].coord_inaltime = x;
                (*fereastra)[*nr_ferestre].coord_latime = y;
                (*fereastra)[*nr_ferestre].corelatie = corelatie;
                (*fereastra)[*nr_ferestre].culoare = culoare;
                (*nr_ferestre)++;
            }
        }
}

void colorare_imagine (pixel* P, unsigned int latime_img, unsigned int latime_sablon, unsigned int inaltime_sablon, detectie fereastra, pixel culoare)
{
    int i, j;

    for (i = fereastra.coord_inaltime; i < fereastra.coord_inaltime + inaltime_sablon; i++)
        P[i * latime_img + fereastra.coord_latime] = culoare;

    for (j = fereastra.coord_latime + 1; j < fereastra.coord_latime + latime_sablon; j++)
        P[fereastra.coord_inaltime * latime_img + j] = culoare;

    for (i = fereastra.coord_inaltime + 1; i < fereastra.coord_inaltime + inaltime_sablon; i++)
        P[i * latime_img + fereastra.coord_latime + latime_sablon - 1] = culoare;

    for (j = fereastra.coord_latime + 1; j < fereastra.coord_latime + latime_sablon - 1; j++)
        P[(fereastra.coord_inaltime + inaltime_sablon - 1) * latime_img + j] = culoare;

}

int cmp_descrescator(const void* a, const void* b)
{
    detectie va = *(detectie*)a;
    detectie vb = *(detectie*)b;
    if (va.corelatie > vb.corelatie)
        return -1;
    if (va.corelatie < vb.corelatie)
        return 1;
    return 0;
}

void sortare_detectii(detectie *D, int nr_detectii)
{
    qsort(D, nr_detectii, sizeof(D[0]), cmp_descrescator);
}

void eliminare_nonmaxime(detectie *D, int nr_detectii, int latime_sablon, int inaltime_sablon)
{
    int i = 0;
    while (i < nr_detectii - 1)
    {
        int j = i + 1;
        while (j < nr_detectii)
        {
            if ((abs(D[i].coord_inaltime - D[j].coord_inaltime) >= inaltime_sablon) || (abs(D[i].coord_latime - D[j].coord_latime) >= latime_sablon))
                j++;
            else
            {
                double arie_intersectie, arie_Di, arie_Dj, arie_reuniune;
                arie_intersectie = (latime_sablon - abs(D[i].coord_latime - D[j].coord_latime)) * (inaltime_sablon - abs(D[i].coord_inaltime - D[j].coord_inaltime));
                arie_Di = arie_Dj = latime_sablon * inaltime_sablon;
                arie_reuniune = arie_Di + arie_Dj - arie_intersectie;
                double suprapunere = arie_intersectie / arie_reuniune;
                if (suprapunere > 0.2)
                {
                    int k;
                    for (k = j; k < nr_detectii - 1; k++)
                        D[k] = D[k + 1];
                    nr_detectii--;
                }
                else
                    j++;
            }
        }
        i++;
    }
}

int main()
{
//    In fisierul text, prima linie contine imaginea initiala, cea de-a doua contine imaginea criptata, cea de-a treia imaginea decriptata,
//    cea de-a patra contine fisierul cu cheia secreta, cea de-a cincea imaginea pe care se face template matching, liniile 6-15 contin
//    imaginile cifrelor.

    FILE* f = fopen("nume_fisiere.txt", "r");
    char *s = (char*)malloc(30 * sizeof(char));
    char** nume_fisiere = NULL;
    int nr = 0;
    while (fgets(s, 30, f) != NULL)
    {
        if (s[strlen(s) - 1] == '\n')
            s[strlen(s) - 1] = '\0';
        int lungime_sir = strlen(s) + 1;
        nume_fisiere = (char**)realloc(nume_fisiere, (nr + 1) * sizeof(char*));
        if (nume_fisiere == NULL)
        {
            printf("Eroare la alocarea memoriei pentru nume_fisiere");
            return 0;
        }
        nume_fisiere[nr] = (char*)malloc(lungime_sir * sizeof(char));
        strcpy(nume_fisiere[nr], s);
        nr++;
    }
    fclose(f);

    pixel *P;
    unsigned int latime_img, inaltime_img;
    incarcare_BMP(nume_fisiere[0], &P, &latime_img, &inaltime_img);
    criptare_BMP(nume_fisiere[0], nume_fisiere[1], nume_fisiere[3], P, latime_img, inaltime_img);
    decriptare_BMP(nume_fisiere[1], nume_fisiere[2], nume_fisiere[3], P, latime_img, inaltime_img);
    printf("\nValori test chi square pentru imaginea initiala (%s): ", nume_fisiere[0]);
    valori_chi_square(nume_fisiere[0], latime_img, inaltime_img);
    printf("\nValori test chi square pentru imaginea criptata (%s): ", nume_fisiere[1]);
    valori_chi_square(nume_fisiere[1], latime_img, inaltime_img);

    free(P);

    pixel *S;
    unsigned int latime_sablon, inaltime_sablon;
    detectie** fereastra;
    fereastra = (detectie**)malloc(10 * sizeof(detectie*));
    int* nr_ferestre = (int*)malloc(10 * sizeof(int));
    double ps = 0.5;

    pixel* culoare = (pixel*)malloc(10 * sizeof(pixel));

    culoare[0].R = 255;
    culoare[0].G = 0;
    culoare[0].B = 0;

    culoare[1].R = 255;
    culoare[1].G = 255;
    culoare[1].B = 0;

    culoare[2].R = 0;
    culoare[2].G = 255;
    culoare[2].B = 0;

    culoare[3].R = 0;
    culoare[3].G = 255;
    culoare[3].B = 255;

    culoare[4].R = 255;
    culoare[4].G = 0;
    culoare[4].B = 255;

    culoare[5].R = 0;
    culoare[5].G = 0;
    culoare[5].B = 255;

    culoare[6].R = 192;
    culoare[6].G = 192;
    culoare[6].B = 192;

    culoare[7].R = 255;
    culoare[7].G = 140;
    culoare[7].B = 0;

    culoare[8].R = 128;
    culoare[8].G = 0;
    culoare[8].B = 128;

    culoare[9].R = 128;
    culoare[9].G = 0;
    culoare[9].B = 0;

    incarcare_BMP(nume_fisiere[4], &P, &latime_img, &inaltime_img);
    grayscale_image(P, latime_img, inaltime_img);

    int i;
    for (i = 0; i <= 9; i++)
    {
        incarcare_BMP(nume_fisiere[i + 5], &S, &latime_sablon, &inaltime_sablon);
        grayscale_image(S, latime_sablon, inaltime_sablon);
        template_matching(P, latime_img, inaltime_img, S, latime_sablon, inaltime_sablon, ps, &fereastra[i], &nr_ferestre[i], i);
    }

    int j;
    int nr_detectii = 0;
    detectie* D = NULL;
    for (i = 0; i < 10; i++)
        for (j = 0; j < nr_ferestre[i]; j++)
        {
            D = (detectie*)realloc(D, (nr_detectii + 1) * sizeof(detectie));
            if (D == NULL)
            {
                printf("Eroare la alocarea memoriei pentru D");
                return 0;
            }
            D[nr_detectii] = fereastra[i][j];
            nr_detectii++;
        }

    sortare_detectii(D, nr_detectii);

    eliminare_nonmaxime(D, nr_detectii, latime_sablon, inaltime_sablon);

    incarcare_BMP(nume_fisiere[4], &P, &latime_img, &inaltime_img);

    for (i = 0; i < nr_detectii; i++)
        colorare_imagine(P, latime_img, latime_sablon, inaltime_sablon, D[i], culoare[D[i].culoare]);

    f = fopen(nume_fisiere[4], "rb");
    unsigned char* header_img;
    header_img = (unsigned char*)malloc(54 * sizeof(unsigned char));
    if (header_img == NULL)
    {
        printf("Eroare la alocarea memoriei pentru header");
        return 0;
    }
    for (i = 0; i < 54; i++)
    {
        unsigned char c;
        fread(&c, sizeof(unsigned char), 1, f);
        header_img[i] = c;
    }
    fclose(f);

    salvare_BMP(nume_fisiere[15], P, latime_img, inaltime_img, header_img);

    for (i = 0; i <= 14; i++)
        free(nume_fisiere[i]);
    free(nume_fisiere);
    free(s);
    free(culoare);
    for (i = 0; i < 10; i++)
        free(fereastra[i]);
    free(nr_ferestre);
    free(D);
    free(P);
    free(header_img);
    return 0;
}
