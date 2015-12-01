//
// Created by lukasa on 11/20/15.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "complex_struct.h"
#include "file_util.h"

int CalculatePoints(Complex *comp);

Simplex *simplex_intersection(Simplex *a, Simplex *b);

bool isOrdered(Simplex *simp);

bool simp_binary_search(Simplex *simp, SimplexElem search);

int main(int argc, char *argv[]) {
    clock_t begin, end;
    double time_spent;
    begin = clock();
    
    char const *const fileName = "/Users/lukasa/Documents/CC.txt";//argv[1];
    FILE *file = fopen(fileName, "r");
    char line[256];
    
    int line_number = 0;
    int fsin = 0;
    
    Complex *A, *B;
    int bPoints;
    
    Complex *poset_matrix = Init_Complex();
    int poset_matrix_size;
    
    while (fgets(line, sizeof(line), file)) {
        line_number++;
        
        if (line_number == 1) {
            const char *AL = strtok(line, "->");
            const char *BL = strtok(NULL, "->");
            
            A = literalToComplex((char *) AL);
            B = literalToComplex((char *) BL);
            
            bPoints = CalculatePoints(B);
            poset_matrix_size = bPoints * bPoints;
            
            for (int i = 0; i < poset_matrix_size + 1; ++i) {
                addSimplex(poset_matrix, Init_Simplex());
            }
            
        } else {
            Complex *fsi = literalToComplex(line);
            
            if (fsi == NULL) {
                continue;
            }
            fsin++;
            
            for (int i = 0; i < fsi->simplexCount; ++i) {
                Simplex *simp = getSimpexAt(fsi, i);
                for (int j = 0; j < simp->elementCount; ++j) {
                    SimplexElem elem = getElementAt(simp, j);
                    
                    Simplex *posetSimp = getSimpexAt(poset_matrix, (i * bPoints) + (elem - 1));
                    
                    addElement(posetSimp, fsin);
                }
            }
            
            Dest_Complex(fsi);
        }
    }
    fclose(file);
    
    FILE *dumpf = fopen("./aij.txt", "w");
    
    for (int i = 0; i < poset_matrix->simplexCount; ++i) {
        Simplex* simp = getSimpexAt(poset_matrix, i);
        fputs("{", dumpf);
        for (int j = 0; j < simp->elementCount; ++j) {
            SimplexElem elem = getElementAt(simp, j);
            char str[10];
            
            sprintf(str, "%d",  elem);
            fputs(str, dumpf);
            if (j != simp->elementIndex) {
                fputs(",", dumpf);
            }
        }
        fputs("}", dumpf);
        if (i != poset_matrix->simplexIndex) {
            fputs("\n", dumpf);
        }
    }
    
    fclose(dumpf);
    
    printf("A(i,j) Dumped!\n");
    fflush(stdout);
    
    clock_t begin1, end1;
    double time_spent1;
    begin1 = clock();
    srand((unsigned int) time(NULL));
    int printCount = 500 + (rand() % 500);
    
    Simplex *max_fsi = Init_Simplex();
    
    for (int i = 1; i <= fsin; ++i) {
        Complex *poset_col = Init_Complex();
        for (int j = 0; j < poset_matrix->simplexCount; ++j) {
            Simplex *poset = getSimpexAt(poset_matrix, j);
            
            if (simp_binary_search(poset, i)) {
                addSimplex(poset_col, poset);
            }
        }
        
        Simplex *intersection = simplex_intersection(getSimpexAt(poset_col, 0), getSimpexAt(poset_col, 1));
        for (int k = 2; k < poset_col->simplexCount; ++k) {
            Simplex *poset = getSimpexAt(poset_col, k);
            
            
            Simplex *preIntersection = simplex_intersection(intersection, poset);
            Dest_Simplex(intersection);
            intersection = preIntersection;
        }
        
        if (intersection->elementCount == 1) {
            SimplexElem max_fsi_elem = getElementAt(intersection, 0);
            addElement(max_fsi, max_fsi_elem);
        }
        Dest_Simplex(intersection);
        
        if (i >= printCount) {
            end1 = clock();
            time_spent1 = (double) (end1 - begin1) / CLOCKS_PER_SEC;
            
            int tmpRand = 500 + (rand() % 500);
            
            printf("\nTime: %f; Amount:%d; | %d/%d | MaxFound: %d\n", time_spent1, tmpRand, i, fsin,
                   max_fsi->elementCount);
            fflush(stdout);
            
            begin1 = clock();
            printCount += tmpRand;
        }
    }
    
    
    printf("\n Max FSI Indexes: %s \n", simplexToLiteral(max_fsi));
    fflush(stdout);
    fsin = 0;
    file = fopen(fileName, "r");
    while (fgets(line, sizeof(line), file)) {
        line_number++;
        
        if (line_number < 2) continue;
        Complex *fsi = literalToComplex(line);
        
        if (fsi == NULL) {
            continue;
        }
        fsin++;
        
        if (simp_binary_search(max_fsi, fsin)) {
            printf("%d: %s\n", fsin, line);
        }
    }
    fclose(file);
    
    end = clock();
    time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    
    printf("\n Time Spend: %f \n", time_spent);
    printf("The End..");
    fflush(stdout);
}


int CalculatePoints(Complex *comp) {
    SimplexElem elem = -1;
    
    for (int i = 0; i < comp->simplexCount; ++i) {
        Simplex *simp = getSimpexAt(comp, i);
        for (int j = 0; j < simp->elementCount; ++j) {
            SimplexElem elemMax = getElementAt(simp, j);
            if (elemMax > elem) {
                elem = elemMax;
            }
        }
    }
    return elem;
}

Simplex *simplex_intersection(Simplex *a, Simplex *b) {
    Simplex *c = Init_Simplex();
    
    if (a->elementCount > b->elementCount) {
        Simplex *tmp = a;
        a = b;
        b = tmp;
    }
    
    for (int i = 0; i < a->elementCount; ++i) {
        SimplexElem search = getElementAt(a, i);
        
        if (simp_binary_search(b, search)) {
            addElement(c, search);
        }
    }
    
    return c;
}

bool simp_binary_search(Simplex *simp, SimplexElem search) {
    SimplexElem first = 0;
    SimplexElem last = simp->elementIndex;
    SimplexElem middle = (first + last) / 2;
    
    while (first <= last) {
        if (getElementAt(simp, middle) < search)
            first = middle + 1;
        else if (getElementAt(simp, middle) == search) {
            return true;
        }
        else
            last = middle - 1;
        
        middle = (first + last) / 2;
    }
    
    return false;
}
