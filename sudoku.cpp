/* ------------------------- PROBLEMA DEL SUDOKU ----------------------- */
#include <ga/GASimpleGA.h> //  Algoritmo Genetico simple
#include <ga/GA1DArrayGenome.h> // Genoma --> array de enteros (dim. 1) alelos
#include <iostream>
#include <fstream>
using namespace std;

// Definición de la plantilla del Sudoku.

struct plantilla{
       int tam;
       int *fijo;
};

int* marcadas;

void leerSudoku(struct plantilla *S,char *nombreF){
   ifstream f(nombreF);

   f>>S->tam;

   S->fijo = new int[S->tam*S->tam];

   for(int i=0;i<S->tam*S->tam;i++)
           f>>S->fijo[i];

   f.close();
}

// Inicialización de un individuo.

void InicioSudoku(GAGenome& g){

     GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;
     struct plantilla * plantilla1;
     plantilla1 = (struct plantilla *) genome.userData();

     int aux[plantilla1->tam];

     for(int f=0;f<plantilla1->tam;f++){

          for(int j=0;j<plantilla1->tam;j++) aux[j]=0;

          for(int j=1;j<=plantilla1->tam;j++){
            int v=GARandomInt(0,plantilla1->tam-1);
            while (aux[v]!=0) v=(v+1)%plantilla1->tam;
            aux[v]=j;
          }

          int i=0;

          while(i<plantilla1->tam){

              while((plantilla1->fijo[(f*plantilla1->tam)+i]==0) && (i<plantilla1->tam)) i++;

              if (i<plantilla1->tam){

                     bool encontrado=false;
                     for(int j=0;(j<plantilla1->tam) && (!encontrado);j++)
                             if (aux[j]==plantilla1->fijo[(f*plantilla1->tam)+i]) {
                                encontrado=true;
                                aux[j]=aux[i];
                             }

                     aux[i]=plantilla1->fijo[(f*plantilla1->tam)+i];
              }
              i++;

          }
          for(int c=0;c<plantilla1->tam;c++)
                  genome.gene((f*plantilla1->tam)+c,aux[c]);
     }
}

// Definición del método de cruce.

int CruceSudoku(const GAGenome& p1,const GAGenome & p2,GAGenome* c1,GAGenome* c2){

    GA1DArrayAlleleGenome<int> & m = (GA1DArrayAlleleGenome<int> &) p1;
    GA1DArrayAlleleGenome<int> & p = (GA1DArrayAlleleGenome<int> &) p2;
    struct plantilla * plantilla1 = (struct plantilla *) m.userData();
    int n=0;

    int punto1=GARandomInt(0,m.length());
    while ((punto1%plantilla1->tam)!=0) punto1++;
    int punto2=m.length()-punto1;

    if (c1){
             GA1DArrayGenome<int> & h1 = (GA1DArrayAlleleGenome<int> &) *c1;
             h1.copy(m,0,0,punto1);
             h1.copy(p,punto1,punto1,punto2);
             n++;
    }

    if (c2){
             GA1DArrayGenome<int> & h2 = (GA1DArrayAlleleGenome<int> &) *c2;
             h2.copy(p,0,0,punto1);
             h2.copy(m,punto1,punto1,punto2);
             n++;
    }

    return n;

}


bool checkColumna(int col[], int * check, int tam){
     bool repe=false;

     for(int i=0;i<tam;i++) check[i]=0;

     for(int i=0;i<tam;i++)
             check[col[i]-1]++;
     for(int i=0;i<tam;i++) if (check[i]>1) repe=true;

     return repe;
}

// Definición del método de mutación.

int MutacionSudoku(GAGenome& g,float pmut){

    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;
    struct plantilla * plantilla1;
    plantilla1 = (struct plantilla *) genome.userData();
    int nmut=0;
    int aux;
    int fil;
    bool fila;

    int caux[plantilla1->tam];
    int *checkC=new int[plantilla1->tam];

    if (pmut<=0.0) return 0;

    for(int f=0; f<plantilla1->tam; f++)
       for(int c=0; c<plantilla1->tam; c++)
          if (plantilla1->fijo[(f*plantilla1->tam)+c]==0){
           if (GAFlipCoin(pmut) ){
                if (GAFlipCoin(0.5)) fila = true;
                else fila = false;

                if (!fila){

                      for(int j=0;j<plantilla1->tam;j++) caux[j]=genome.gene((j*plantilla1->tam)+c);
                      if (checkColumna(caux,checkC,plantilla1->tam)){
                         int v1 = GARandomInt(0,plantilla1->tam-1);
                         while (checkC[v1]<=1) v1=(v1+1)%plantilla1->tam;
                         v1++;
                         int v2 = GARandomInt(0,plantilla1->tam-1);
                         while (checkC[v2]!=0) v2=(v2+1)%plantilla1->tam;
                         v2++;

                         bool encontrado = false;
                         for(int j=0;j<plantilla1->tam && !encontrado;j++)
                                 if ((plantilla1->fijo[j*(plantilla1->tam)+c]==0)&&(genome.gene(j*(plantilla1->tam)+c)==v1)){
                                    encontrado = true;
                                    genome.gene((j*plantilla1->tam)+c,v2);
                                    fil = j;
                                 }

                         int col=(c+1)%plantilla1->tam;
                         while(genome.gene((fil*plantilla1->tam)+col)!=v2) col=(col+1)%plantilla1->tam;
                         if (plantilla1->fijo[(fil*plantilla1->tam)+col]==0) {
                                nmut++;
                                genome.gene((fil*plantilla1->tam)+col,v1);
                         }
                         else {
                              genome.gene((fil*plantilla1->tam)+c,v1);
                         }

                      }

                }
                else{
                   int v1 = (c + 1) %plantilla1->tam;
                   while ((plantilla1->fijo[(f*plantilla1->tam)+v1]!=0)) v1=(v1+1)%plantilla1->tam;
                   aux = genome.gene((f*plantilla1->tam)+c);
                   genome.gene((f*plantilla1->tam)+c,genome.gene((f*plantilla1->tam)+v1));
                   genome.gene((f*plantilla1->tam)+v1,aux);
                   nmut++;
                }
           }
          }

    return nmut;
}

float Objective(GAGenome &); // Funcion objetivo --> al final

GABoolean Termina(GAGeneticAlgorithm &); // Funcion de terminacion --> al final

int main(int argc, char **argv)
{
// Declaramos variables para los parametros del GA y las inicializamos

    char* namef = argv[1];
    int popsize = atoi(argv[2]);
    int type_sel = atoi(argv[3]);
    float pcross = atof(argv[4]);
    float pmut = atof(argv[5]);

// Inicializamos la plantilla del sudoku y la matriz de posiciones marcadas

    struct plantilla* P = (struct plantilla*) malloc(sizeof(struct plantilla*));
    leerSudoku(P, namef);
    marcadas = new int[P->tam*P->tam];

// Conjunto enumerado de alelos --> valores posibles de cada gen del genoma

    GAAlleleSet<int> alelos;
    for(int i=0;i<=9;i++) alelos.add(i);

// Creamos el genoma y definimos operadores de inicio, cruce y mutación

    GA1DArrayAlleleGenome<int> genome(P->tam*P->tam,alelos,Objective,P);
    genome.initializer(InicioSudoku);
    genome.crossover(CruceSudoku);
    genome.mutator(MutacionSudoku);

// Creamos el algoritmo genetico

    GASimpleGA ga(genome);

// Inicializamos - minimizar funcion objetivo, tamaño poblacion, nº generaciones,
// pr. cruce y pr. mutacion, selección y le indicamos que evolucione.

    ga.minimaxi(-1);
    ga.populationSize(popsize);
    ga.nGenerations(12000);
    ga.pCrossover(pcross);
    ga.pMutation(pmut);
    switch (type_sel) {
        case '1' : {
            GARouletteWheelSelector selector;
            ga.selector(selector);
            }
        case '2' : {
            GATournamentSelector selector;
            ga.selector(selector);
            }
    }
    ga.terminator(Termina);
    ga.evolve(1);

// Imprimimos el mejor individuo que encuentra el GA y su valor fitness

    cout << ga.statistics().minEver() << "\n" << ga.statistics().bestIndividual() << endl;
    return 0;
}

// Funcion objetivo.

float Objective(GAGenome& g) {

    //Inicializamos las variables auxiliares que usaremos más tarde.
    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;
    for (int n = 0; n < 81; n++) marcadas[n] = 0;
    float repetidos = 0;
    int tamSud = 9;
    int aux[9];
    int marAux[9];
    int cont = 0;

    // Repetidos en la misma fila.
    for(int k=0; k<tamSud; k++)                                     // Para cada fila del sudoku comparamos cada elemento con los de su fila para ver
        for(int i=k*tamSud; i<(k+1)*tamSud; i++)                    // repeticiones. Saltamos en el array tamSud posiciones para cambiar de fila, ya que
           for(int j=i+1; j<(k+1)*tamSud; j++)                      // el array es lineal y no bidimensional.
                if ((genome.gene(i)==genome.gene(j)) && (marcadas[i] == 0)) {
                        repetidos++;                                // Cómo hemos encontrado una repetición la contamos y marcamos esa posición
                        marcadas[i] = 1;                            // como ya contada.
                }

    // Repetidos en la misma columna.
    for(int k=0; k<tamSud; k++)                                     // Para cada columna del sudoku comparamos cada elemento con los de su columna para ver
        for(int i=k; i<=tamSud*(tamSud-1)+k; i+=tamSud)             // repeticiones. Para recorrer las columnas avanzamos tamSud posiciones para cada elemento
           for(int j=i+tamSud; j<=tamSud*(tamSud-1)+k; j+=tamSud)   // del array lineal.
                if ((genome.gene(i)==genome.gene(j)) && (marcadas[i] == 0)) {
                        repetidos++;                                // Cómo hemos encontrado una repetición la contamos y marcamos esa posición
                        marcadas[i] = 1;                            // como ya contada.
                }

    // Repetidos en la misma subcuadrícula.
    for(int k=0; k<=6; k+=3) {                                      // 3 iteraciones para recorrer las tres columnas de la matriz de subcuadrículas donde
                                                                    // las primeras posiciones son 0, 3 y 6.
        for(int i=k; i<=k+54; i+=27) {                              // 3 iteraciones más para recorrer las 3 subcuadrículas que hay dentro de las columnas
                                                                    // anteriormente descritas.
           for(int j=i; j<=i+18; j+=9)                              // 3 iteraciones para recorrer las filas que hay dentro de cada subcuadrícula.
                for(int l=j; l<=j+2; l++) {                         // Últimas 3 iteraciones para recorrer los elementos de una fila.
                    aux[cont] = genome.gene(l);                     // Rellenamos un array auxiliar con los valores de la subcuadrícula y otro
                    marAux[cont] = marcadas[l];                     // con las posiciones marcadas de esa subcuadrícula.
                    cont++;
                }
           for(int a=0; a<9; a++)                                   // Recorremos el array con dos bucles y comprobamos si se repiten elementos tal y como
                for(int b=a+1; b<9; b++)                            // haríamos para una fila normal.
                   if ((aux[a] == aux[b]) && (marAux[a] == 0)) {
                        repetidos++;                                // Cómo hemos encontrado una repetición la contamos y marcamos esa posición
                        marAux[a] = 1;                              // como ya contada.
                   }
           cont = 0;
        }
    }
    return repetidos;
}

// Funcion de terminacion

GABoolean Termina(GAGeneticAlgorithm & ga){
    if ( (ga.statistics().minEver()==0) ||                          // La condición de parada es encontrar la solución óptima
        (ga.statistics().generation()==ga.nGenerations())) return gaTrue; // o llegar a 12000 generaciones.
    else return gaFalse;
}
