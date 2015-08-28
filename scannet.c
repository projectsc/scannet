/* 
 * File:   scannet.c
 * Author: austeclino
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

// Permitindo a sua declaração como um tipo qualquer:
typedef struct CN{
    int** adjacencyMatrix;
    int** neighborhoodMatrix;
    int** colorsMatrix;
    int** dendrogramMatrix;
    float** dendrogramGraphic;
    int numberOfCommunities;
    int qtdOfEdges;
    int numberRemovedEdges;
    int numberDendrogramLines;
    int* realIndex;
    int* realIndexDendrogramGraphic;
} complexNetwork;

typedef struct C_EDGE{
    int V1;
    int V2;
}edge;

// Prototypes
void loadSimilarityMatrix();
void createComplexNetworks(void);
void executeNga(complexNetwork);
void dijkstra (int, int**, int**, int);
void processSimilarityMatrix(int**);
void createAdjacencyMatrix(int, complexNetwork*);
int** createDynamicMatrix(int, int);
float** createDynamicMatrixFloat(int, int);
int* createDynamicArray(int);
void createNeighborhoodMatrix(complexNetwork);
void DestroyDynamicMatrix(int**, int );
void clearMatrix(int**);
void clearMatrixBySize(int**, int );
int isThereEdge(int**);
void calculateBetweenness(int, complexNetwork, float**);
void removeEdgesWithBiggerBetweenness(int[][2], float**, int**, int**);
void processDendrogram(int[][2],complexNetwork*, int**);
void processGraphicDendrogram(int[][2],complexNetwork*, int**);
void startDendrogram(complexNetwork*);
void copyMatrixOfFloat(int**, float**);
void copyMatrix(int**, int**);
void writeLog(int*, int*);
void conveterArrayInMatrix(int*, int**, int, int);
void conveterMatrixInArray(int**, int*, int , int );
int findInArray(int*, int, int);
void sortCommunitiesList(complexNetwork*, int);
void sortGraphicDendrograma(complexNetwork*);
void createColorsMatrix(complexNetwork*);
void setupEixoYDendrogram(complexNetwork*, int);
void saveMatrizInFile();
void saveMatrizInFileFloat();
int findRealIndex(int* realIndex, int);
void adjustNeighborhoodMatrix(int[][2], int**, int**, int);
void plotGraphic();
double getTime();
complexNetwork createNetworkTeste();
// Global vars
// Criando a enumeração:
enum boolean {
    true = 1, false = 0
};
// Permitindo a sua declaração como um tipo qualquer:
typedef enum boolean bool;

int sizeSimilarityMatrix = 1;
int** similarityMatrix;
int NCutPoints = 101;
int numThreads = 4;
int isInCluster = 1;
int N_ARESTAS = 100;
int indexMelhorRede = 0;
int qtdDendrogramLines = 0;
complexNetwork *complexNetworks;
complexNetwork SelectBetterNetworkToFindCommunities();
edge *removedEdges;
int calculateDiameter();
void showtime();
float calculateDistance();
void worker();
int qtdNodeWithEdgeRemoved =0;
char current_path_entrada[1000];
char current_path_saida[1000];
int rank;
int size;

//Variáveis de Tempo
double start_t_RC, end_t_RC, tempo_total_RC, start_t_Distancia, end_t_Distancia, 
       tempo_total_Distancia, start_t_NGA, end_t_NGA, tempo_total_NGA, start_t_MC, 
       end_t_MC, tempo_total_MC, start_t_TOTAL, end_t_TOTAL, tempo_total_TOTAL ;

double start_t_Paralelo, end_t_Paralelo, total_t_Paralelo;
double * tempos_paralelo;
FILE *arqLog;


int main(int argc, char** argv) {
        complexNetwork betterNetwork;
        showtime("Inicio da execucao");
        
         printf("\n\n Informe a quantidade de vertices da rede analisada:  ");
         scanf("%d", &sizeSimilarityMatrix);

        // -------------- PARAMETROS DE CONFIGURAÇÃO ---------------
       //  printf("\n\n Informe o diret�rio onde se encontra o arquivo: (Ex: C:/Users/Leitura)  ");
      //   scanf("%s", current_path_entrada);

         strcpy(current_path_entrada,"C:/Users/Austeclino/Desktop/Mestrado/leitura/");
        // ---------------------------------------------------------
        
        similarityMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
        loadSimilarityMatrix("resultBlast.txt", similarityMatrix);   
        start_t_TOTAL = getTime();
      	start_t_RC = getTime();      
        //      ---------------- CRIAÇÃO DAS REDES COMPLEXAS ---------------------------
        createComplexNetworks();
        end_t_RC = getTime();
        tempo_total_RC = end_t_RC - start_t_RC;

        //      ---------------- SELEÇÃO DA REDE CR�?TICA -------------------------------
        
        start_t_Distancia = getTime();
        betterNetwork = SelectBetterNetworkToFindCommunities();
        end_t_Distancia = getTime();
        tempo_total_Distancia = end_t_Distancia - start_t_Distancia;
        printf("\n Quantidade Total de Arestas: %d", betterNetwork.qtdOfEdges);
        
        //      ---------------- PREPARAÇÃO DO NGA -------------------------------
        betterNetwork.colorsMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);      
        betterNetwork.dendrogramMatrix = createDynamicMatrix(betterNetwork.qtdOfEdges+2, sizeSimilarityMatrix);
        //##
        betterNetwork.dendrogramGraphic = createDynamicMatrixFloat(betterNetwork.qtdOfEdges*2, sizeSimilarityMatrix+1);
        betterNetwork.realIndex = createDynamicArray(sizeSimilarityMatrix);
        betterNetwork.realIndexDendrogramGraphic = createDynamicArray(sizeSimilarityMatrix);

        int i = 0;
        for (i = 0; i < sizeSimilarityMatrix; i++){ 
            betterNetwork.realIndex[i] = i;
            betterNetwork.realIndexDendrogramGraphic[i] = i;
        }
        copyMatrix(betterNetwork.neighborhoodMatrix, betterNetwork.colorsMatrix);

        //      ---------------- EXECUÇÃO DO NGA -------------------------------
        executeNga(betterNetwork);
        //      ---------------- GERAÇÃO DAS SA�?DAS -------------------------------
  	end_t_TOTAL = getTime();
        tempo_total_TOTAL = end_t_TOTAL - start_t_TOTAL;
        
        for (i = 0; i < sizeSimilarityMatrix; i++){ 
           betterNetwork.dendrogramMatrix[betterNetwork.qtdOfEdges+1][i] = betterNetwork.realIndex[i]+1;
           betterNetwork.dendrogramGraphic[qtdDendrogramLines+1][i] = betterNetwork.realIndexDendrogramGraphic[i]+1;     
        }
        
        setupEixoYDendrogram(&betterNetwork, qtdDendrogramLines);
        
        saveMatrizInFile(betterNetwork.colorsMatrix, "mColors.dat", sizeSimilarityMatrix, sizeSimilarityMatrix );
        saveMatrizInFile(betterNetwork.dendrogramMatrix, "mDendrogram.dat",betterNetwork.qtdOfEdges+2, sizeSimilarityMatrix);
        saveMatrizInFileFloat(betterNetwork.dendrogramGraphic, "mDendrogramGraphic.dat",qtdDendrogramLines+2, sizeSimilarityMatrix+1);
  
        showtime("Fim da execucao");  
        free(similarityMatrix);
        free(complexNetworks);
        
        char path[1000];
        strcpy(path,current_path_entrada);
        strcat(path,"logexec.txt");
        FILE *arqLog = fopen(path, "wt");
        fprintf(arqLog,"\n\n ------- Resumo dos tempos de execucao(ms) --------- \n");
        fprintf(arqLog,"\n Quantidade de Arestas: %d", betterNetwork.qtdOfEdges);
        fprintf(arqLog,"\n Tempo gasto na criacao das 101 redes complexas: %.5lf", tempo_total_RC); 
        fprintf(arqLog,"\n Tempo gasto na selecao da rede complexa ideal: %.5lf", tempo_total_Distancia); 
        fprintf(arqLog,"\n A rede ideal selecionada possui %d arestas", betterNetwork.qtdOfEdges);
        fprintf(arqLog,"\n Tempo gasto na NGA: %.5lf", tempo_total_NGA);
        fprintf(arqLog,"\n Tempo gasto na execucao da matriz de cores: %.5lf", tempo_total_MC);
        fprintf(arqLog,"\n Tempo total de execucao: %.5lf \n", tempo_total_TOTAL);
        fprintf(arqLog,"\n --------------------------------------------------- \n\n\n");
        fclose(arqLog); 
        system("pause"); 
        return 0;
}


// Traz a �ltima coluna(eixo Y) para a primeira coluna, pois o origim considera a primeira 
// coluna como sendo o eixo Y para plotar o dendrograma. 
void setupEixoYDendrogram(complexNetwork *complexNetwork, int numberDendrogramLines){         
    int i,j;    
    for (i = 0; i < numberDendrogramLines; i++){ 
        float value = (*complexNetwork).dendrogramGraphic[i][sizeSimilarityMatrix];
        for(j = sizeSimilarityMatrix; j > 0; j--){
            (*complexNetwork).dendrogramGraphic[i][j] = (*complexNetwork).dendrogramGraphic[i][j-1];
        } 
        (*complexNetwork).dendrogramGraphic[i][0] = value;
    }
}

void conveterArrayInMatrix(int* array, int** matrix, int n, int m){
    int i = 0;  
    int j = 0;
    int index = 0;
    int size = n*m;
    for (index = 0; index < size; index++){ 
        matrix[i][j] = array[index];
        
        if(j == (m-1) ){
            j = 0;
            i++;
        }else{
            j++;
        } 
    }
}

void conveterMatrixInArray(int** matrix, int* array, int n, int m){
    int i = 0;  
    int j = 0;
    int index = 0;
    
    for (i = 0; i < n; i++){ 
        for (j = 0; j < m; j++){ 
            array[index] = matrix[i][j];
            index++;
        }
    }
}

void showtime(char texto[30]){
    time_t currentTime;
    struct tm *timeinfo;
    currentTime = time(NULL);
    timeinfo = localtime(&currentTime);
    printf("\n %s: %02d:%02d:%02d \n",texto, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);   
}

void showtimeWithParameter(char texto[30], int value){
    time_t currentTime;
    struct tm *timeinfo;
    currentTime = time(NULL);
    timeinfo = localtime(&currentTime);
    printf("\n%s: %d  Hora: %02d:%02d:%02d \n",texto, value, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);   
}

double getTime(){
    time_t t;
    t = time(NULL);
    return (double)t;
}

int** createDynamicMatrix(int sizeLine, int sizeColumn){
    int** matrix = malloc(sizeLine * sizeof(int*));    
    int i;
    
    if(matrix == NULL){
         printf("\n ERRO NA ALOCAÇÃO DE MEMÓRIA. \n");
         exit(1);
    }
        
    for(i=0; i < sizeLine; i++){
        matrix[i] = calloc(sizeColumn, sizeof(int));
        if(matrix[i] == NULL){
            printf("\n ERRO NA ALOCAÇÃO DE MEMÓRIA. \n");
            exit(1);
        }
    }
    
    return matrix;
}

int* createDynamicArray(int size){
    int* matrix = calloc(size, sizeof(int));     
  
    if(matrix == NULL){
         printf("\n ERRO NA ALOCAÇÃO DE MEMÓRIA. \n");
         exit(1);
    }
    
    return matrix;
}

void reiniciaMatriz(int sizeLine, int sizeColumn, int** matrix){
    int i, j;
    
    for(i=0; i < sizeLine; i++){
        for(j=0; j < sizeColumn; j++){
            matrix[i][j] = 0;
        }
    }
}

void createAdAndNM(int sizeLine, int sizeColumn, int** adjacencyMatrix, int ** neighborhoodMatrix){
    int i;
    
    for(i=0; i < sizeLine; i++){
        adjacencyMatrix[i] = calloc(sizeColumn, sizeof(int));
        neighborhoodMatrix[i] = calloc(sizeColumn, sizeof(int));
    }
}


float** createDynamicMatrixFloat(int sizeLine, int sizeColumn){
    float ** matrix = malloc(sizeLine * sizeof(float*));    
    int i;
    
    if(matrix == NULL){
         printf("\n ERRO NA ALOCAÇÃO DE MEMÓRIA. \n");
         exit(1);
    }
    
   // #pragma omp parallel for private(i) num_threads(numThreads)
    for(i=0; i < sizeLine; i++){
        matrix[i] = calloc(sizeColumn, sizeof(float));
         if(matrix[i] == NULL){
            printf("\n ERRO NA ALOCAÇÃO DE MEMÓRIA. \n");
            exit(1);
         }
    }
    
    return matrix;
}

void DestroyDynamicMatrix(int** matrix, int sizeLine){
    int i;
    
    for(i=0; i < sizeLine; i++){
        free(matrix[i]); 
    }
    free(matrix);
}

void DestroyDynamicMatrixFloat(float** matrix, int sizeLine){
    int i;
    
    for(i=0; i < sizeLine; i++){
        free(matrix[i]); 
    }
    free(matrix);
}

void loadSimilarityMatrix(char nameOfFile[30], int** matrix){
    FILE *arqBLAST; 
    char *endLinha; 
    int size = sizeSimilarityMatrix * 4; // Estimativa do tamanho da linha do arquivo.
    char linha[size];

    char *token = NULL;
    int n = 0;
    int indexLine = 0;
    int indexCollumn = 0;
  
    char path[1000];
    strcpy(path,current_path_entrada);  
    strcat(path,nameOfFile);

    arqBLAST = fopen(path, "rt");
    
    if (arqBLAST == NULL){
        printf("\n Arquivo de resultados do BLAST Nao pode ser aberto\n");
        exit (EXIT_FAILURE);
    }

    while (indexLine < sizeSimilarityMatrix) {
        endLinha = fgets(linha, size, arqBLAST);
        if (endLinha){ 
            token= strtok(linha, " " );
            n = atol(token);
            matrix[indexLine][indexCollumn] = n;
           
            
            indexCollumn++;
            
            while (indexCollumn < sizeSimilarityMatrix) {
                 token= strtok(NULL, " " ); 
                 n = atol(token);
                 matrix[indexLine][indexCollumn] = n;
                 indexCollumn++;
            }             
        }
        indexCollumn = 0;
        indexLine++;
    }
         
   fclose(arqBLAST); // Fechamento do arquivo.
} 


void saveMatrizInFile(int** matriz, char nameOfFile[30], int n, int m){         
    FILE *arqMATRIZ;
    int w,k,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
   
    if (arqMATRIZ == NULL){ // Se não conseguiu criar
        printf("Problemas na criacao do arquivo da matriz de similaridade.\n");
        exit (EXIT_FAILURE);
    }

    for (w = 0; w < n; w++){
        for (k = 0; k < m; k++){
            if(matriz[w][k] < 10){
                result = fprintf(arqMATRIZ,"%d   ", matriz[w][k]);
            }else{
                result = fprintf(arqMATRIZ,"%d  ", matriz[w][k]);
            }
            if (result == EOF) printf("Erro na gravacao do arquivo\n");
        }

        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("Erro na gravacao do arquivo da matriz de similaridade\n");
    }

    fclose(arqMATRIZ);

    printf("\n O arquivo %s foi salvo com sucesso! \n", nameOfFile);
}
    
void saveMatrizInFileFloat(float** matriz, char nameOfFile[30], int n, int m){         
    FILE *arqMATRIZ;
    int w,k,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
    
    if (arqMATRIZ == NULL){ // Se não conseguiu criar
        printf("Problemas na criacao do arquivo da matriz de similaridade.\n");
        exit (EXIT_FAILURE);
    } 

    for (w = 0; w < n; w++){
        for (k = 0; k < m; k++){
            result = fprintf(arqMATRIZ,"%.2f ", matriz[w][k]);
            if (result == EOF) printf("Erro na gravacao do arquivo\n");
        }

        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("Erro na gravacao do arquivo da matriz de similaridade\n");
    }

    fclose(arqMATRIZ);

    printf("\n O arquivo %s foi salvo com sucesso! \n", nameOfFile);
}

void saveArrayInFileFloat(float* array, char nameOfFile[30], int size){         
    FILE *arqMATRIZ;
    int w,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
   
    if (arqMATRIZ == NULL){ // Se não conseguiu criar
        printf("Problemas na criacao do arquivo de dist�ncias\n");
        exit (EXIT_FAILURE);
    }

    for (w = 0; w < size; w++){
        result = fprintf(arqMATRIZ,"%.2f", array[w]);
        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("Erro na gravacao do arquivo\n");
    }

    fclose(arqMATRIZ);

    printf("\n O arquivo %s foi salvo com sucesso! \n", nameOfFile);
}

void createComplexNetworks(){
    complexNetworks = (complexNetwork *)malloc(101 * sizeof(complexNetwork));
    if(complexNetworks == NULL){
         printf("\n ERRO NA ALOCAÇÃO DE MEMÓRIA. \n");
         exit(1);
    }

    int i;
    
    for (i = 0;  i < NCutPoints; i++){
        complexNetworks[i].adjacencyMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
        complexNetworks[i].neighborhoodMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
        complexNetworks[i].numberOfCommunities = 1;
        complexNetworks[i].numberRemovedEdges = 0;
        complexNetworks[i].numberDendrogramLines = 0;
        complexNetworks[i].qtdOfEdges = 0;
                
        createAdjacencyMatrix(i, complexNetworks);
        createNeighborhoodMatrix(complexNetworks[i]);
    }
}   

void createAdjacencyMatrix(int cutPoint, complexNetwork *complexNetworks){
    int i,j;    
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for(j = (i+1); j < sizeSimilarityMatrix; j++){
            if(similarityMatrix[i][j] >= cutPoint){
                complexNetworks[cutPoint].adjacencyMatrix[i][j] = 1;
                complexNetworks[cutPoint].adjacencyMatrix[j][i] = 1;
                complexNetworks[cutPoint].qtdOfEdges++;
            }
        } 
    }
}

void createNeighborhoodMatrix(complexNetwork complexNetwork){
    int i;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){
         dijkstra(i, complexNetwork.adjacencyMatrix, complexNetwork.neighborhoodMatrix, sizeSimilarityMatrix);
    }
}

// Calcula as distâncias de 'Vi' a todos os outros vértices de um grafo com 'V' vértices e armazena-as em dis[]
void dijkstra (int Vi, int** adjacencyMatrix, int** neighborhoodMatrix, int maxSize)
{       
    // Armazena a distância mínima partindo de um vértice 'i' até todos os outros vértices
    // dis[j] representa a menor distância de 'i' a 'j'.
    int dis[maxSize];

    // vis[i] informa se o vértice 'i' já foi visitado/analisado ou não (inicialmente nenhum vértice foi)
    char vis[maxSize];
    memset (vis, 0, sizeof (vis));

    // Inicialmente afirmamos que a menor distância encontrada entre Vi e qualquer outro vértice (exceto o próprio Vi) é infinita
    memset (dis, 0x7f, sizeof (dis));
    dis[Vi] = 0;

    while (1){
        int i, n = -1;
        for (i = 0; i < maxSize; i++){
            if (!vis[i] && (n < 0 || dis[i] < dis[n])){
                    n = i;
            }        
        }

        if (n < 0){
           break;
        }

        vis[n] = 1;
        for (i = 0; i < maxSize; i++){
            if (adjacencyMatrix[n][i] && dis[i] > dis[n] + adjacencyMatrix[n][i]){
                dis[i] = dis[n] + adjacencyMatrix[n][i];
                neighborhoodMatrix[Vi][i] = dis[i];
                neighborhoodMatrix[i][Vi] = dis[i];
            }        
        }        
    }
}

complexNetwork SelectBetterNetworkToFindCommunities(){
      float distance = 0;
      float biggerDistance = 0;
      int i = 0;
      complexNetwork betterNetwork;
      float distances[101];
      int i_distance = 0;
      // ESSE TRECHO DEVE SER DESCOMENTADO CASO QUEIRA SE FAZER UMA AN�LISE DA REDE CR�?TICA    
      for (i = 0; i < (NCutPoints-1); i++){ 
          distance =  calculateDistance(complexNetworks[i], complexNetworks[i+1]);    
          printf("\n Distancia entre a rede %d e %d = %f", i, i+1,distance);
          distances[i_distance] = distance;
          i_distance = i_distance + 1;
          if(distance > biggerDistance){
              biggerDistance = distance;
              indexMelhorRede = i;
          }
      }
      saveArrayInFileFloat(distances, "mDistanceGraphic.dat",i_distance);
      
      plotGraphic();
      
      printf("\n\n\n A rede %d possui a maior distancia: %f  ", indexMelhorRede,  biggerDistance);
      
      printf("\n\n\n Escolha o pico que deseja utilizar:  ");
      scanf("%d", &indexMelhorRede);
      betterNetwork = complexNetworks[indexMelhorRede];
      return betterNetwork;
}


void plotGraphic(){
    char path[1024];
    char syscommand[1024];
    strcpy(path,current_path_entrada);
    strcat(path,"mDistanceGraphic.dat");
    sprintf(syscommand, "start gnuplot -persist -e \"plot '%s' with lines title 'SCANNET - Gr�fico de Dist�ncias';pause -1 'Pressione Enter para sair. . .';\"\"",path);
    system(syscommand);

}

float calculateDistance(complexNetwork networkOne, complexNetwork networkTwo){
    int i = 0;
    int j = 0;
    float t1 = 1;
    float t2 = (sizeSimilarityMatrix * sizeSimilarityMatrix);
    float a = t1/t2;
    float r = 0;
    float total = 0;
     
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = (i+1); j < sizeSimilarityMatrix; j++){ 
            r = (networkOne.neighborhoodMatrix[i][j]) - (networkTwo.neighborhoodMatrix[i][j]);
            r = r * r;
            total = total + r;        
        }
     }
    
    return a * total;
}

int calculateDiameter(complexNetwork network){
    int i = 0;
    int j = 0;
    int bigger = 0;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
         for (j = (i+1); j < sizeSimilarityMatrix; j++){ 
             if(network.neighborhoodMatrix[i][j] > bigger){
                 bigger = network.neighborhoodMatrix[i][j];
             }   
         }
    }

    return bigger;
}

void executeNga(complexNetwork complexNetwork){
    start_t_NGA = getTime();
     
    float** betweennessMatrix = createDynamicMatrixFloat(sizeSimilarityMatrix, sizeSimilarityMatrix);
    int** neighborhoodMatrixPrevious = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
    int diameter = 0, i = 1;
    int parcialArestasRemovidas = 0, totalArestasRemovidas = 0; 
    startDendrogram(&complexNetwork);
    while( isThereEdge(complexNetwork.adjacencyMatrix) ){
        diameter = calculateDiameter(complexNetwork);
        calculateBetweenness(diameter, complexNetwork, betweennessMatrix); 
        int nodeWithEdgeRemoved[i][2];
        qtdNodeWithEdgeRemoved = 0;
        copyMatrix(complexNetwork.neighborhoodMatrix ,neighborhoodMatrixPrevious);
        removeEdgesWithBiggerBetweenness(nodeWithEdgeRemoved, betweennessMatrix, complexNetwork.adjacencyMatrix, complexNetwork.neighborhoodMatrix);  
        complexNetwork.numberRemovedEdges++;  
        adjustNeighborhoodMatrix(nodeWithEdgeRemoved, complexNetwork.adjacencyMatrix, complexNetwork.neighborhoodMatrix, diameter);
        processDendrogram(nodeWithEdgeRemoved, &complexNetwork, neighborhoodMatrixPrevious);  
        processGraphicDendrogram(nodeWithEdgeRemoved, &complexNetwork, neighborhoodMatrixPrevious);  
        writeLog(&parcialArestasRemovidas, &totalArestasRemovidas);
    }

    end_t_NGA = getTime();
    tempo_total_NGA = end_t_NGA - start_t_NGA;

    start_t_MC = getTime();
    createColorsMatrix(&complexNetwork);
    end_t_MC = getTime();
    tempo_total_MC = end_t_MC - start_t_MC;
    
    DestroyDynamicMatrixFloat(betweennessMatrix, sizeSimilarityMatrix);
    DestroyDynamicMatrix(neighborhoodMatrixPrevious, sizeSimilarityMatrix);
    
    qtdDendrogramLines = complexNetwork.numberDendrogramLines;
}

void startDendrogram(complexNetwork *complexNetwork){
    int i,j, index = 0;
    int index2 = 0;
    int* newCommunityItens = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix); 
    int* newCommunityItens2 = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);  
    int pontoMaximo = 0;
    float dividendo = 2.0;
    int realIndex = 0;
    int hasNewCommunity = 0;
    int repeatLine = 1;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){
        (*complexNetwork).dendrogramGraphic[0][i] = (sizeSimilarityMatrix+1) / dividendo;
    }
    (*complexNetwork).numberDendrogramLines = 1;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){
        if(findInArray(newCommunityItens2, i, index2) == -1){
            index = 0;
            newCommunityItens[index] = i;
            newCommunityItens2[index2] = i;
            realIndex = findRealIndex((*complexNetwork).realIndex, i);
            (*complexNetwork).dendrogramMatrix[0][realIndex] = (*complexNetwork).numberOfCommunities;
            index2++;
            index++;
                           
             for (j = 0; j < sizeSimilarityMatrix; j++){
                 if( (*complexNetwork).neighborhoodMatrix[i][j] != 0 ){
                      hasNewCommunity = 1;
                      if( (findInArray(newCommunityItens, j, index) == -1) ){
                          realIndex = findRealIndex((*complexNetwork).realIndex, j);
                          (*complexNetwork).dendrogramMatrix[0][realIndex] = (*complexNetwork).numberOfCommunities;
                           newCommunityItens[index] = j;
                           index++;
                           newCommunityItens2[index2] = j;
                           index2++;
                      }
                 }
            }
            
            if(hasNewCommunity && repeatLine){
                repeatLine = 0;
                //##
                for (j = 0; j < sizeSimilarityMatrix+1; j++){
                    (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][j] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines - 1][j];
                }
                (*complexNetwork).numberDendrogramLines = (*complexNetwork).numberDendrogramLines + 1;
            }
            
            for (j = 0; j < index; j++){
               realIndex = findRealIndex((*complexNetwork).realIndexDendrogramGraphic, newCommunityItens[j]);
               (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][realIndex] = pontoMaximo + ((index+1) / dividendo);  
            }
            //##
            (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][sizeSimilarityMatrix] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1][sizeSimilarityMatrix] + 1;
            
            sortCommunitiesList(complexNetwork, (*complexNetwork).numberOfCommunities); 
            sortGraphicDendrograma(complexNetwork);
            
            pontoMaximo = pontoMaximo + index;
 
            (*complexNetwork).numberOfCommunities++;         
        } 
    }
    (*complexNetwork).numberDendrogramLines =  (*complexNetwork).numberDendrogramLines + 1;
    
    (*complexNetwork).numberOfCommunities--;
    free(newCommunityItens);
    free(newCommunityItens2);
}

int findRealIndex(int* realIndex, int virtualIndex){
    int j;
    for (j = 0; j < sizeSimilarityMatrix; j++){
        if(realIndex[j] == virtualIndex ){
            return j;
        }
    }  
    return -1;
}


void processGraphicDendrogram(int nodeWithEdgeRemoved[1][2], complexNetwork *complexNetwork, int** neighborhoodMatrixPrevious){
    int i,j, index = 0, index2 = 0, index3 = 0;
    int* newCommunityItens = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);
    int* itensOldCommunity = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);  
    int* AuxCommunity = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);  
    
    int hasNewCommunity = 0;
    int finded1 = 0; 
    int finded2 = 0; 
    
    int indiceJ = nodeWithEdgeRemoved[0][0];
    int indiceI = nodeWithEdgeRemoved[0][1];
    
    float dividendo = 2.0;
    int realIndex = 0;

    if( (*complexNetwork).neighborhoodMatrix[indiceI][indiceJ] == 0) {
        hasNewCommunity = 1;
        newCommunityItens[index] = indiceJ;
        index++;
        itensOldCommunity[index2] = indiceI;
        index2++; 
        for (j = 0; j < sizeSimilarityMatrix; j++){
            if(((*complexNetwork).neighborhoodMatrix[indiceI][j] == 0) && (neighborhoodMatrixPrevious[indiceI][j] != 0) ){
                 if( (findInArray(newCommunityItens, j, index) == -1)  && ((*complexNetwork).neighborhoodMatrix[indiceJ][j] != 0) ){
                     newCommunityItens[index] = j;
                     index++;
                 }
            }else if ((*complexNetwork).neighborhoodMatrix[indiceI][j] != 0) {
                if (findInArray(itensOldCommunity, j, index2) == -1) {
                    itensOldCommunity[index2] = j;
                    index2++;
                }
            }
        } 
    }
    
    // ADICIONADO PARA MANTER AS COMUNIDADES COM MENOR N�MERO SEMPRE NA ESQUERDA NO DENDROGRAMA
    if(index2 > index){
        for (j = 0; j < index2; j++){
            AuxCommunity[j] = itensOldCommunity[j];
        }
        
        for (j = 0; j < index; j++){
            itensOldCommunity[j] = newCommunityItens[j];
        }
        
        for (j = 0; j < index2; j++){
            newCommunityItens[j] = AuxCommunity[j];
        } 
        index3 = index2;
        index2 = index;
        index = index3;
    }
    
    //##
    for (j = 0; j < sizeSimilarityMatrix; j++){
        (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][j] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1][j];             
    }
    (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][sizeSimilarityMatrix] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1][sizeSimilarityMatrix] ;             
   
    (*complexNetwork).numberDendrogramLines = (*complexNetwork).numberDendrogramLines+1;
    i = (*complexNetwork).numberDendrogramLines;
        
    if(hasNewCommunity == 1){
        float positionBefore = 0.0;
        (*complexNetwork).numberOfCommunities = (*complexNetwork).numberOfCommunities + 1;
        
        for (j = 0; j < sizeSimilarityMatrix; j++){
            finded1 = findInArray(newCommunityItens, j, index);
            finded2 = findInArray(itensOldCommunity, j, index2);
            realIndex = findRealIndex((*complexNetwork).realIndexDendrogramGraphic, j);
            positionBefore = (*complexNetwork).dendrogramGraphic[i-1][realIndex];
            
            if(finded2 != -1){    
                (*complexNetwork).dendrogramGraphic[i][realIndex] = positionBefore - (index/dividendo);
            }else if(finded1 != -1){
                (*complexNetwork).dendrogramGraphic[i][realIndex] = positionBefore + (index2/dividendo);
            }else{
                (*complexNetwork).dendrogramGraphic[i][realIndex] = (*complexNetwork).dendrogramGraphic[i-1][realIndex];             
            }                    
        }
   
        sortGraphicDendrograma(complexNetwork);
        //##
        (*complexNetwork).dendrogramGraphic[i][sizeSimilarityMatrix] = (*complexNetwork).dendrogramGraphic[i-1][sizeSimilarityMatrix] + 1;             
        (*complexNetwork).numberDendrogramLines = (*complexNetwork).numberDendrogramLines+1;
    }  
            
    free(newCommunityItens);
    free(AuxCommunity);
    free(itensOldCommunity);   
}



void processDendrogram(int nodeWithEdgeRemoved[1][2], complexNetwork *complexNetwork, int** neighborhoodMatrixPrevious){
    int i,j, index = 0, index2 = 0;
    int* newCommunityItens = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);
    int* itensOldCommunity = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);   
    int hasNewCommunity = 0;
    int finded = 0; 
    int indiceJ = nodeWithEdgeRemoved[0][0];
    int indiceI = nodeWithEdgeRemoved[0][1];

    if( (*complexNetwork).neighborhoodMatrix[indiceI][indiceJ] == 0) {
        hasNewCommunity = 1;
        newCommunityItens[index] = indiceJ;
        index++;
        for (j = 0; j < sizeSimilarityMatrix; j++){
            if(((*complexNetwork).neighborhoodMatrix[indiceI][j] == 0) && (neighborhoodMatrixPrevious[indiceI][j] != 0) ){
                 if( (findInArray(newCommunityItens, j, index) == -1)  && ((*complexNetwork).neighborhoodMatrix[indiceJ][j] != 0) ){
                      newCommunityItens[index] = j;
                      index++;
                 }
            }else if ((*complexNetwork).neighborhoodMatrix[indiceI][j] != 0) {
                if (findInArray(itensOldCommunity, j, index2) == -1) {
                    itensOldCommunity[index2] = j;
                    index2++;
                }
            }
        } 
    }
    
    // ADICIONADO PARA MANTER AS COMUNIDADES COM MENOR N�MERO SEMPRE NA ESQUERDA NO DENDROGRAMA
    if(index2 > index){
        for (j = 0; j < index2; j++){
            newCommunityItens[j] = itensOldCommunity[j];
        } 
        index = index2;
    }
    
    i = (*complexNetwork).numberRemovedEdges;
    if(hasNewCommunity == 1){
        (*complexNetwork).numberOfCommunities = (*complexNetwork).numberOfCommunities + 1;
        int item = newCommunityItens[0];
        int realIndex = findRealIndex((*complexNetwork).realIndex, item);
        int level = (*complexNetwork).dendrogramMatrix[i-1][realIndex] + 1;
       
        for (j = 0; j < sizeSimilarityMatrix; j++){
            finded = findInArray(newCommunityItens, j, index);
            realIndex = findRealIndex((*complexNetwork).realIndex, j);
            if(finded == -1){
                if((*complexNetwork).dendrogramMatrix[i-1][realIndex] >= level){
                    (*complexNetwork).dendrogramMatrix[i][realIndex] = (*complexNetwork).dendrogramMatrix[i-1][realIndex] + 1;
                }else{
                    (*complexNetwork).dendrogramMatrix[i][realIndex] = (*complexNetwork).dendrogramMatrix[i-1][realIndex];
                }                
            }else{
                level = (*complexNetwork).dendrogramMatrix[i-1][realIndex] + 1; 
                (*complexNetwork).dendrogramMatrix[i][realIndex] = level;             
            } 
        }
        sortCommunitiesList(complexNetwork, level);
    } else{
        for (j = 0; j < sizeSimilarityMatrix; j++){
            int realIndex = findRealIndex((*complexNetwork).realIndex, j);
           (*complexNetwork).dendrogramMatrix[i][realIndex] = (*complexNetwork).dendrogramMatrix[i-1][realIndex];
        }
    }  
    
    free(itensOldCommunity); 
    free(newCommunityItens);
}

void sortCommunitiesList(complexNetwork *complexNetwork, int level){
    int i = (*complexNetwork).numberRemovedEdges;
    int j, l, m, aux, realIndexL, realIndexM = 0;
     
    for (j = (sizeSimilarityMatrix-1); j > -1; j--){
       if((*complexNetwork).dendrogramMatrix[i][j] == level){
            m = j; 
            for (l = (j+1); l < sizeSimilarityMatrix; l++){ 
                 if((*complexNetwork).dendrogramMatrix[i][l] < (*complexNetwork).dendrogramMatrix[i][m]){
                     aux = (*complexNetwork).dendrogramMatrix[i][m];
                     (*complexNetwork).dendrogramMatrix[i][m] = (*complexNetwork).dendrogramMatrix[i][l];
                     (*complexNetwork).dendrogramMatrix[i][l] = aux;
                     aux = (*complexNetwork).realIndex[m];
                     (*complexNetwork).realIndex[m] =  (*complexNetwork).realIndex[l];
                     (*complexNetwork).realIndex[l] = aux;
                     realIndexL = (*complexNetwork).realIndex[l];                          
                     realIndexM = (*complexNetwork).realIndex[m];
                     m++;
                 }
             }
        } 
    }
}



void sortGraphicDendrograma(complexNetwork *complexNetwork){
    int i = (*complexNetwork).numberDendrogramLines;
    int j, l, m, auxIndice, realIndexL, realIndexM = 0;
    float aux = 0.0;
    
    for (j = (sizeSimilarityMatrix-1); j > -1; j--){
        m = j; 
        for (l = (j+1); l < sizeSimilarityMatrix; l++){ 
             if((*complexNetwork).dendrogramGraphic[i][l] < (*complexNetwork).dendrogramGraphic[i][m]){
                 aux = (*complexNetwork).dendrogramGraphic[i][m];
                 (*complexNetwork).dendrogramGraphic[i][m] = (*complexNetwork).dendrogramGraphic[i][l];
                 (*complexNetwork).dendrogramGraphic[i][l] = aux;
                 auxIndice = (*complexNetwork).realIndexDendrogramGraphic[m];
                 (*complexNetwork).realIndexDendrogramGraphic[m] =  (*complexNetwork).realIndexDendrogramGraphic[l];
                 (*complexNetwork).realIndexDendrogramGraphic[l] = auxIndice;
                 realIndexL = (*complexNetwork).realIndexDendrogramGraphic[l];                          
                 realIndexM = (*complexNetwork).realIndexDendrogramGraphic[m];
                 m++;
             }
         }
    }
}



// Esse método é responsável pelas trocas de linhas e colunas das matrizes para 
// geração da metriz de cores e dendrograma.
void createColorsMatrix(complexNetwork *complexNetwork){
    int index = 0, index1, index2 = 0, i;
    int value = 0;
    int * arrayToSort = createDynamicArray(sizeSimilarityMatrix);
  
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        arrayToSort[i] = i;
    } 

    for (index = 0; index < sizeSimilarityMatrix; index++){
        index1 =(*complexNetwork).realIndex[index];
        if(index1 != index){       
            
            // Procura onde está a coluna real referente a coluna index1. 
            // (considere que trocas estão sendo feitas por isso essa necessidade)
            index1 = findRealIndex(arrayToSort, index1);
            value = arrayToSort[index1];
            arrayToSort[index1] = arrayToSort[index];
            arrayToSort[index] = value;
            
            // troca colunas
            for (index2 = 0; index2 < sizeSimilarityMatrix; index2++){              
                value = (*complexNetwork).colorsMatrix[index2][index1];
                (*complexNetwork).colorsMatrix[index2][index1] = (*complexNetwork).colorsMatrix[index2][index];
                (*complexNetwork).colorsMatrix[index2][index] = value;
            }

            // troca linhas
            for (index2 = 0; index2 < sizeSimilarityMatrix; index2++){              
                value = (*complexNetwork).colorsMatrix[index1][index2];
                (*complexNetwork).colorsMatrix[index1][index2] = (*complexNetwork).colorsMatrix[index][index2];
                (*complexNetwork).colorsMatrix[index][index2] = value;
            }          
        }     
    }  
}

int findInArray(int* array, int item, int sizeOfArray){    
     int i;
     for (i = 0; i < sizeOfArray; i++){ 
         if(array[i] == item){
             return item;
         }
     }
     return -1;
}


void writeLog(int* parcialArestasRemovidas, int* totalArestasRemovidas){
    *parcialArestasRemovidas = *parcialArestasRemovidas + qtdNodeWithEdgeRemoved;
     *totalArestasRemovidas = qtdNodeWithEdgeRemoved + *totalArestasRemovidas;
     if(*parcialArestasRemovidas >= N_ARESTAS){
        showtimeWithParameter(" Arestas removidas ", *totalArestasRemovidas);
        *parcialArestasRemovidas = 0;
     }
}

int isThereEdge(int** adjacencyMatrix){
    int true = 1;
    int false = 0;
    int i,j;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = 0; j < sizeSimilarityMatrix; j++){
            if(adjacencyMatrix[i][j] == 1){
                return true;
            }
        }
     }

    return false;
}


void removeEdgesWithBiggerBetweenness(int nodeWithEdgeRemoved[1][2], float** betweennessMatrix, int** adjacencyMatrix, int** neighborhoodMatrix){
    int i,j;
    float biggerValue = 0;
    int already[sizeSimilarityMatrix*100][2];
    int sizeOfAlready = 0;
    int indexToRemove;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i+1; j < sizeSimilarityMatrix; j++){
            if(betweennessMatrix[i][j] > biggerValue){
               biggerValue = betweennessMatrix[i][j];
            }
        }
    }

    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i+1; j < sizeSimilarityMatrix; j++){
            if(betweennessMatrix[i][j] == biggerValue){  
                 already[sizeOfAlready][0] = i;
                 already[sizeOfAlready][1] = j;
                 sizeOfAlready++;
            }
        }
    }  

    srand(time(NULL));
    srand(100);
    indexToRemove = rand() % sizeOfAlready;
    i = already[indexToRemove][0];
    j = already[indexToRemove][1];
    adjacencyMatrix[i][j] = 0;
    adjacencyMatrix[j][i] = 0;
    nodeWithEdgeRemoved[qtdNodeWithEdgeRemoved][0] = i;
    nodeWithEdgeRemoved[qtdNodeWithEdgeRemoved][1] = j;
    qtdNodeWithEdgeRemoved++;
}

void adjustNeighborhoodMatrix(int nodeWithEdgeRemoved[1][2], int** adjacencyMatrix, int** neighborhoodMatrix, int diameter){
    int indiceI, indiceJ, p, i, d, v1,v2, pathSize, item;
    int vertex[sizeSimilarityMatrix];
    int sizeVertex = 0, removeEdge = 1, sizeAux = 0, sizeEdges = 0, hasPath = 0;
    int aux, finded;
    
    indiceI = nodeWithEdgeRemoved[0][0];
    indiceJ = nodeWithEdgeRemoved[0][1];
    neighborhoodMatrix[indiceI][indiceJ] = 0;
    neighborhoodMatrix[indiceJ][indiceI] = 0;
    vertex[sizeVertex] = indiceI;
    sizeVertex++;
    vertex[sizeVertex] = indiceJ;
    sizeVertex++;

    removedEdges = (edge *)malloc( (sizeSimilarityMatrix * 300) * sizeof(edge));
    
    if(removedEdges == NULL){
         printf("\n ERRO NA CRIAÇÃO DA LISTA DE ARESTAS. PROBLEMA NA ALOCAÇÃO DE MEMÓRIA. \n");
         exit(1);
    }
    removedEdges[sizeEdges].V1 = indiceI;
    removedEdges[sizeEdges].V2 = indiceJ;
    sizeEdges++;
    
    // IDENTIFICA OS VERTICES QUE TIVERAM SUAS DISTANCIAS AFETADAS COM A REMOÇÃO.     
    for (d = 1; d < diameter; d++){
        sizeAux = 0;
        for (i = 0; i < sizeVertex; i++){ 
            item = vertex[i];
            for (v1 = 0; v1 < sizeSimilarityMatrix; v1++){
                if(neighborhoodMatrix[item][v1] == (d+1) ){
                    for (v2 = 0; v2 < sizeSimilarityMatrix; v2++){
                         removeEdge = 1;
                         if( ( (neighborhoodMatrix[item][v2] == 1)&&(neighborhoodMatrix[v1][v2] == d) ) || ( (neighborhoodMatrix[v1][v2] == 1)&&(neighborhoodMatrix[item][v2] == d) ) )  {
                             removeEdge = 0;
                             break;
                         }
                     }
                     if(removeEdge){
                         neighborhoodMatrix[item][v1] = 0;
                         neighborhoodMatrix[v1][item] = 0;
                         finded = 0;
                        
                         for(aux = 0; aux < (sizeVertex-1)+sizeAux; aux++){
                             if(vertex[aux] == v1){
                                 finded = 1;
                                 break;
                             }
                         }
                      //   printf("\n teste -1.3.2, finded = %d, soma=%d, v1 = %d, vertex = %d", finded, (sizeVertex-1)+(sizeAux+1), v1, vertex[(sizeVertex-1)+(sizeAux+1)]);
                         if(finded == 0){
                             sizeAux++;
                             vertex[(sizeVertex-1)+sizeAux] = v1;
                         }
                    //     printf("\n teste -1.3.3 | item = %d, v1 = %d, sizeEdges = %d, removedEdges[sizeEdges] = %d", item, v1, sizeEdges,removedEdges[sizeEdges] );
                      
                         if(item < v1){
                              removedEdges[sizeEdges].V1 = item;
                              removedEdges[sizeEdges].V2 = v1;
                         }else{
                              removedEdges[sizeEdges].V1 = v1;
                              removedEdges[sizeEdges].V2 = item;
                         }
                         sizeEdges++;   
                      //   printf("\n Size: %d", sizeEdges);
                     } 
                }
            }
        }
        sizeVertex = sizeVertex + sizeAux;
    }   

    // CALCULA A NOVA DISTANCIA ENTRE OS VERTICES AFETADOS.
    for (d = 2; d <= (2*diameter); d++){
        pathSize = 0;
        for (i = 0; i < sizeEdges; i++){ 
            v1 =  removedEdges[i].V1;
            v2 =  removedEdges[i].V2;            
            if(neighborhoodMatrix[v1][v2] == 0 ){
                pathSize = 1;
                hasPath = 0;
                while( (hasPath==0) && (pathSize < d) ){
                    for (p = 0; p < sizeSimilarityMatrix; p++){                        
                         if( ( (neighborhoodMatrix[v1][p] == pathSize)&&(neighborhoodMatrix[v2][p] == (d-pathSize)) ) || ( (neighborhoodMatrix[v1][p] == (d-pathSize))&&(neighborhoodMatrix[v2][p] == pathSize) ) )  {
                             hasPath = 1;
                             break;
                         }                        
                    }
                 
                    if(hasPath){
                        neighborhoodMatrix[v1][v2] = d;
                        neighborhoodMatrix[v2][v1] = d;
                    }
                    
                    pathSize++;
                }             
            }
        }
        // Todos os caminhos já foram recalculados.
        if(pathSize == 0){
            break;
        }
    }           
    free(removedEdges);
}

void calculateBetweenness(int diameter, complexNetwork complexNetwork, float** betweennessMatrix){
    int i, j, l, k, m;
    copyMatrixOfFloat(complexNetwork.adjacencyMatrix, betweennessMatrix);
    int** D = complexNetwork.neighborhoodMatrix;
    
    int requirementOK[sizeSimilarityMatrix];
    memset(requirementOK, 0x7f, sizeof(requirementOK));
    int reqCount = 0;
    
    for (l = diameter; l > 1; l--){
        for (i = 0; i < sizeSimilarityMatrix; i++){
              for (j = (i+1); j < sizeSimilarityMatrix; j++){
                  
                  if( D[i][j] == l ){
                      reqCount = 0;
                      for (k = 0; k < sizeSimilarityMatrix; k++){
                           if( ( (D[i][k] == 1) && (D[j][k] == (l-1))  ) || ( (D[i][k] == (l-1)) && (D[j][k] == 1) ) ){
                               requirementOK[reqCount] = k;
                               reqCount++;
                           }
                       }
                       
                       for (m = 0; m < reqCount; m++){
                           k = requirementOK[m];
                           betweennessMatrix[i][k] = betweennessMatrix[i][k] + ((betweennessMatrix[i][j] + 1)/reqCount);
                           betweennessMatrix[j][k] = betweennessMatrix[j][k] + ((betweennessMatrix[i][j] + 1)/reqCount);
                           betweennessMatrix[k][i] = betweennessMatrix[i][k];
                           betweennessMatrix[k][j] = betweennessMatrix[j][k]; 
                       }
                  }
              }
         }
    }

    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i; j < sizeSimilarityMatrix; j++){
            if (complexNetwork.adjacencyMatrix[i][j] == 0){
                 betweennessMatrix[i][j] = 0;
                 betweennessMatrix[j][i] = 0;
            }         
        }
    }
}


void copyMatrixOfFloat(int** sourceMatrix, float** targetMatrix){
    int i,j;
    
   // #pragma omp parallel for private(i,j) num_threads(numThreads)
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i; j < sizeSimilarityMatrix; j++){
            targetMatrix[i][j] = (float)sourceMatrix[i][j];
            targetMatrix[j][i] = (float)sourceMatrix[j][i];
        }
    }   
}

void copyMatrix(int** sourceMatrix, int** targetMatrix){
    int i,j;
    
  //  #pragma omp parallel for private(i,j) num_threads(numThreads)
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i+1; j < sizeSimilarityMatrix; j++){
            targetMatrix[i][j] = sourceMatrix[i][j];
            targetMatrix[j][i] = sourceMatrix[j][i];
        }
    }    
}

void clearMatrix(int** matrix){
    int i,j;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i; j < sizeSimilarityMatrix; j++){
            matrix[i][j] = 0;
            matrix[j][i] = 0;
        }
    }    
}

void clearMatrixBySize(int** matrix, int size){
    int i,j;
    
    for (i = 0; i < size; i++){ 
        for (j = i; j < size; j++){
            matrix[i][j] = 0;
            matrix[j][i] = 0;
        }
    }
}

void clearVetor(int* vetor){
    int i;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        vetor[i] = 0;
    }    
}
