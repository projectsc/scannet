/* 
    O SCANNET visa simplificar a utilização do método PSN (Protein Similarity Networks)
	a partir da unificação dos seus passos.
	
    Copyright (C) 2015 Austeclino Magalhães Barros Júnior

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

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
void loadSequenceList(char[], char**);
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
void removeEdgesWithBiggerBetweenness(int[][2], float**, int**);
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
int countNumberLinesOfFile(char[30]);
void plotGraphic();
double getTime();
complexNetwork createNetworkTeste();
void setupOutputNames(char*, int, char[10]);
void saveNetwork(int**, char[30]);
void orderSequences(complexNetwork*);
void saveArrayInFile(char**, char[], int);
char** createDynamicCharMatrix(int);

//Global vars
enum boolean {
    true = 1, false = 0
};

typedef enum boolean bool;

int sizeSimilarityMatrix = 1;
int** similarityMatrix;
int NCutPoints = 101;
int numThreads = 4;
int isInCluster = 1;
int N_EDGES = 100;
int indexMelhorRede = 0;
int qtdDendrogramLines = 0;
complexNetwork *complexNetworks;
complexNetwork SelectBetterNetworkToFindCommunities(FILE*);
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
int optionLoadSequenceFile = 0;
char** sequenceArray;

double start_t_RC, end_t_RC, tempo_total_RC, start_t_Distancia, end_t_Distancia, 
       tempo_total_Distancia, start_t_NGA, end_t_NGA, tempo_total_NGA, start_t_MC, 
       end_t_MC, tempo_total_MC, start_t_TOTAL, end_t_TOTAL, tempo_total_TOTAL ;

double start_t_Paralelo, end_t_Paralelo, total_t_Paralelo;
double * tempos_paralelo;
char filename[21];
char sequenceFilename[21];
FILE *arqLog;


int main(int argc, char** argv) {
        complexNetwork betterNetwork;
        int option = 1;       
		
		printf(" Copyright (C) 2015 Austeclino Magalhaes Barros Junior");
		printf("\n\n This program is free software: you can redistribute it and/or modify");
		printf("\n it under the terms of the GNU General Public License as published by");
		printf("\n the Free Software Foundation, either version 3 of the License, or");
		printf("\n (at your option) any later version.");
     	printf("\n\n This program is distributed in the hope that it will be useful,");
		printf("\n but WITHOUT ANY WARRANTY; without even the implied warranty of");
		printf("\n MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the");
		printf("\n GNU General Public License for more details.");
		printf("\n\n You should have received a copy of the GNU General Public License");
		printf("\n along with this program.  If not, see <http://www.gnu.org/licenses/>.");
		
		printf("\n\n\n --------------------------- SCANNET v1.0 ------------------------------ \n\n");
		
        while (option == 1 || option == 2){
        
            if(option == 1){
                
                printf("\n\n Enter the similarity file name (Ex: input.txt):  ");
                scanf("%s", filename);

                sizeSimilarityMatrix = countNumberLinesOfFile(filename);
				getcwd(current_path_entrada, 255);
               
                #ifdef WINDOWS
					strcat(current_path_entrada, "\\");
				#else
					strcat(current_path_entrada, "/");
				#endif
				
                similarityMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
                loadSimilarityMatrix(filename, similarityMatrix);   
                start_t_TOTAL = getTime();
                optionLoadSequenceFile = 0;
				
                while((optionLoadSequenceFile != 1)&&(optionLoadSequenceFile != 2)){
                    printf("\n\n Do you want to load the initial order of proteins (1-Yes/2-No):  ");
                    scanf("%d", &optionLoadSequenceFile);

                    if(optionLoadSequenceFile == 1){
                         printf("\n\n Enter the order of protein file name (Ex: initialOrder.csv):  ");
                         scanf("%s", sequenceFilename);
                         sequenceArray = createDynamicCharMatrix(sizeSimilarityMatrix);
                         loadSequenceList(sequenceFilename, sequenceArray);   
                    }else if((optionLoadSequenceFile != 1)&&(optionLoadSequenceFile != 2)){
                         printf("\n\n Invalid Option  ");
                    }                 
                 }
                
                start_t_RC = getTime();  
                printf("\n\n The complex networks are being created, please wait... \n");

                createComplexNetworks();
                end_t_RC = getTime();
                tempo_total_RC = end_t_RC - start_t_RC;
                
                char path[1000];
                strcpy(path,current_path_entrada);
                char log[30] = "log_";
                setupOutputNames(log, 0, ".txt");
                strcat(path,log);
                arqLog = fopen(path, "wt");
                fprintf(arqLog,"\n\n ------- Summary execution --------- \n");
            }

		    start_t_Distancia = getTime();
            betterNetwork = SelectBetterNetworkToFindCommunities(arqLog);
            end_t_Distancia = getTime();
            tempo_total_Distancia = end_t_Distancia - start_t_Distancia;
            printf("\n Number of edges: %d ", betterNetwork.qtdOfEdges);

			int** neighborhoodMatrixCopy = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix); 
			copyMatrix(betterNetwork.neighborhoodMatrix, neighborhoodMatrixCopy);
			
            betterNetwork.colorsMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);      
            betterNetwork.dendrogramMatrix = createDynamicMatrix(betterNetwork.qtdOfEdges+2, sizeSimilarityMatrix);
            betterNetwork.dendrogramGraphic = createDynamicMatrixFloat(betterNetwork.qtdOfEdges*2, sizeSimilarityMatrix+1);
            betterNetwork.realIndex = createDynamicArray(sizeSimilarityMatrix);
            betterNetwork.realIndexDendrogramGraphic = createDynamicArray(sizeSimilarityMatrix);

            int i = 0;
            for (i = 0; i < sizeSimilarityMatrix; i++){ 
                betterNetwork.realIndex[i] = i;
                betterNetwork.realIndexDendrogramGraphic[i] = i;
            }
			
			char networkName[30] = "network_";
            setupOutputNames(networkName, 1, ".csv");
            saveNetwork(betterNetwork.neighborhoodMatrix, networkName);
			
            copyMatrix(betterNetwork.neighborhoodMatrix, betterNetwork.colorsMatrix);
            printf("\n\n The Girvan-Newman algorithm is running, please wait ...\n\n");
            executeNga(betterNetwork);
            end_t_TOTAL = getTime();
            tempo_total_TOTAL = end_t_TOTAL - start_t_TOTAL;

            for (i = 0; i < sizeSimilarityMatrix; i++){ 
               betterNetwork.dendrogramMatrix[betterNetwork.qtdOfEdges+1][i] = betterNetwork.realIndex[i]+1;
               betterNetwork.dendrogramGraphic[qtdDendrogramLines+1][i] = betterNetwork.realIndexDendrogramGraphic[i]+1;     
            }

            setupEixoYDendrogram(&betterNetwork, qtdDendrogramLines);

            char mColorName[30] = "mColors_";
            setupOutputNames(mColorName, 1, ".dat");
            saveMatrizInFile(betterNetwork.colorsMatrix, mColorName, sizeSimilarityMatrix, sizeSimilarityMatrix );

            char mDendrogramGraphic[30] = "mDendrogramGraphic_";
            setupOutputNames(mDendrogramGraphic, 1, ".dat");
            saveMatrizInFileFloat(betterNetwork.dendrogramGraphic, mDendrogramGraphic,qtdDendrogramLines+2, sizeSimilarityMatrix+1);

            if(optionLoadSequenceFile == 1){
                char sequenceOutput[30] = "orderedSequence_";
                setupOutputNames(sequenceOutput, 1, ".csv");
                saveArrayInFile(sequenceArray, sequenceOutput,sizeSimilarityMatrix);
            }
            
            showtime("End of execution");  
            fprintf(arqLog,"\n\n The selected optimal network has %d edges", betterNetwork.qtdOfEdges);
            fprintf(arqLog,"\n Time spent on the creation of 101 complex networks (ms): %.5lf", tempo_total_RC); 
            fprintf(arqLog,"\n Time spent in the ideal complex network selection (ms): %.5lf", tempo_total_Distancia); 
            fprintf(arqLog,"\n Time spent on NGA (ms): %.5lf", tempo_total_NGA);
            fprintf(arqLog,"\n Time spent in generating the color matrix (ms): %.5lf", tempo_total_MC);
            fprintf(arqLog,"\n Running time (ms): %.5lf \n", tempo_total_TOTAL);
            fprintf(arqLog,"\n --------------------------------------------------- \n\n\n");

            printf("\n\n What do you want to do now?  ");
            printf("\n\n 1 - Analyze new sequence  ");
            printf("\n 2 - Analyze this sequence with other similarity value ");
            printf("\n 3 - Exit \n");
            scanf("%d", &option);
            if(option == 1 || option == 3){
                 fclose(arqLog);
            }else{
				copyMatrix(neighborhoodMatrixCopy, betterNetwork.neighborhoodMatrix);
			}
            
        } 
        free(sequenceArray);
        free(similarityMatrix);
        free(complexNetworks);
        return 0;
}


void saveNetwork(int** matrix , char networkName[30]){
      char path[1000];
      int i, j;
	  FILE *arqNetwork;
      strcpy(path,current_path_entrada);
      strcat(path,networkName);
      arqNetwork = fopen(path, "wt");
      fprintf(arqNetwork,"source;target;type\n");
      
       for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for(j = (i+1); j < sizeSimilarityMatrix; j++){
           if(matrix[i][j] == 1){
               fprintf(arqNetwork,"%d;%d;undirected\n",i+1,j+1);
           } 
        } 
    }

    fclose(arqNetwork);
    printf("\n\n The file %s was successfully saved! \n", networkName);
}

void setupOutputNames(char * nameOfFile, int concatNetwork, char format[10] ){
    
    int sizeFileName = strlen (filename);
    char fileNameWithoutExtension[30] = "";
    char buffer [33];
    
    int i = 0;
    for (i = 0; i < sizeFileName; i++){
        if(filename[i] == '.'){
            break;
        }else{
            sprintf(buffer, "%c", filename[i]);
            strcat(fileNameWithoutExtension, buffer);
        }
    }
    strcat(nameOfFile, fileNameWithoutExtension); 
    
    if(concatNetwork == 1){
        strcat(nameOfFile, "_");
        char index[30];
        sprintf(index, "%d", indexMelhorRede);
        strcat(nameOfFile, index);
    }
    strcat(nameOfFile, format);
}

int countNumberLinesOfFile(char nameOfFile[30]){   
    FILE *fileptr;
    int numberLines = 1;
    int chr;
	int chrBefore;
    fileptr = fopen(nameOfFile, "r");
    chr = getc(fileptr);
    while ((chr != EOF) && (chr != -1)){
        if (chr == '\n' && chrBefore != '\n'){
            numberLines = numberLines + 1;
        }
		chrBefore = chr;
        chr = getc(fileptr);
    }
	if (chrBefore == '\n'){
		numberLines = numberLines - 1;
    }
    fclose(fileptr);
    return numberLines;
}

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
    printf("\n%s: %d  Time: %02d:%02d:%02d \n",texto, value, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);   
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
         printf("\n  Error in memory allocation. \n");
         exit(1);
    }
        
    for(i=0; i < sizeLine; i++){
        matrix[i] = calloc(sizeColumn, sizeof(int));
        if(matrix[i] == NULL){
            printf("\n Error in memory allocation. \n");
            exit(1);
        }
    }
    
    return matrix;
}

int* createDynamicArray(int size){
    int* matrix = calloc(size, sizeof(int));     
  
    if(matrix == NULL){
         printf("\n Error in memory allocation. \n");
         exit(1);
    }
    
    return matrix;
}

char** createDynamicCharMatrix(int size){
    char** matrix = malloc(size*sizeof(char*));
    int i;
  
    if(matrix == NULL){
         printf("\n Error in memory allocation. \n");
         exit(1);
    }
        
    for(i=0; i < size; i++){
        matrix[i] = calloc(1000, sizeof(int));
        if(matrix[i] == NULL){
            printf("\n Error in memory allocation. \n");
            exit(1);
        }
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
         printf("\n Error in memory allocation. \n");
         exit(1);
    }
    
    for(i=0; i < sizeLine; i++){
        matrix[i] = calloc(sizeColumn, sizeof(float));
         if(matrix[i] == NULL){
            printf("\n Error in memory allocation. \n");
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
    int size = sizeSimilarityMatrix * 4;
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
        printf("\n\n\n The file with the similarity matrix can not be opened \n");
        system("pause");
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
         
   fclose(arqBLAST);
} 

void loadSequenceList(char nameOfFile[30], char** array){
    FILE *arqSequence; 
    char ch;
    int n = 0;
    int indexLine = 0;
  
    char path[1000];
    strcpy(path,current_path_entrada);  
    strcat(path,nameOfFile);
    
    arqSequence = fopen(path, "rt");
    
    if (arqSequence == NULL){
        printf("\n\n\n The file with the sequence list can not be opened \n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    while (indexLine < sizeSimilarityMatrix) {
        ch = getc(arqSequence);
       
	if(ch == '\n'){
            indexLine++;
            n = 0;  
        }else{
            array[indexLine][n] = ch;
            n++;
        }       
    }
         
   fclose(arqSequence);
} 


void saveMatrizInFile(int** matriz, char nameOfFile[30], int n, int m){         
    FILE *arqMATRIZ;
    int w,k,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
   
    if (arqMATRIZ == NULL){
        printf("\n Problems in the creation of the similarity matrix file.\n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    for (w = 0; w < n; w++){
        for (k = 0; k < m; k++){
            if(matriz[w][k] < 10){
                result = fprintf(arqMATRIZ,"%d   ", matriz[w][k]);
            }else{
                result = fprintf(arqMATRIZ,"%d  ", matriz[w][k]);
            }
            if (result == EOF) printf("\n Error in file recording\n");
        }

        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("\n Error in the recording of the similarity matrix file \n");
    }

    fclose(arqMATRIZ);

    printf("\n\n\n The file %s was successfully saved! \n", nameOfFile);
}
    
void saveMatrizInFileFloat(float** matriz, char nameOfFile[30], int n, int m){         
    FILE *arqMATRIZ;
    int w,k,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
    
    if (arqMATRIZ == NULL){ 
        printf("\n Problems in the creation of the similarity matrix file.\n");
        system("pause");
        exit (EXIT_FAILURE);
    } 

    for (w = 0; w < n; w++){
        for (k = 0; k < m; k++){
            result = fprintf(arqMATRIZ,"%.2f ", matriz[w][k]);
            if (result == EOF) printf("Error in file recording\n");
        }

        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("\n Error in the recording of the similarity matrix file\n");
    }

    fclose(arqMATRIZ);

    printf("\n\n The file %s was successfully saved! \n", nameOfFile);
}

void saveArrayInFileFloat(float* array, char nameOfFile[30], int size){         
    FILE *arqMATRIZ;
    int w,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
   
    if (arqMATRIZ == NULL){
        printf("\n Problems in the creation of the distances file \n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    for (w = 0; w < size; w++){
        result = fprintf(arqMATRIZ,"%.2f", array[w]);
        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("\n Error in file recording\n");
    }

    fclose(arqMATRIZ);

    printf("\n\n\n The file %s was successfully saved! \n", nameOfFile);
}

void saveArrayInFile(char** array, char nameOfFile[30], int size){         
    FILE *arqSequence;
    int i, result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqSequence = fopen(path, "wt");
   
    if (arqSequence == NULL){
        printf("\n Problems in the creation of the sequence file \n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    for (i = 0; i < size; i++){
        for (i = 0; i < size; i++){
            result = fprintf(arqSequence,"%s", array[i]);
            result = fprintf(arqSequence,"\n");
            if (result == EOF) printf("\n Error in file recording \n");
        }
        
    }

    fclose(arqSequence);
    printf("\n\n\n The file %s was successfully saved! \n", nameOfFile);
}

void createComplexNetworks(){
    complexNetworks = (complexNetwork *)malloc(101 * sizeof(complexNetwork));
    if(complexNetworks == NULL){
         printf("\n Error in memory allocation. \n");
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

void dijkstra (int Vi, int** adjacencyMatrix, int** neighborhoodMatrix, int maxSize)
{       
    int dis[maxSize];
    char vis[maxSize];
    memset (vis, 0, sizeof (vis));
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

complexNetwork SelectBetterNetworkToFindCommunities(FILE* arqLog){
      float distance = 0;
      float biggerDistance = 0;
      int i = 0;
      complexNetwork betterNetwork;
      float distances[101];
      int i_distance = 0;
      
      for (i = 0; i < (NCutPoints-1); i++){ 
          distance =  calculateDistance(complexNetworks[i], complexNetworks[i+1]);  
          fprintf(arqLog,"\n Distance between the network %d and %d = %f", i, i+1,distance);
          printf("\n Distance between the network %d and %d = %f", i, i+1,distance);
          distances[i_distance] = distance;
          i_distance = i_distance + 1;
          if(distance > biggerDistance){
              biggerDistance = distance;
              indexMelhorRede = i;
          }
      }
      
      char mDistanceGraphic[30] = "mDistanceGraphic_";
      setupOutputNames(mDistanceGraphic, 0, ".dat");
      saveArrayInFileFloat(distances, mDistanceGraphic ,i_distance);
      
      plotGraphic();
      
      fprintf(arqLog, "\n\n The similarity value %d has the largest distance: %f  ", indexMelhorRede,  biggerDistance);
      printf("\n\n The similarity value %d has the largest distance: %f  ", indexMelhorRede,  biggerDistance);
      
      printf("\n\n\n Select the similarity value you want to use:  ");
      scanf("%d", &indexMelhorRede);
      
      fprintf(arqLog, "\n\n Similarity value selected: %d  ", indexMelhorRede);
      
      betterNetwork = complexNetworks[indexMelhorRede];
      return betterNetwork;
}


void plotGraphic(){

    char path[1024];
    char syscommand[1024];
    strcpy(path,current_path_entrada);
    char mDistanceGraphic[100] = "mDistanceGraphic_";
    setupOutputNames(mDistanceGraphic, 0, ".dat");
    strcat(path,mDistanceGraphic);
	
    #ifdef WIN32
        sprintf(syscommand, "start gnuplot -persist -e \"plot '%s' with lines title 'SCANNET - Distances chart'; set xlabel 'Similarity Value'; set ylabel 'Distance'; set xtics 2; set ytics 0.2; pause -1 'Press enter to exit. . .';\"\"",path);
	#else
		sprintf(syscommand, "gnuplot -p -e \"plot '%s' with lines title 'SCANNET - Distances chart';set xlabel 'Similarity Value'; set ylabel 'Distance'; set xtics 2; set ytics 0.2;\"",path);
	#endif
	
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
        removeEdgesWithBiggerBetweenness(nodeWithEdgeRemoved, betweennessMatrix, complexNetwork.adjacencyMatrix);  
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
    if(optionLoadSequenceFile == 1){
        orderSequences(&complexNetwork);
    }
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
              
                for (j = 0; j < sizeSimilarityMatrix+1; j++){
                    (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][j] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines - 1][j];
                }
                (*complexNetwork).numberDendrogramLines = (*complexNetwork).numberDendrogramLines + 1;
            }
            
            for (j = 0; j < index; j++){
               realIndex = findRealIndex((*complexNetwork).realIndexDendrogramGraphic, newCommunityItens[j]);
               (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][realIndex] = pontoMaximo + ((index+1) / dividendo);  
            }
         
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
            
            index1 = findRealIndex(arrayToSort, index1);
            value = arrayToSort[index1];
            arrayToSort[index1] = arrayToSort[index];
            arrayToSort[index] = value;
            
            for (index2 = 0; index2 < sizeSimilarityMatrix; index2++){              
                value = (*complexNetwork).colorsMatrix[index2][index1];
                (*complexNetwork).colorsMatrix[index2][index1] = (*complexNetwork).colorsMatrix[index2][index];
                (*complexNetwork).colorsMatrix[index2][index] = value;
            }

            for (index2 = 0; index2 < sizeSimilarityMatrix; index2++){              
                value = (*complexNetwork).colorsMatrix[index1][index2];
                (*complexNetwork).colorsMatrix[index1][index2] = (*complexNetwork).colorsMatrix[index][index2];
                (*complexNetwork).colorsMatrix[index][index2] = value;
            }          
        }     
    }  
}

void orderSequences(complexNetwork *complexNetwork){
    int index = 0, index1, i, valueInt;
    char* value = 0;
    
    int * arrayToSort = createDynamicArray(sizeSimilarityMatrix);
  
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        arrayToSort[i] = i;
    } 
     
    for (index = 0; index < sizeSimilarityMatrix; index++){
        index1 = (*complexNetwork).realIndexDendrogramGraphic[index];
        if(index1 != index){     
            index1 = findRealIndex(arrayToSort, index1);
            
            valueInt = arrayToSort[index1];
            arrayToSort[index1] = arrayToSort[index];
            arrayToSort[index] = valueInt;
            
            value = sequenceArray[index];
            sequenceArray[index] = sequenceArray[index1];
            sequenceArray[index1] = value;  
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
     if(*parcialArestasRemovidas >= N_EDGES){
        showtimeWithParameter(" Edges removed", *totalArestasRemovidas);
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


void removeEdgesWithBiggerBetweenness(int nodeWithEdgeRemoved[1][2], float** betweennessMatrix, int** adjacencyMatrix){
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

    removedEdges = (edge *)malloc( (sizeSimilarityMatrix * 500) * sizeof(edge));
    
    if(removedEdges == NULL){
         printf("\n ERROR IN CREATING THE EDGES LIST. PROBLEM IN MEMORY ALLOCATION. \n");
         exit(1);
    }
    removedEdges[sizeEdges].V1 = indiceI;
    removedEdges[sizeEdges].V2 = indiceJ;
    sizeEdges++;
       
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

                         if(finded == 0){
                             sizeAux++;
                             vertex[(sizeVertex-1)+sizeAux] = v1;
                         }
                      
                         if(item < v1){
                              removedEdges[sizeEdges].V1 = item;
                              removedEdges[sizeEdges].V2 = v1;
                         }else{
                              removedEdges[sizeEdges].V1 = v1;
                              removedEdges[sizeEdges].V2 = item;
                         }
                         sizeEdges++;   
                     } 
                }
            }
        }
        sizeVertex = sizeVertex + sizeAux;
    }   

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
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i; j < sizeSimilarityMatrix; j++){
            targetMatrix[i][j] = (float)sourceMatrix[i][j];
            targetMatrix[j][i] = (float)sourceMatrix[j][i];
        }
    }   
}

void copyMatrix(int** sourceMatrix, int** targetMatrix){
    int i,j;
    
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
