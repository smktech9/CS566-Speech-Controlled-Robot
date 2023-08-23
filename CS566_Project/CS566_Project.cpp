// CS566_Project.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <conio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <windows.h>
#include<string>

using namespace std;

#define S 320 // Samples in one frame
#define W 80
#define P 12
#define maxF 160 // Count of stable frames
#define pi 3.1415926535
#define D 0.00001
#define epsilon 0.03
#define M 32
#define N 5
#define maxT 160

// array for storing normalized data, energy, steady frames and average Ci
long double normal[100000], energy[100000], steadyFrames[maxF][S];
long int nSize, eSize, rIndex = 0;
long int F;
long int universeSize;

// for codeBook
long double input[100000][P];	// Array to store input vectors
long double codeBook[M][P];	// CodeBook

// array for stroring Ri, Ai & Ci
double Ri[maxF][P+1], Ai[maxF][P+1], Ci[maxF][P+1];

// Tokhura Weights
double tokhuraWeights[P] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

// for HMM
long int T;
long double A[N+1][N+1], B[N+1][M+1], PI[N+1];
long int O[maxT+1], Q[maxT+1];
long double alpha[maxT+1][N+1], beta[maxT+1][N+1], gamma[maxT+1][N+1], delta[maxT+1][N+1], psi[maxT+1][N+1], xi[maxT+1][N+1][N+1];
long double new_PI[N+1], new_A[N+1][N+1], new_B[N+1][M+1];

void Normalize(char *filename)
{
	FILE *myfile;
	myfile = fopen(filename, "r");
	long double read, sum = 0, max = INT_MIN;
	int count = 0;
	long double dcshift;
	FILE *fp;
	while(fscanf(myfile, "%Lf", &read) != EOF)  // Finding max value and DC Shift
	{
			if(abs(read) > max)
			{
				max = abs(read);
			}
			normal[count] = read;
			sum += read;
			count++;
	}
	dcshift = sum / count;  // DC Shift
	nSize = count;
	int index = 0, eMaxI = -1;
	long double fsum = 0, eMax = -1;
	int i = 0, j = 0;
	for(int x = 0; x < count; x++)
	{
		normal[x] = (normal[x] - dcshift) * 5000 / max;
		//printf("%Lf\n", normal[x]);
	}
	while(i < count) // Normalizing input & calculating energy for each frame
	{
		fsum += normal[i] * normal[i];
		i++;
		j++;
		if(j % S == 0)
		{
			energy[index] = fsum / S; // Energy
			//printf("%Lf\n", energy[index]);
			fsum = 0;
			if(energy[index] > eMax) // Max energy
			{
				eMax = energy[index];
				eMaxI = i;
			}
			index++;
			i = i - S + W;
			j = 0;
		}
	}
	if(fsum != 0)
	{
		energy[index] = fsum / j;
		fsum = 0;
		if(energy[index] > eMax)
		{
			eMax = energy[index];
			eMaxI = i;
		}
		index++;
	}
	fclose(myfile);
	eSize = index;
	int k = 0, f = 0;
	int start = eMaxI - S - (maxF / 2) * W;
	int end = eMaxI + (maxF / 2) * W;
	if(start < 0)
	{
		start = 0;
		end = S + W * (maxF - 1);
	}
	else if(end >= count)
	{
		end = count - 1;
		start = count - 1 - (S + W * (maxF - 1));
	}
	if(start < 0 || end >= count)
	{
		start = 0;
		end = count - 1;
	}
	for(int i = start; i <= end; i++) // Storing Steady Frmaes
	{
		steadyFrames[f][k++] = normal[i];
		//printf("%Lf ", steadyFrames[f][k-1]);
		if(k == S)
		{
			k = 0;
			f++;
			i = i - S + W;
			//printf("\n");
		}
		if(f == maxF)
		{
			break;
		}
	}
	T = f;
	F = f;
	fclose(myfile);
}

void RiAiCi()
{
	// Calculating Ri

	for(int i = 0; i < F; i++)
	{
		for(int m = 0; m <= P; m++)
		{
			Ri[i][m] = 0;
			for(int n = 0; n < N - m; n++)
			{
				Ri[i][m] += steadyFrames[i][n] * steadyFrames[i][n + m];
			}
		}
	}

	// Calculating Ai using Durbin's Algorithm

	double E[P+1], sum = 0, A[P+1][P+1], K[P+1];
	for(int k = 0; k < F; k++)
	{
		E[0] = Ri[k][0];
		for(int i = 1; i <= P; i++)
		{
			for(int j = 1; j <= i - 1; j++)
			{
				sum += A[i-1][j] * Ri[k][i - j];
			}
			K[i] = (Ri[k][i] - sum) / E[i-1];
			sum = 0;
			A[i][i] = K[i];
			for(int j = 1; j <= i-1; j++)
			{
				A[i][j] = A[i-1][j] - K[i] * A[i-1][i-j];
			}
			E[i] = (1 - K[i] * K[i]) * E[i-1];
		}
		for(int i = 1; i <= P; i++)
		{
			Ai[k][i] = A[P][i];
		}
	}

	// Calculating Ci

	for(int i = 0; i < F; i++)
	{
		Ci[i][0] = log(Ri[i][0] * Ri[i][0]);
		for(int m = 1; m <= P; m++)
		{
			sum = 0;
			for(int k = 1; k <= m-1; k++)
			{
				sum += (k / m) * Ci[i][k] * Ai[i][m-k];
			}
			Ci[i][m] = Ai[i][m] + sum;
		}
	}

	// Applying raised-sin window

	for(int i=0;i<F;i++)
	{
		for(int j=1;j<=P;j++)
		{
			Ci[i][j] *= (1 + (P / 2) * sin((j * pi) / P));
			//printf("%Lf ", Ci[i][j]);
		}
	}

}

// Function to create Universe
void createUniverse()
{
	char fname[15];
	char move[5] = {'L', 'R', 'U', 'D'};
	char folder[2] = {'S', 'A'};

	FILE* myfile;
	myfile = fopen("Universe/universe.txt", "w");
	for(int i = 0; i < 4; i++)
	{
		for(int k = 0; k < 2; k++)
		{
			for(int j = 1; j <= 20; j++)
			{	
				sprintf(fname, "%c/%c%d.txt", folder[k], move[i], j);
				Normalize(fname); // Normalize data
				RiAiCi(); // Calculate Ri, Ai and Ci 
				printf("%s\n", fname);
				for(int k = 0; k < F; k++)
				{
					for(int l = 1; l <= P; l++)
					{
						fprintf(myfile, "%Lf ", Ci[k][l]);
						//printf("%Lf ", Ci[k][l]);
					}
					//printf("\n");
					fprintf(myfile, "\n");
				}
			}
		}
	}
	fclose(myfile);
}

// Function to read the input file
void readFileLBG()
{
	FILE *myfile;
	myfile = fopen("Universe/universe.txt", "r");
	long double read;
	int i = 0, j = 0;
	while(fscanf(myfile, "%Lf", &read)!=EOF)
	{
		input[i][j] = read;
		//printf("%Lf ", input[i][j]);
		j++;
		if(j == P)
		{
			//printf("\n");
			j = 0;
			i++;
		}
	}
	universeSize = i;
	//printf("%d", universeSize);
	fclose(myfile);
}

// Function to implement k-means algorithm
void kMeans(int size)
{
	long double currDist = 1, prevDist = 0, totalDist = 0; 
	int iter = 0;
	while(fabs(currDist - prevDist) > D)
	{
		int clusterSize[M] = {0};
		long double clusterSum[M][P] = {0}, max_cluster = INT_MIN;
		int minI = -1;
		totalDist = 0;
		int maxI = -1;
		for(int i = 0; i < universeSize; i++)
		{
			long double min = INT_MAX;
			for(int j = 0; j < size; j++)
			{
				long double tokDist[M] = {0};
				for(int k = 0; k < P; k++)
				{
					tokDist[j] += tokhuraWeights[k] * (input[i][k] - codeBook[j][k]) * (input[i][k] - codeBook[j][k]);	// Calculating Tokhura Distance
				}
				if(tokDist[j] < min)
				{
					min = tokDist[j];	// Finding minimum distortion
					minI = j;			// Minimum distortion index
				}
			}
			clusterSize[minI]++;	// Updating cluster size
			for(int l = 0; l < P; l++)
			{
				clusterSum[minI][l] += input[i][l];	// Updating sum of vectors in cluster
			}
			totalDist += min;
		}
		for(int x = 0; x < size; x++)
		{
			if(clusterSize[x] > max_cluster)
			{
				max_cluster = clusterSize[x];
				maxI = x;
			}
		}
		for(int m = 0; m < size; m++)
		{
			for(int n = 0; n < P; n++)
			{
				if(clusterSize[m] == 0)
				{
					time_t t1;
					srand ( (unsigned) time (&t1));
					codeBook[m][n] = (clusterSum[maxI][n] / max_cluster) * ((rand() % 10) + 5) / 10;
				}
				else
				{
					codeBook[m][n] = clusterSum[m][n]/clusterSize[m];	// Updating CodeBook
				}
			}
		}
		prevDist = currDist;
		currDist = totalDist / universeSize; 
		iter++;
	}
}

// Function to implement LBG algorithm
void LBG()
{
	readFileLBG();

	// Initialize codebook with average of the universe vectors

	//printf("\n");
	long double sum[P] = {0};
	for(int i = 0; i < P; i++)
	{	
		for(int j = 0; j < universeSize; j++)
		{
			sum[i] += input[j][i];
		}
		codeBook[0][i] = sum[i] / universeSize;
		//printf("%Lf ", codeBook[0][i]);
	}
	//printf("\n\n");
	int m = 1;	// Initial codebook size (m = 1)
	while(m != M)
	{
		// Double codeBook size and enter new values
		for(int i = 0; i < m; i++)	
		{
			for(int j = 0; j < P; j++)
			{
				codeBook[i + m][j] = codeBook[i][j] * (1 + epsilon);
				codeBook[i][j] = codeBook[i][j] * (1 - epsilon);
			}
		}
		m = m * 2;
		// Run k-means algorithm on current codeBook
		kMeans(m);
		for(int l = 0; l < m; l++)
		{
			for(int n = 0; n < P; n++)
			{
				printf("%Lf ", codeBook[l][n]);
			}
			printf("\n");
		}
		printf("\n");
	}
	FILE * myfile;
	myfile = fopen("Codebook/codebook.txt", "w");
	for(int i = 0; i < M; i++)
	{
		for(int j = 0; j < P; j++)
		{
			fprintf(myfile, "%Lf ", codeBook[i][j]);
		}
		fprintf(myfile, "\n");
	}
	fclose(myfile);
}

void tokhuraDistance()
{
	for(int f = 0; f < T; f++)
	{
		long double input;
		FILE* myfile;
		myfile = fopen("Codebook/codebook.txt", "r");
		long double tempD = 0, finalD = INT_MAX;
		int p = 1, i = 1;
		while(fscanf(myfile, "%Lf ", &input)!=EOF)
		{
			tempD += tokhuraWeights[p-1]*(Ci[f][p] - input)*(Ci[f][p] - input);
			p++;
			if(p > 12)
			{
				p = 1;
				if(tempD <= finalD)
				{
					finalD = tempD;
					O[f+1] = i;
					//printf("%ld ", O[f+1]);
				}
				i++;
				tempD = 0;
			}
		}
		fclose(myfile);
	}
}

void tokhuraDistanceDigit()
{
	for(int f = 0; f < T; f++)
	{
		long double input;
		FILE* myfile;
		myfile = fopen("Codebook/codebook2.txt", "r");
		long double tempD = 0, finalD = INT_MAX;
		int p = 1, i = 1;
		while(fscanf(myfile, "%Lf ", &input)!=EOF)
		{
			tempD += tokhuraWeights[p-1]*(Ci[f][p] - input)*(Ci[f][p] - input);
			p++;
			if(p > 12)
			{
				p = 1;
				if(tempD <= finalD)
				{
					finalD = tempD;
					O[f+1] = i;
					//printf("%ld ", O[f+1]);
				}
				i++;
				tempD = 0;
			}
		}
		fclose(myfile);
	}
}

void readFileHMM()
{
	FILE* myfile;
	myfile = fopen("InitialModel/B.txt", "r");
	int i = 1, j = 1;
	long double inp;
	while(fscanf(myfile, "%Lf", &inp)!=EOF)
	{
		B[i][j] = inp;
		//printf("%24.17Lf ", B[i][j]);
		j++;
		if(j > M)
		{
			j = 1;
			i++;
			//printf("\n");
		}
	}
	fclose(myfile);
	myfile = fopen("InitialModel/A.txt", "r");
	i = 1, j = 1;
	//printf("\n");
	//printf("\n");
	while(fscanf(myfile, "%Lf", &inp)!=EOF)
	{
		A[i][j] = inp;
		//printf("%24.17Lf ", A[i][j]);
		j++;
		if(j > N)
		{
			j = 1;
			i++;
			//printf("\n");
		}
	}
	fclose(myfile);
	myfile = fopen("InitialModel/PI.txt", "r");
	i = 1;
	//printf("\n");
	//printf("\n");
	while(fscanf(myfile, "%Lf", &inp)!=EOF)
	{
		PI[i] = inp;
		//printf("%24.17Lf ", PI[i]);
		i++;
		//printf("\n");
	}
	fclose(myfile);
}

// Function to solve problem 1
long double problem1()
{
	long double p = 0;

	// Calculating value of alpha
	
	//printf("Values of alpha - \n");
	for(int i=1;i<=N;i++)	// Initialization
	{
		//printf("%24.17E ", B[i][O[1]]);
		alpha[1][i] = PI[i] * B[i][O[1]];
		//printf("%24.17E ", alpha[1][i]);
	}
	//printf("\n");
	for(int t=1;t<T;t++)	// Iteration
	{
		for(int j=1;j<=N;j++)
		{
			long double sum = 0;
			for(int i=1;i<=N;i++)
			{
				sum += alpha[t][i] * A[i][j];
				//printf("%24.17E ", A[i][j]);
			}
			alpha[t+1][j] = sum * B[j][O[t+1]];
			//printf("\n");
			//printf("%24.17E ", alpha[t+1][j]);
		}
		//printf("\n");
	}
	
	// Calculating Probability of O given lambda
	
	for(int i=1;i<=N;i++)
	{
		p += alpha[T][i];
	}
	//printf("\n\n");
	//printf("\t\t\t\t\t\t****************************************************************************************************\n\n");
	//printf("\t\t\t\t\t\t\t\t\tProbability of (O/lambda) : %E\n\n", p);
	//printf("\t\t\t\t\t\t****************************************************************************************************\n");
	
	// Calculating beta values
	
	for(int i=1;i<=N;i++)	//Initialization
	{
		beta[T][i] = 1;
	}
	for(int t=T-1; t>0; t--)	// Iteration
	{
		for(int i=1;i<=N;i++)
		{
			long double sum = 0;
			for(int j = 1; j<=N; j++)
			{
				sum += A[i][j] * B[j][O[t+1]] * beta[t+1][j];
			}
			beta[t][i] = sum;
		}
	}
	//printf("\n\nValues of beta - \n");
	//for(int t=1;t<=T;t++)	// Printing beta
	//{
	//	for(int i=1;i<=N;i++)
	//	{
	//		printf("%24.17E ", beta[t][i]);
	//	}
	//	printf("\n");
	//}
	return p;
}

// Function to solve problem 2
void problem2()
{
	// Calculating delta and psi
	
	for(int i=1;i<=N;i++)
	{
		delta[1][i] = PI[i]*B[i][O[1]];
		psi[1][i] = 0;
	}
	for(int j=1;j<=N;j++)
	{
		for(int t=2;t<=T;t++)
		{
			long double max = INT_MIN;
			int index = 0;
			for(int i=1;i<=N;i++)
			{
				if(delta[t-1][i]*A[i][j] > max)
				{
					max = delta[t-1][i]*A[i][j];
					index = i;
				}
			}
			delta[t][j] = max*B[j][O[t]];
			psi[t][j] = index;
		}
	}
	
	// Calculating P(star) and finding state sequence
	
	long double P_s = INT_MIN;
	for(int i=1;i<=N;i++)
	{
		if(delta[T][i]>P_s)
		{
			P_s = delta[T][i];
			Q[T] = i;
		}
	}
	//printf("\n\nValues of delta - \n");
	//for(int t=1;t<=T;t++)
	//{
	//	for(int i=1;i<=N;i++)
	//	{
	//		printf("%24.17E ", delta[t][i]);
	//	}
	//	printf("\n");
	//}
	//printf("\n\nValues of psi - \n");
	//for(int t=1;t<=T;t++)
	//{
	//	for(int i=1;i<=N;i++)
	//	{
	//		printf(" %1.Lf\t\t", psi[t][i]);
	//	}
	//	printf("\n");
	//}
	for(int t=T-1; t>0;t--)
	{
		Q[t] = psi[t+1][Q[t+1]];
	}
	//printf("\n\n");
	//printf("\t\t\t\t\t\t****************************************************************************************************\n\n");
	//printf("\t\t\t\t\t\t\t\t\t\tP (star) : %E\n\n", P_s);
	//printf("\t\t\t\t\t\t****************************************************************************************************\n");
	//printf("\n\nState Sequence :\n");
	//for(int i=1;i<=T;i++)
	//{
	//	printf("%ld ", Q[i]);
	//}
}

// Function to solve problem 3
void problem3()
{
	// Calculating values of Xi
	
	for(int t=1; t<=T-1; t++)
	{
		long double sum = 0;
		for(int i=1; i<=N; i++)
		{
			for(int j=1; j<=N; j++)
			{
				sum += alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j];
			}
		}
		for(int i=1; i<=N; i++)
		{
			for(int j=1; j<=N; j++)
			{
				xi[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j]) / sum;
			}
		}
	}
	//printf("\n\nValues of xi - \n");
	//for(int t=1; t<=T-1; t++)
	//{
	//	for(int i=1; i<=N; i++)
	//	{
	//		for(int j=1; j<=N; j++)
	//		{
	//			printf("%24.17E ", xi[t][i][j]);
	//		}
	//		printf("\n");
	//	}
	//	printf("\n");
	//}
	
	// Calculating gamma using Xi
	
	//printf("\n\nValues of gamma - \n");
	for(int t=1; t<=T-1; t++)
	{
		for(int i=1; i<=N; i++)
		{
			long double sum = 0;
			for(int j=1; j<=N; j++)
			{
				sum += xi[t][i][j];
			}
			gamma[t][i] = sum;
			//printf("%24.17E ", gamma[t][i]);
		}
		//printf("\n");
	}

	// Calculating Lamba_bar

	//printf("\n\nNew Pi - \n");
	for(int i=1; i<=N; i++)		// PI_BAR
	{
		PI[i] = gamma[1][i];
	//	printf("%24.17E ", PI[i]);
	}

	//printf("\n\nNew A - \n");
	for(int i=1; i<=N; i++)		// A_BAR
	{
		long int index = 0;
		long double max = INT_MIN;
		long double sum = 0;
		for(int j=1; j<=N; j++)
		{
			long double sumXi = 0, sumG = 0;
			for(int t=1; t<=T-1; t++)
			{
				sumXi += xi[t][i][j];
				sumG += gamma[t][i];
			}
			A[i][j] = sumXi / sumG;
			//printf("%24.17E ", A[i][j]);
			if(A[i][j] < pow(10.0, -30.0))
			{
				A[i][j] = pow(10.0, -30.0);
			}
			sum += A[i][j];
			if(A[i][j] > max)
			{
				index = j;
				max = A[i][j];
			}
		}
		if(sum > 1 || sum < 1)	// Making sure that sum of row is 1
		{
			A[i][index] = A[i][index] + 1 - sum;
		}
		//printf("\n");
	}

	long double sum[N+1] = {0};
	//printf("\n\nNew B - \n");
	for(int j=1; j<=N; j++)		// B_BAR
	{ 
		long int index = 0;
		long double max = INT_MIN;
		for(int k=1; k<=M; k++)
		{
			long double sumG=0, sumGK=0;
			for(int t=1; t<=T; t++)
			{
				sumG += gamma[t][j];
				if(k == O[t])
				//if(k == O[t] && Q[t] == j)
				{
					sumGK += gamma[t][j];
				}
			}
			B[j][k] = sumGK / sumG;
			if(B[j][k] < pow(10.0, -30.0))
			{
				B[j][k] = pow(10.0, -30.0);
			}
			sum[j] += B[j][k];
			if(B[j][k] > max)
			{
				max = B[j][k];
				index = k;
			}
		}
		if(sum[j] > 1 || sum[j] < 1)		// Making sure that sum of row is 1
		{
			B[j][index] = B[j][index] + 1 - sum[j];
		}
	}
	for(int j=1; j<=N; j++)
	{
		for(int k=1; k<=M; k++)
		{
			//printf("%24.17E ", B[j][k]);
		}
		//printf("\n");
	}
}

void reset()
{
	for(int i = 0; i <= maxT; i++)
	{
		for(int j = 0; j <= N; j++)
		{
			alpha[i][j] = beta[i][j] = gamma[i][j] = delta[i][j] = psi[i][j] = 0;
		}
	}
}

void resetNew()
{
	for(int i = 1; i <= N; i++)
	{
		for(int j = 1; j <= N; j++)
		{
			new_A[i][j] = 0;
		}
	}
	for(int i = 1; i <= N; i++)
	{
		for(int j = 1; j <= M; j++)
		{
			new_B[i][j] = 0;
		}
	}
	for(int j = 1; j <= N; j++)
	{
		new_PI[j] = 0;
	}
}

void training()
{
	char fname[25];
	char move[5] = {'L', 'R', 'U', 'D'};
	char folder[5] = {'S', 'A'};

	for(int i = 0; i < 4; i++)
	{
		FILE* myfile;
		resetNew();
		for(int k = 0; k < 2; k++)
		{
			for(int j = 1; j <= 20; j++)
			{
				sprintf(fname, "%c/%c%d.txt", folder[k], move[i], j);
				printf(fname);
				readFileHMM();
				Normalize(fname); // Normalize data
				RiAiCi(); // Calculate Ri, Ai and Ci 
				tokhuraDistance();
				//for(int k = 0; k < T; k++)
				//{
				//	printf("%d ", O[k]);
				//}
				//printf("%d\n", T);
				for(int it = 0; it < 10; it++)
				{
					//printf("%d %d ", j, it);
					reset();
					long double val = problem1();
					problem2();
					problem3();
					//printf("%24.17E -> ", val);
				}
				printf("\n");
				//char s[50];
				//sprintf(s, "Log/224101017_E_%d_%d.txt", i, j);
				//myfile = fopen(s, "w");
				//fprintf(myfile, "A matrix\n");
				for(int l = 1; l <= N; l++)
				{
					for(int r = 1; r <= N; r++)
					{
						new_A[l][r] += A[l][r];
						//fprintf(myfile, "%24.17E ", A[l][r]);
					}
					//fprintf(myfile, "\n");
				}
				//fprintf(myfile, "\n");
				//fprintf(myfile, "B matrix\n");
				for(int l = 1; l <= N; l++)
				{
					for(int r = 1; r <= M; r++)
					{
						new_B[l][r] += B[l][r];
						//fprintf(myfile, "%24.17E ", B[l][r]);
					}
					//fprintf(myfile, "\n");
				}
				//fclose(myfile);
			}
		}
		//printf("\n");
		char s[50];
		sprintf(s, "FinalModel/%c/A.txt", move[i]);
		myfile = fopen(s, "w");
		for(int l = 1; l <= N; l++)
		{
			for(int r = 1; r <= N; r++)
			{
				fprintf(myfile, "%24.17E ", new_A[l][r] / 40);
			}
			fprintf(myfile, "\n");
		}
		fclose(myfile);
		sprintf(s, "FinalModel/%c/B.txt", move[i]);
		myfile = fopen(s, "w");
		for(int l = 1; l <= N; l++)
		{
			for(int r = 1; r <= M; r++)
			{
				fprintf(myfile, "%24.17E ", new_B[l][r] / 40);
			}
			fprintf(myfile, "\n");
		}
		fclose(myfile);
	}
}

void readImprove(char ch)
{
	FILE* myfile;
	char fname[30];
	sprintf(fname, "FinalModel/%c/B.txt", ch);
	myfile = fopen(fname, "r");
	int i = 1, j = 1;
	long double inp;
	while(fscanf(myfile, "%Lf", &inp)!=EOF)
	{
		B[i][j] = inp;
		//printf("%24.17Lf ", B[i][j]);
		j++;
		if(j > M)
		{
			j = 1;
			i++;
			//printf("\n");
		}
	}
	fclose(myfile);
	sprintf(fname, "FinalModel/%c/A.txt", ch);
	myfile = fopen(fname, "r");
	i = 1, j = 1;
	//printf("\n");
	//printf("\n");
	while(fscanf(myfile, "%Lf", &inp)!=EOF)
	{
		A[i][j] = inp;
		//printf("%24.17Lf ", A[i][j]);
		j++;
		if(j > N)
		{
			j = 1;
			i++;
			//printf("\n");
		}
	}
	fclose(myfile);
	myfile = fopen("InitialModel/PI.txt", "r");
	i = 1;
	//printf("\n");
	//printf("\n");
	while(fscanf(myfile, "%Lf", &inp)!=EOF)
	{
		PI[i] = inp;
		//printf("%24.17Lf ", PI[i]);
		i++;
		//printf("\n");
	}
	fclose(myfile);
}

void improveModel()
{
	char fname[25];
	char move[5] = {'L', 'R', 'U','D'};
	char file[2] = {'A','S'};
	for(int i = 0; i < 4; i++)
	{
		for(int k=0;k<2;k++)
		{
		FILE* myfile;
		resetNew();
		for(int j = 1; j <= 20; j++)
		{
			sprintf(fname, "%c/%c%d.txt", file[k], move[i], j);
			printf(fname);
			readImprove(move[i]);
			Normalize(fname); // Normalize data
			RiAiCi(); // Calculate Ri, Ai and Ci 
			tokhuraDistance();
			//for(int k = 0; k < T; k++)
			//{
			//	printf("%d ", O[k]);
			//}
			//printf("%d\n", T);
			for(int it = 0; it < 20; it++)
			{
				//printf("%d %d ", j, it);
				reset();
				long double val = problem1();
				problem2();
				problem3();
				//printf("%24.17E -> ", val);
			}
			printf("\n");
			//char s[50];
			//sprintf(s, "Log/224101017_E_%d_%d.txt", i, j);
			//myfile = fopen(s, "w");
			//fprintf(myfile, "A matrix\n");
			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= N; r++)
				{
					new_A[l][r] += A[l][r];
					//fprintf(myfile, "%24.17E ", A[l][r]);
				}
				//fprintf(myfile, "\n");
			}
			//fprintf(myfile, "\n");
			//fprintf(myfile, "B matrix\n");
			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= M; r++)
				{
					new_B[l][r] += B[l][r];
					//fprintf(myfile, "%24.17E ", B[l][r]);
				}
				//fprintf(myfile, "\n");
			}
			//fclose(myfile);
		}
		//printf("\n");
		char s[50];
		sprintf(s, "FinalModel/%c/A.txt", move[i]);
		myfile = fopen(s, "w");
		for(int l = 1; l <= N; l++)
		{
			for(int r = 1; r <= N; r++)
			{
				fprintf(myfile, "%24.17E ", new_A[l][r] / 20);
			}
			fprintf(myfile, "\n");
		}
		fclose(myfile);
		sprintf(s, "FinalModel/%c/B.txt", move[i]);
		myfile = fopen(s, "w");
		for(int l = 1; l <= N; l++)
		{
			for(int r = 1; r <= M; r++)
			{
				fprintf(myfile, "%24.17E ", new_B[l][r] / 20);
			}
			fprintf(myfile, "\n");
		}

		fclose(myfile);
		}
	}
}

void testing()
{
	char fname[25];
	char move[5] = {'L', 'R', 'U', 'D'};
	int count = 0;
	for(int l = 0; l < 4; l++)
	{
		FILE* myfile;
		resetNew();
		int countD = 0;
		for(int m = 6; m <= 15; m++)
		{
			char test[5] = {'L', 'R', 'U', 'D'};
			sprintf(fname, "S/%c%d.txt", move[l], m);
			Normalize(fname); // Normalize data
			RiAiCi(); // Calculate Ri, Ai and Ci 
			tokhuraDistance();
			long double max = INT_MIN, val = 0;
			int index = -1;
			for(int k = 0; k < 4; k++)
			{
				char s[50];
				sprintf(s, "FinalModel/%c/B.txt", test[k]);
				myfile = fopen(s, "r");
				int i = 1, j = 1;
				long double inp;
				while(fscanf(myfile, "%Lf", &inp)!=EOF)
				{
					B[i][j] = inp;
					//printf("%24.17E ", B[i][j]);
					j++;
					if(j > M)
					{
						j = 1;
						i++;
						//printf("\n");
					}
				}
				fclose(myfile);
				sprintf(s, "FinalModel/%c/A.txt", test[k]);
				myfile = fopen(s, "r");
				i = 1, j = 1;
				//printf("\n");
				//printf("\n");
				while(fscanf(myfile, "%Lf", &inp)!=EOF)
				{
					A[i][j] = inp;
					//printf("%24.17E ", A[i][j]);
					j++;
					if(j > N)
					{
						j = 1;
						i++;
						//printf("\n");
					}
				}
				fclose(myfile);
				myfile = fopen("InitialModel/PI.txt", "r");
				i = 1;
				//printf("\n");
				//printf("\n");
				while(fscanf(myfile, "%Lf", &inp)!=EOF)
				{
					PI[i] = inp;
					//printf("%24.17Lf ", PI[i]);
					i++;
					//printf("\n");
				}	
				fclose(myfile);
				val = problem1();
				if(val >= max)
				{
					max = val;
					index = k;
				}
				printf("%c -> %24.17E\n", test[k], val);
			}
			printf("\n%c%d -> %c\n", move[l], m, move[index]);
		}
	}
}

char liveTesting()
{
	system("Recording_Module.exe 2 input_file.wav input_file.txt  >nul 2>nul");
	Normalize("input_file.txt"); // Normalize data
	RiAiCi(); // Calculate Ri, Ai and Ci 
	tokhuraDistance();
	char move[4] = {'L', 'R', 'U', 'D'};
	FILE* myfile;
	long double max = INT_MIN, val = 0;
	int index = -1;
	for(int k = 0; k < 4; k++)
	{
		char s[50];
		sprintf(s, "FinalModel/%c/B.txt", move[k]);
		myfile = fopen(s, "r");
		int i = 1, j = 1;
		long double inp;
		while(fscanf(myfile, "%Lf", &inp)!=EOF)
		{
			B[i][j] = inp;
			//printf("%24.17E ", B[i][j]);
			j++;
			if(j > M)
			{
				j = 1;
				i++;
				//printf("\n");
			}
		}
		fclose(myfile);
		sprintf(s, "FinalModel/%c/A.txt", move[k]);
		myfile = fopen(s, "r");
		i = 1, j = 1;
		//printf("\n");
		//printf("\n");
		while(fscanf(myfile, "%Lf", &inp)!=EOF)
		{
			A[i][j] = inp;
			//printf("%24.17E ", A[i][j]);
			j++;
			if(j > N)
			{
				j = 1;
				i++;
				//printf("\n");
			}
		}
		fclose(myfile);
		myfile = fopen("InitialModel/PI.txt", "r");
		i = 1;
		//printf("\n");
		//printf("\n");
		while(fscanf(myfile, "%Lf", &inp)!=EOF)
		{
			PI[i] = inp;
			//printf("%24.17Lf ", PI[i]);
			i++;
			//printf("\n");
		}	
		fclose(myfile);
		val = problem1();
		//printf("%c -> %24.13E\n", move[k], val);
		if(val >= max)
		{
			max = val;
			index = k;
		}
	}
	//printf("\nMove is -> %c", move[index]);
	return move[index] ;
}

void liveTraining()
{
	char name[100];
	int u_i = 1 , d_i = 1, l_i = 1, r_i = 1;
	int ans;
	char txtfileName[200], wavfileName[200], command[300];
	long double data;
	int word_no;
	char fileA[50], fileB[50];
	int flag = 0;

	printf("\n\tWelcome to live training.\n");

	printf("\nEnter your first name : ");

	scanf("%s", &name);

	FILE *fpra;
	FILE *fprb;
	FILE *fpwa;
	FILE *fpwb;

	do
	{
		printf("\nType index of word you want to record : \n");
		printf(" 1.Up \n 2.Down \n 3.Left \n 4.Right \n");
		
		scanf("%d", &word_no);

		flag = 0;

		if(word_no == 1)
		{
			printf("\nSay \"Up\"\n\n");

			sprintf(txtfileName, "LiveTrain/%s_U%d.txt", name, u_i);
			sprintf(wavfileName, "LiveTrain/%s_U%d.wav", name, u_i);

			sprintf(command, "Recording_Module.exe 3 %s %s", wavfileName, txtfileName);

			u_i++;

			sprintf(fileA, "FinalModel/U/A.txt");
			sprintf(fileB, "FinalModel/U/B.txt");
		}
		else if(word_no == 2)
		{
			printf("\nSay \"Down\"\n\n");

			sprintf(txtfileName, "LiveTrain/%s_D%d.txt", name, u_i);
			sprintf(wavfileName, "LiveTrain/%s_D%d.wav", name, u_i);

			sprintf(command, "Recording_Module.exe 3 %s %s", wavfileName, txtfileName);

			d_i++;

			sprintf(fileA, "FinalModel/D/A.txt");
			sprintf(fileB, "FinalModel/D/B.txt");
		}
		else if(word_no == 3)
		{
			printf("\nSay \"Left\"\n\n");

			sprintf(txtfileName, "LiveTrain/%s_L%d.txt", name, u_i);
			sprintf(wavfileName, "LiveTrain/%s_L%d.wav", name, u_i);

			sprintf(command, "Recording_Module.exe 3 %s %s", wavfileName, txtfileName);

			l_i++;

			sprintf(fileA, "FinalModel/L/A.txt");
			sprintf(fileB, "FinalModel/L/B.txt");
		}
		else if(word_no == 4)
		{
			printf("\nSay \"Right\"\n\n");

			sprintf(txtfileName, "LiveTrain/%s_R%d.txt", name, u_i);
			sprintf(wavfileName, "LiveTrain/%s_R%d.wav", name, u_i);

			sprintf(command, "Recording_Module.exe 3 %s %s", wavfileName, txtfileName);

			r_i++;

			sprintf(fileA, "FinalModel/R/A.txt");
			sprintf(fileB, "FinalModel/R/B.txt");
		}
		else
		{
			printf("\nPlease enter a valid choice\n");
			flag = 1;
		}

		if(flag == 0)
		{
			system(command);
		
			readFileHMM();
			Normalize(txtfileName); // Normalize data
			RiAiCi(); // Calculate Ri, Ai and Ci 
			tokhuraDistance();

			for(int it = 0; it < 10; it++)
			{
				reset();
				long double val = problem1();
				problem2();
				problem3();
				//printf("%24.17E -> ", val);
			}
			printf("\n");

			//read final model of A and B
			fpra = fopen(fileA, "r");

			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= N; r++)
				{
					fscanf(fpra, "%Lf", &new_A[l][r]);
					new_A[l][r] *= 0.75;
					//printf("%Lf", new_A[l][r]);
				}
			}

			fprb = fopen(fileB, "r");
			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= M; r++)
				{
					fscanf(fprb, "%Lf", &new_B[l][r]);
					new_B[l][r] *= 0.75;
					//printf("%Lf", new_B[l][r]);
				}
			}

			//modify final models to include new training data
			
			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= N; r++)
				{
					new_A[l][r] += 0.25 * A[l][r];	
					//printf("%24.17e", A[l][r]);
				}
			}
				
			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= M; r++)
				{
					new_B[l][r] += 0.25 * B[l][r];
					//printf("%24.17e", B[l][r]);
				}
			}

			//store modified model back in files
		
			fpwa = fopen(fileA, "w");
		
			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= N; r++)
				{
					fprintf(fpwa, "%24.17E ", new_A[l][r]);
				}
				fprintf(fpwa, "\n");
			}
			fclose(fpwa);

			fpwb = fopen(fileB, "w");
			for(int l = 1; l <= N; l++)
			{
				for(int r = 1; r <= M; r++)
				{
					fprintf(fpwb, "%24.17E ", new_B[l][r]);
				}
				fprintf(fpwb, "\n");
			}
			fclose(fpwb);

		}
		printf("\nDo you want to record again ?\n 1.Yes \n 2.No : ");

		scanf("%d", &ans);
	}while(ans == 1 || ans == 0);


}

int liveNTesting()
{
	system("Recording_Module.exe 2 input_file.wav input_file.txt >nul 2>nul");
	Normalize("input_file.txt"); // Normalize data
	RiAiCi(); // Calculate Ri, Ai and Ci 
	tokhuraDistanceDigit();
	FILE* myfile;
	long double max = INT_MIN, val = 0;
	int index = -1;
	for(int k = 1; k <= 3; k++)
	{
		char s[50];
		sprintf(s, "FM/%d/B.txt", k);
		myfile = fopen(s, "r");
		int i = 1, j = 1;
		long double inp;
		while(fscanf(myfile, "%Lf", &inp)!=EOF)
		{
			B[i][j] = inp;
			//printf("%24.17E ", B[i][j]);
			j++;
			if(j > M)
			{
				j = 1;
				i++;
				//printf("\n");
			}
		}
		fclose(myfile);
		sprintf(s, "FM/%d/A.txt", k);
		myfile = fopen(s, "r");
		i = 1, j = 1;
		//printf("\n");
		//printf("\n");
		while(fscanf(myfile, "%Lf", &inp)!=EOF)
		{
			A[i][j] = inp;
			//printf("%24.17E ", A[i][j]);
			j++;
			if(j > N)
			{
				j = 1;
				i++;
				//printf("\n");
			}
		}
		fclose(myfile);
		myfile = fopen("InitialModel/PI.txt", "r");
		i = 1;
		//printf("\n");
		//printf("\n");
		while(fscanf(myfile, "%Lf", &inp)!=EOF)
		{
			PI[i] = inp;
			//printf("%24.17Lf ", PI[i]);
			i++;
			//printf("\n");
		}	
		fclose(myfile);
		val = problem1();
		//printf("%d -> %24.13E\n", k, val);
		if(val >= max)
		{
			max = val;
			index = k;
		}
	}
	//printf("\nDigit is -> %d", index);
	return index;
}


void game()
{
	system("cls");
	char c ='Y';
	char s[50];
	FILE * myfile;
	sprintf(s, "moves.txt");
	myfile = fopen(s, "w");
				
	// enumeration (an enum is a set of named integer constants)
	// MazeObject - user (custom) data type
	enum MazeObject { HALL, WALL, COIN, ENEMY, BORDER , GATE };
	enum Color { DARKGREEN = 2, YELLOW = 14, RED = 12, BLUE = 9, WHITE = 15, DARKYELLOW = 6, DARKRED = 4 };
	enum KeyCode { ENTER = 13, ESCAPE = 27, SPACE = 32, LEFT = 75, RIGHT = 77, UP = 72, DOWN = 80 };

	HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);

	// console font settings
	CONSOLE_FONT_INFOEX font; // https://docs.microsoft.com/en-us/windows/console/console-font-infoex
	font.cbSize = sizeof(font);
	font.dwFontSize.Y = 50;
	font.FontFamily = FF_DONTCARE;
	font.FontWeight = FW_NORMAL;
	wcscpy_s(font.FaceName, 9, L"Consolas");
	SetCurrentConsoleFontEx(h, 0, &font);

	// hiding the blinking cursor
	CONSOLE_CURSOR_INFO cursor;
	cursor.bVisible = false; // hide cursor
	cursor.dwSize = 1; 
	SetConsoleCursorInfo(h, &cursor);

	system("title Maze");
	MoveWindow(GetConsoleWindow(), 0, 0, 1920, 1080, true);
	//0 - padding to the left of the left border of the desktop to the left border of the console window (in pixels!)
	//0 - top margin from the top border of the desktop to the top border of the console window
	//1920 - console window width in pixels
	//1080 - console window height
	//true - redraw the window after moving
	MessageBoxA(0, "FIND THE EXIT TO WIN \rSPEAK LEFT TO MOVE LEFT \rSPEAK RIGHT TO MOVE RIGHT \rSPEAK UP TO MOVE UP \rSPEAK DOWN TO MOVE DOWN","GOOD LUCK", 0);

	srand(time(0));

	const int WIDTH = 40; // maze width
	const int HEIGHT = 10; // maze height

	int maze[HEIGHT][WIDTH] = {}; // maze

	// array filling algorithm
	for (int y = 0; y < HEIGHT; y++) //line enumeration
	{
		for (int x = 0; x < WIDTH; x++) //iterate over columns
		{
			maze[y][x] = rand() % 4; // 4 types of objects in the game

			if (maze[y][x] == MazeObject::ENEMY) // if an enemy is generated in the labyrinth
			{
				int probability = rand() % 15; // 0...14, if 0 is rolled - the enemy will remain, only 1/15 of the enemies will remain
				if (probability != 0) // remove the enemy
				{
					maze[y][x] = MazeObject::HALL; // put a hall in place of the enemy
				}
			}

			if (maze[y][x] == MazeObject::WALL) // if a wall is generated in the maze
			{
				int probability = rand() % 2; // 0...1, if 0 is rolled, the wall will remain, only half of the walls will remain
				if (probability != 0) // remove the wall
				{
					maze[y][x] = MazeObject::HALL; // put a hall in place of the wall
				}
			}

			if (x == 0 || y == 0 || x == WIDTH - 1 || y == HEIGHT - 1) maze[y][x] = MazeObject::BORDER; // white frame

			if (x == 0 && y == 1 || x == 1 && y == 2 || x == 2 && y == 2) maze[y][x] = MazeObject::HALL; // Entry

			if (x == WIDTH - 1 && y == HEIGHT - 3 ||
				x == WIDTH - 2 && y == HEIGHT - 3 ||
				x == WIDTH - 3 && y == HEIGHT - 3) maze[y][x] = MazeObject::GATE; // Exit

		}
	}

	// labyrinth display
	for (int y = 0; y < HEIGHT; y++) // line enumeration
	{
		for (int x = 0; x < WIDTH; x++) // iterate over columns
		{
			switch (maze[y][x])
			{
			case MazeObject::HALL: // hall 
				cout << " ";
				break;

			case MazeObject::WALL: // wall
				SetConsoleTextAttribute(h, Color::DARKGREEN);
				cout << (char)178;
				break;

			case MazeObject::BORDER: // border
				SetConsoleTextAttribute(h, Color::WHITE);
				cout << (char)178;
				break;

			case MazeObject::GATE: // border
				SetConsoleTextAttribute(h, Color::YELLOW);
				cout << (char)178;
				break;

			case MazeObject::COIN: // coin
				SetConsoleTextAttribute(h, Color::YELLOW);
				cout << ".";
				break;
			
			
			}
		}
		cout << "\n";
	}

	/////////////////////////////////////////////////////////////////////

	// placement of the main character (GG)
	COORD position; // our character's coordinates
	position.X = 0;
	position.Y = 2;
	SetConsoleCursorPosition(h, position);
	SetConsoleTextAttribute(h, Color::BLUE);
	cout << (char)2;

	int coins = 0; // collected coins counter
	int health = 100; // the number of health points of the main character

	/////////////////////////////////////////////////////////////////////
	// information on all indicators
	
	COORD infobox;
	infobox.X = WIDTH + 1;
	infobox.Y = 1;
	SetConsoleCursorPosition(h, infobox);
	SetConsoleTextAttribute(h, Color::DARKYELLOW);
	cout << "COINS: ";
	SetConsoleTextAttribute(h, Color::YELLOW);
	cout << coins << "\n"; // 0
	/*
	infobox.Y = 2;
	SetConsoleCursorPosition(h, infobox);
	SetConsoleTextAttribute(h, Color::DARKRED);
	cout << "HEALTH: ";
	SetConsoleTextAttribute(h, Color::RED);
	cout << health << "\n";
	*/
	while (true)
	{
		MessageBoxA(0, "Speak NOW","GOOD LUCK", 0);
		c = liveTesting();
		char buff[100];
		int i = 0;
		string name ;
		name += c;
		int answer = IDNO ;
		while( i == 0 )
		{
			sprintf_s(buff,"You said %s , do you want to move 1-3 steps in this direction",name.c_str());
			int answer = MessageBoxA(0, buff ,"GOOD LUCK", MB_YESNO);
			if (answer == IDNO)
			{
				MessageBoxA(0, "sorry to hear that , Speak again","GOOD LUCK", 0);
				c = liveTesting();
				name[0] = c;
			}
			else
				i = 1;
		}

		i=0;
		int p = 0;
		MessageBoxA(0, "Speak how many steps (1-3) you want to move","GOOD LUCK", 0);
		i = liveNTesting();
		c = '0' + i;
		name += c;

		while( p == 0 )
		{
			sprintf_s(buff,"You said %s ",name.c_str());
			int answer = MessageBoxA(0, buff ,"GOOD LUCK", MB_YESNO);
			if (answer == IDNO)
			{
				MessageBoxA(0, "sorry to hear that , Speak again","GOOD LUCK", 0);
				i = liveNTesting();
				c = '0' + i;
				name[1] = c;
			}
			else
				p = 1;
		}

		fprintf(myfile, "%c ", c);
		if (_kbhit() || c) // if there was a keystroke by the user
		{
			int code = _getch(); //get character, getting the code of the pressed key
			if (code == 224) { // if it's an arrow
				code = _getch(); // get specific arrow code
			}
															// erasing the character in the old position
			SetConsoleCursorPosition(h, position);
			cout << " ";
			if ((code == KeyCode::LEFT || name[0] == 'L'))
				while( i > 0 && maze[position.Y][position.X - 1] != MazeObject::WALL && maze[position.Y][position.X - 1] != MazeObject::BORDER)
				{
					position.X--;
					i--;
					if (maze[position.Y][position.X] == MazeObject::COIN)
					{
						coins++; // collected more by one coin
						infobox.X = WIDTH + 1;
						infobox.Y = 1;
						SetConsoleCursorPosition(h, infobox);
						SetConsoleTextAttribute(h, Color::DARKYELLOW);
						cout << "COINS: ";
						SetConsoleTextAttribute(h, Color::YELLOW);
						cout << coins << "\n";
						maze[position.Y][position.X] = MazeObject::HALL; // remove the coin from the maze
					}
				}

			else if ((code == KeyCode::RIGHT || name[0] == 'R'))
				while( i > 0 && maze[position.Y][position.X + 1] != MazeObject::WALL && maze[position.Y][position.X + 1] != MazeObject::BORDER)
				// and at the same time in the maze on the same line (where GG) and
				// a little (one cell) to the right, 1 column from GG
				{
					position.X++;
					i--;
					if (maze[position.Y][position.X] == MazeObject::COIN)
					{
						coins++; // collected more by one coin
						infobox.X = WIDTH + 1;
						infobox.Y = 1;
						SetConsoleCursorPosition(h, infobox);
						SetConsoleTextAttribute(h, Color::DARKYELLOW);
						cout << "COINS: ";
						SetConsoleTextAttribute(h, Color::YELLOW);
						cout << coins << "\n";
						maze[position.Y][position.X] = MazeObject::HALL; // remove the coin from the maze
					}
				}
			else if ((code == KeyCode::UP || name[0] == 'U'))
				while( i > 0  && maze[position.Y - 1][position.X] != MazeObject::WALL && maze[position.Y - 1][position.X] != MazeObject::BORDER)
				{
					position.Y--;
					i--;
					if (maze[position.Y][position.X] == MazeObject::COIN)
					{
						coins++; // collected more by one coin
						infobox.X = WIDTH + 1;
						infobox.Y = 1;
						SetConsoleCursorPosition(h, infobox);
						SetConsoleTextAttribute(h, Color::DARKYELLOW);
						cout << "COINS: ";
						SetConsoleTextAttribute(h, Color::YELLOW);
						cout << coins << "\n";
						maze[position.Y][position.X] = MazeObject::HALL; // remove the coin from the maze
					}
				}
			else if ((code == KeyCode::DOWN || name[0] == 'D'))
				while(i > 0 && maze[position.Y + 1][position.X] != MazeObject::WALL && maze[position.Y + 1][position.X] != MazeObject::BORDER)
				{
					position.Y++;
					i--;
					if (maze[position.Y][position.X] == MazeObject::COIN)
					{
						coins++; // collected more by one coin
						infobox.X = WIDTH + 1;
						infobox.Y = 1;
						SetConsoleCursorPosition(h, infobox);
						SetConsoleTextAttribute(h, Color::DARKYELLOW);
						cout << "COINS: ";
						SetConsoleTextAttribute(h, Color::YELLOW);
						cout << coins << "\n";
						maze[position.Y][position.X] = MazeObject::HALL; // remove the coin from the maze
					}
				}

			// showing GG in a new position
			SetConsoleCursorPosition(h, position);
			SetConsoleTextAttribute(h, Color::BLUE);
			cout << (char)2;

			////////////////////////////////////////////////////////////////
			// intersection with array elements
			if (position.Y == HEIGHT - 3 &&
				position.X == WIDTH - 1)
			{
				MessageBoxA(0, "you find exit!", "WIN!", 0);
				system("cls");
				exit(0);
				game(); // to start first, but on a different random location
			}

			// crossing with coins
			// if there is a coin in the labyrinth at the GG position (under it)
			
			if (maze[position.Y][position.X] == MazeObject::COIN)
			{
				coins++; // collected more by one coin
				infobox.X = WIDTH + 1;
				infobox.Y = 1;
				SetConsoleCursorPosition(h, infobox);
				SetConsoleTextAttribute(h, Color::DARKYELLOW);
				cout << "COINS: ";
				SetConsoleTextAttribute(h, Color::YELLOW);
				cout << coins << "\n";
				maze[position.Y][position.X] = MazeObject::HALL; // remove the coin from the maze
				}
		}
	}
	fclose(myfile);
	return;
}


int menu()
{
    system("cls");
    string Menu[3] = { "New User", "Play Game", "Exit" };
    int pointer = 0;
    bool bMainMenu = true;

    while (bMainMenu)
    {
        system("cls");

        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 15);
        cout << "Main Menu\n\n";

        for (int i = 0; i < 3; ++i)
        {
            if (i == pointer)
            {
                SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 11);
                cout << Menu[i] << endl; 
            }
            else
            {
                SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 15);
                cout << Menu[i] << endl ;
            }
        }

        while (bMainMenu)
        {
            if (GetAsyncKeyState(VK_UP)&1)
            {
                pointer = pointer - 1;
                if (pointer == -1)
                {
                    pointer = 2;
                }
                break;
            }
            else if (GetAsyncKeyState(VK_DOWN)&1)
            {
                pointer += 1;
                if (pointer == 3)
                {
                    pointer = 0;
                }
                break;
            }
            else if (GetAsyncKeyState(VK_RETURN)&1)
            {
                switch (pointer)
                {
                case 0:
					{
						system("cls");
						liveTraining();
						Sleep(1000);
						bMainMenu = false;
					}
                case 1:
					{
						system("cls");
						game();
					   // Sleep(1000);
						bMainMenu = false;
						break;
					}
                case 2:
					{
						//thank_you();
						system("cls");
						std::cout << "Thank you , have a Good Day \n";
						Sleep(1000);
						bMainMenu = false;
						break;
					}
                default:
					{
						cout << "Invalid Input! ";
					}
                }
            }
        }
        Sleep(150);
    }
    return 0;
}


int _tmain(int argc, _TCHAR* argv[])
{
	int choice;
	menu();
	/*
	scanf("%d",&choice);
	switch (choice)
	{
		case 1: liveTraining();
				//break;

		case 2: system("cls");
				game();
				break ;
	}
	//createUniverse();
	//LBG();
	//training();
	/*for(int i = 0; i < 15; i++)
	{
		improveModel();
	}
	*/
	//testing();
	//liveTesting();
	//liveTraining();
	return 0;
}

