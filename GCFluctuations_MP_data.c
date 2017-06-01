/****************************************************************************************************************
 %    FLUCTUATIONS IN PERCOLATION OF SPARSE COMPLEX NETWORKS
 %
 %
 % This code implements the message passing algorithm to evaluate fluctuation in percolation of sparse complex networks
 % and compairs the result to brute force average over different realizations of the initial damage.
 % 
 % This code can be redistributed and/or modified
 % under the terms of the GNU General Public License as published by
 % the Free Software Foundation, either version 3 of the License, or (at
 % your option) any later version.
 %  
 % This program is distributed ny the authors in the hope that it will be 
 % useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 %
 %  
 % If you use this code please cite the following  paper:
 %
 % G. Bianconi, arXiv preprint, arXiv:1703.05528 (2017) 
 %
 % (c) G. Bianconi (email: ginestra.bianconi@gmail.com ) 
 %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 Using the notation of the paper our code used node-independent probabilities p^{(1)}=p^{(2)}=p 
 that a node is not initially damaged in each of the two realization of the percolation problem
 and a probability 
 p^{[11]} that the node is not initially damaged in both realizations of the percolation problem
 Here we take
 p^{[11]}=p^a with a>=1 to be an imput parameter of the code
 
 ****************************************************************************************************************
 Input:
 1) filename  -name of the file including the edge list of weighted undirected network 
 the form of 
 node1	node2  weight 
 (weight not ot be used by the program).
 2) Nrunmax- total number of pair of random realization of the initial damage
 3) a- parameter determining the correlations between the two realizations of the initial disorder:
		p^{[11]}=p^a
		-a=2 uncorrelated realizations
		-a>=1 a<2 positive correlations
		-a>2      negative correlations
 
 Output: 
 1)file with name "GC_MP_a%f_filename" has five columns
 
 p S_{[10]} S_{[11]} S_{(1)} X p11
 
 p-probability that initially a node is not damaged in the first or the second realization of the damage
 S_{[10]}-fraction of nodes in the giant component of the first but not of the second realization of the percolation problem
 S_{[11]}-fraction of nodes in the giant component of both realizations of the percolation problem
 S_{(1)}-fraction of nodes in the giant component in the first realization of the percolation problem
 p11- probability that initially a node is not damaged in both realizations of the damage
 
  2)file with name "GC_sim_a%f_filename" indicates the same quantities in the same order 
 but this time obtained from Nrunmax realizations of pairs of initial damage.
 
 
 *********************************************************************************************************************/




#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#define err_tol 0.01

int *vis1,*size_cluster1,**knn1,*k1,c1,*occ,*dam1;


/**************************************************************************
 Recurrence is a subrutine for calulating the giant component in each layer
 **************************************************************************/
int Recurrence( int i , int cluster_size, int ncluster){
	int j, n3,aus_cluster_size;
	
	cluster_size++;		
		vis1[i]=ncluster;
		for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
			if((vis1[j]==0)&&(dam1[j]==1)){
				aus_cluster_size=Recurrence(j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
			}
	
		}
	
		return cluster_size;
}

int main(int argc, char** argv){
	int i,j,n,it,ncluster1, GMCC, cluster_size,m1,m1_aus,c1_aus,nrun,nc,np;
	int s1,Nc1,aus,aus3,**adj1,N,N0,**s,nrun2,w,Nrunmax,nc2;
	float p,x,*xd1,**n1,**n2,MGC1,MGC2,MGC3,MGC4,nsum1,nsumold1,aus1,aus2,sigma11,sigma10,**sigma11m,**sigma10m,p11,a;
	char filename[60],string1[60],string2[60];
	
	FILE *gp2,*gp,*ffp;
  

 
	srand48(time(NULL));

		
	
	printf("Insert file name\n");
	scanf("%s" , filename) ;
	ffp=fopen(filename, "r");
	printf("Insert Nrunmax\n");
	scanf("%d",&Nrunmax);
	printf("Insert a\n");
	scanf("%f",&a);
	sprintf(string1,"GC_MP_a%lf_%s",a, filename);
	sprintf(string2,"GC_sim_a%lf_%s",a, filename);
	gp=fopen(string1,"w");
	gp2=fopen(string2,"w");
	
	N=0;
	while(!feof(ffp)){
		fscanf(ffp,"%d %d %d", &i ,&j,&w);
		if(N<i+1) N=i+1;
		if(N<j+1)	N=j+1;
	}
	fclose(ffp);

	
	//Sm=(float*)calloc(200,sizeof(float));
	vis1=(int*)calloc(N,sizeof(int));
	
	occ=(int*)calloc(N,sizeof(int));
	k1=(int*)calloc(N,sizeof(int));
	xd1=(float*)calloc(N,sizeof(float));
	dam1=(int*)calloc(N,sizeof(int));
	knn1=(int**)calloc(N,sizeof(int*));
	n1=(float**)calloc(N,sizeof(float*));
	n2=(float**)calloc(N,sizeof(float*));

	adj1=(int**)calloc(N,sizeof(int*));
	s=(int**)calloc(N,sizeof(int*));
		for(i=0;i<N;i++){
			knn1[i]=(int*)calloc(N,sizeof(int));
			adj1[i]=(int*)calloc(N,sizeof(int));			
			n1[i]=(float*)calloc(N,sizeof(float));
			n2[i]=(float*)calloc(N,sizeof(float));
			s[i]=(int*)calloc(2,sizeof(int));
			
		}
	size_cluster1=(int*)calloc(N,sizeof(int));
	sigma11m=(float**)calloc(100,sizeof(float*));
	sigma10m=(float**)calloc(100,sizeof(float*));
	
	for(i=0;i<100;i++){
		sigma11m[i]=(float*)calloc(100,sizeof(float));
		sigma10m[i]=(float*)calloc(100,sizeof(float));
	}
	
	for(i=0;i<N;i++){
		k1[i]=0;
	}


	
	ffp=fopen(filename, "r");
	
	while(!feof(ffp)){
		fscanf(ffp,"%d %d %d", &i ,&j,&w);
				k1[i]++;
				k1[j]++;
				knn1[i][k1[i]-1]=j;
				knn1[j][k1[j]-1]=i;
				adj1[i][j]=1;
				adj1[j][i]=1;
	}
	N0++;
	for (i=0;i<N;i++){
		if(k1[i]>0){
			N0++;
		}
	}
			
	

		/***********************************************************************************
		 ***Message passing algorithm
		 ***********************************************************************************/


	for(p=1.;p>0.;p-=0.01){
		p11=pow(p,a);
		
		nsum1=0;
		if((2*p-p11)<1){
			for(i=0;i<N;i++){
				for(n=0;n<k1[i];n++){
					j=knn1[i][n];
					n1[i][j]=(drand48());
					n1[j][i]=(drand48());
					n2[i][j]=(drand48());
					n2[j][i]=(drand48());
					nsum1+=(n1[i][j]+n1[j][i]);
					nsum1+=(n2[i][j]+n2[j][i]);
				}
			}
			nsumold1=1000000;

			while(fabs(nsum1-nsumold1)>err_tol){	
				
				nsumold1=nsum1;
				nsum1=0;
				
				for (it=0;it<N;it++){
					i=(int)(N*drand48());
					for(n=0;n<k1[i];n++){
						j=knn1[i][n];
						aus1=1.;
						aus2=1.;
						for(np=0;np<k1[i];np++){
							if(np!=n){
								aus1=aus1*(1.-n1[knn1[i][np]][i]);
								aus2=aus2*(1.-2.*n1[knn1[i][np]][i]+n2[knn1[i][np]][i]);
							}
						}
						nsum1-=n1[i][j]+n2[i][j];
						
						n1[i][j]=p*(1-aus1);
						n2[i][j]=p11*(1-2*aus1+aus2);
						
						nsum1+=n1[i][j]+n2[i][j];						
					}
				}
			}
			MGC1=0;
			MGC2=0;
			MGC3=0;
			MGC4=0;
			for (i=0;i<N;i++){
				if(k1[i]>0){
					aus1=1;
					for(n=0;n<k1[i];n++){
						aus1=aus1*(1-n1[knn1[i][n]][i]);}
					
					aus2=1;
					for(n=0;n<k1[i];n++){
						aus2=aus2*(1-2*n1[knn1[i][n]][i]+n2[knn1[i][n]][i]);
					}
					
					MGC1+=p*(1-aus1);
					MGC2+=p11*(1-2*aus1+aus2);
					MGC3+=p*(1-aus1)-p11*(1-2*aus1+aus2);
					MGC4+=p11*(-aus1*aus1+aus2);
				}
			}
			
			
			fprintf(gp,"%lf %lf %lf %lf %lf \n",p,(float)MGC3/((float)N0),(float)MGC2/((float)N0),(float)(MGC1/((float)N0)),p11);
			
		}	
	}
	
	fclose(gp);
	printf("uno\n");
	getchar();
	for (nrun=0;nrun<Nrunmax;nrun++){	
		
		
		for(i=0;i<N;i++){
			xd1[i]=drand48();
		}
		
		for(nc=0;nc<100;nc++){
			
			p=0.01*nc;
			nc2=0;
			p11=pow(p,a);
			if((2*p-p11)<1){
				for (nrun2=0;nrun2<2;nrun2++){
					if(nrun2==0){
						for(i=0;i<N;i++){
							if(xd1[i]<p){
								dam1[i]=1;}
							else {dam1[i]=0;}
						}
					}
					if(nrun2==1){
						for(i=0;i<N;i++){
							if((xd1[i]>(p-p11))&&(xd1[i]<(2*p-p11))){
								dam1[i]=1;}
							else {dam1[i]=0;}
						}
					}
					
					
					ncluster1=0;		
					for(i=0;i<N;i++){
						vis1[i]=0;
					}
					m1=0;
					ncluster1=0;
					for(n=0;n<N;n++){
						if((vis1[n]==0)&&(dam1[n]==1)){
							cluster_size=0;
							ncluster1++;
							cluster_size=Recurrence(n, cluster_size, ncluster1);
							size_cluster1[ncluster1]=cluster_size;
							if(cluster_size>m1){m1=cluster_size;
								c1=ncluster1;}				
						}			
					}
					
					Nc1=c1;
					//m2=0;
					for(i=0;i<N;i++){
						s[i][nrun2]=0;
						if(vis1[i]==Nc1){
							s[i][nrun2]=1;
						}
					}
	
				}
				
				sigma11=0.;
				sigma10=0.;
				for(i=0;i<N;i++){
					sigma11+=(float)(s[i][1]*s[i][0])/((float)N0);
					sigma10+=(float)(s[i][1]*(1-s[i][0]))/((float)N0);
				}
				//printf("%d %lf %lf\n",nc,sigma11,sigma10);
				sigma11m[nc][nc2]+=sigma11/((float)Nrunmax);
				sigma10m[nc][nc2]+=sigma10/((float)Nrunmax);
				
				
			}
		}
	}
	
	
	for(n=0;n<100;n++){
		p=(float)n*0.01;
		p11=pow(p,a);
		if((2*p-p11)<1){
			fprintf(gp2,"%lf %lf %lf %lf %lf\n",p,sigma10m[n][nc2],sigma11m[n][nc2],sigma10m[n][nc2]+sigma11m[n][nc2],p11);
		}
	}
	
	return 0;	
}


	
	
