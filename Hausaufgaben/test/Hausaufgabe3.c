#include<stdio.h>
#include<stdlib.h>
#include<math.h>


double s(double x){
	return(0.0);
}

double g(double x, double E, double eps){
	return(2*(E-60*pow(cos(M_PI*x),16)+eps*x));
}


double **numerov(double x0, double x1, double y0, double y1, double xend, double E, double eps){

	double h,faktun1,faktun0,faktunM,faktS;
	int steps;

	h=x1-x0;
	steps=ceil((xend-x0)/h);

	double **y=malloc(sizeof(double*)*steps);
	for(int i=0;i<=steps;i++){
		y[i]=malloc(2*sizeof(double));	//Speicherallozierung für die Tupel
	}

	y[0][0]=x0;  /* belege erste Funktionswerte mit den Startwerten  */
	y[0][1]=y0;
	y[1][0]=x1;
	y[1][1]=y1;

	for(int i=1; i<steps; i++){		//betrachte in den Schritten y(i-1) und y(i) um y(i+1) zu berechnen
		faktun1=(1+(h*h/12.0)*g((i+1)*h,E,eps));
		faktun0=(1-(5*h*h/12.0)*g(i*h,E,eps));
		faktunM=(1+(h*h/12.0)*g((i-1)*h,E,eps));
		faktS=h*h/12.0*(s((i+1)*h)+10.0*s(i*h)+s((i-1)*h));
		//printf("%lf	%lf	%lf	%lf	%lf	%lf\n",faktun1,faktun0,faktunM,faktS,y[i][1],y[i-1][1]);
		
		y[i+1][1]=(faktS+2.0*faktun0*y[i][1]-faktunM*y[i-1][1])/faktun1;
		y[i+1][0]=(i+1)*h;
	}
	return(y);
}

double integrate(double **f, int len, double stepsize) {
	double erg;
	int i=1;

	erg=f[0][1]*f[0][1]+f[len][1]*f[len][1];
	while(i<=len-1){
		erg+=2*f[i][1]*f[i][1];
		i++;;
	}
	return (0.5*erg*stepsize);
}

double fend0(double x,double h, double eps){
	double z,i,**y;
	int j;

	j=ceil(8/h);
	y=numerov(0,h,0,0.1,8,x,eps);
	z=y[j][1];
	free(y);
	return(z);
}

double Nullstelle(double (*f)(double,double,double), double x0, double x1, int N, double h){
	double xn;
	int i=0;

	while(i<=N && !isnan((x1-x0)/(f(x1,h,0)-f(x0,h,0)))){
		xn=x1-f(x1,h,0)*(x1-x0)/(f(x1,h,0)-f(x0,h,0)); 
		x0=x1;
		x1=xn;
		i++;
	}
	
	return(xn);
}

void Ausgabe(double **g, int len,char* savename){
	FILE *fp;

	fp = fopen(savename, "w+");						//Erstellen der Ausgabedatei
	fprintf(fp, "#x-Wert		y-Wert\n");					//Bennenung der Zeilen
	for(int i=0; i<=len; i++){
		fprintf(fp, "%lf	%lf\n", g[i][0], g[i][1]);		//Schreiben der Werte
	}
	fclose(fp);
	free(g);
}

void Eigenwertausgabe(double *erg, int len2, char* name){
	double **g=malloc(sizeof(double*)*len2);		//Speicherallozierung für Ausgabearray
	for(int i=0;i<len2;i++){
		g[i]=malloc(2*sizeof(double));
		g[i][1]=10;								//Einlesen in Ausgabearray
		g[i][0]=erg[i];
	}
	Ausgabe(g,len2-1,name);								//Aufruf der Ausgabefunktion mit Ausgabearray
}

void EnergieEigen(double n,double h, char* name, double eps){
	int zaehler,j;
	double *y=malloc(sizeof(double)*(n+1));

	j=0;
	for(double i=0;i<60;i+=(60/n)){
		y[j]=fend0(i,h,eps);
		j++;
	}
	zaehler=24;
	double *z=malloc(sizeof(double*)*zaehler);
	zaehler=0;
	for(j=0;j<=n;j++){
		if((y[j]>0.0 && y[j+1]<0.0)||(y[j]<0.0 && y[j+1]>0.0)||(y[j]==0.0)){	//Festellung des Nulldurchgangs entweder durch exakte Nullstelle oder durch übergang von Negativ auf Positiv zwischen zwei Werten
			z[zaehler]=j*(60/n);				//Inkrementierung des Nullstellenzählers
			zaehler++;
		}
	}
	Eigenwertausgabe(z,zaehler,name);
	free(y);
	free(z);
}

void EpsilonEnergieEigen(){
	
}

void graph(double E, double h, char* name, double eps){
	double z,**y;
	int len;
	len=ceil(8/h);

	y=numerov(0,h,0,0.1,8,E,eps);
	z=integrate(y,len-1,h);
	for(int i=0;i<len;i++){
		y[i][1]=y[i][1]/sqrt(z);
	}
	Ausgabe(y,len,name);
}

int main(){
	double **y,para;
	char name1[]={"ergebniseps2.tsv"};
	char name2[]={"graph4.tsv"};
	

	//para=Nullstelle(&fend0,0,0.1,100,0.1);
	//y=numerov(0,0.1,0,para,8);
	//Ausgabe(y,80);
	//printf("%lf",para);
	EnergieEigen(60000,0.001,name1,2);
	//graph(39.739,0.001,name2);
}