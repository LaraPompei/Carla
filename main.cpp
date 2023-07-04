#include <iostream>
#include <cmath>
#include <fstream>

#define eq 12
#define t_store 1000 //Intervalo de pontos sendo salvos

//Parametros:

#define pi_v 0.68           //Taxa de replicacao viral
#define c_v1 2.63           //Taxa de clareamento viral maximo pelo sistema inato
#define c_v2 0.6            //Constante de meia saturacao
#define k_v1 0.0000482      //Taxa de neutralizacao do virus por unidade anticorpos neutralizantes
#define k_v2 0.000000748    //Taxa de eliminacao do virus por unidade de celulas T CD8+
#define alpha_ap 0.0025     //Taxa de hosmeostase das APCs imaturas
#define beta_ap 0.55        // Taxa de maturacao das APCs 
#define c_ap1 0.8           //Taxa de maturacao maxima das APCs
#define c_ap2 40            //Constante de meia ativacao
#define delta_apm 0.538     //Taxa de morte das APCs maduras
#define alpha_th 0.000217   //Taxa de gineistase das celulas T CD4+
#define beta_th 0.0000001   //Taxa de replicacao das celulas T CD4+ naive
#define pi_th 0.00000001    //Taxa de replicacao das celulas T CD4+ efetoras
#define delta_th 0.22       //Taxa de morte das celulas T CD4+ efetoras
#define alpha_tk 0.000217   //Taxa de homeostase das celulas T CD8+
#define beta_tk 0.000001    //Taxa de ativacao das celulas T CD8+ naive
#define pi_tk 0.00000001    //Taxa de replicacao das celulas T CD8+ efetoras
#define delta_tk 0.0003     //Taxa de morte das celulas T CD8+ efetoras
#define alpha_b 6.0         //Taxa de homeostase das celulas B
#define pi_b1 0.00000483    //Taxa de ativacao das celulas B T-independente
#define pi_b2 0.0000000127  //Taxa de ativacao das celulas B T-dependentes
#define beta_ps 0.000672    //Taxa de diferenciacao das celulas B ativas em plasmocitos de vida curta
#define beta_pl 0.00000561  //Taxa de diferenciacao das celulas B ativas em plasmocitos de vida longa
#define beta_bm 0.000001    //Taxa de diferenciacao das celulas B ativas em celulas B de memoria
#define delta_ps 2.0        //Taxa de morte dos plasmocitos de vida curta
#define delta_pl 0.00024    //Taxa de morte dos plasmocitos de vida longa
#define gama_bm 0.000975    //Taxa de diferenciacao das celulas B de memoria em plasmocitos de vida longa
#define pi_bm1 0.00001      //Taxa de proliferacao das celulas B de memoria
#define pi_bm2 2500			//Constante de crescimento maximo
#define pi_ps 0.002			//Taxa de secrecao de anticorpos por unidade de plasmocitos de vida curta
#define pi_pl 0.00068       //Taxa de secrecao de anticorpos por unidade de plasmocitos ded vida longa
#define delta_a 0.04        //Taxa de morte de anticorpos

//Condicoes iniciais
#define V0 724
#define Ap0 1000000
#define Apm0 0
#define Thn0 1000000
#define The0 0
#define Tkn0 500000
#define Tke0 0
#define B0 250000
#define Ps0 0
#define Pl0 0
#define Bm0 0
#define A0 150

using namespace std;

void Sistema(double t, double *y, double* dydt){
    //V
    dydt[0] = pi_v*y[0] - (c_v1*y[0])/(c_v2+y[0]) - k_v1*y[0]*y[11] - k_v2*y[0]*y[6];
    //Ap
    dydt[1] = alpha_ap*(Ap0 - y[1]) - beta_ap*y[1]*(c_ap1*(y[0])/(c_ap2 + y[0]));
    //Apm
    dydt[2] = beta_ap*y[1]*(c_ap1*(y[0])/(c_ap2 + y[0])) - delta_apm*y[2];
    //Thn
    dydt[3] = alpha_th*(Thn0 - y[3]) - beta_th*y[2]*y[3];
    //The
    dydt[4] = beta_th*y[2]*y[3] + pi_th*y[2]*y[4] - delta_th*y[4];
    //Tkn
    dydt[5] = alpha_tk*(Tkn0 - y[5]) - beta_tk*y[2]*y[5];
    //Tke
    dydt[6] = beta_tk*y[2]*y[5] + pi_tk*y[2]*y[6] - delta_tk*y[6];
    //B
    dydt[7] = alpha_b*(B0 - y[7]) + pi_b1*y[0]*y[7] + pi_b2*y[4]*y[7] - beta_ps*y[2]*y[7] - beta_pl*y[4]*y[7] - beta_bm*y[4]*y[7];
    //Ps
    dydt[8] = beta_ps*y[2]*y[7] - delta_ps*y[8];
    //Pl
    dydt[9] = beta_pl*y[4]*y[7] - delta_pl*y[9] + gama_bm*y[10];
    //Bm
    dydt[10] = beta_bm*y[4]*y[7] + pi_bm1*y[10]*(1 - y[10]/(pi_bm2)) - gama_bm*y[10];
    //A
    dydt[11] = pi_ps*y[8] + pi_pl*y[9] - delta_a*y[11];
}

void saveData (double** y, double* t, double h, int pont, int inter){
//    cout<<"Saving Data..."<<endl;
    ofstream outputFile("output.csv",ios_base::app);
    for(int i=0; i<pont; i++){
        outputFile<<(inter-(pont*h)+i*h);
        for(int j=0; j<eq;j++){
            outputFile<<","<<y[i][j];
        }
        outputFile<<endl;
    }
    outputFile.close();
}


void RK5(double* t, double h, double** y, int inter){
    double* k1 = new double[eq]; 
    double* k2 = new double[eq]; 
    double* k3 = new double[eq]; 
    double* k4 = new double[eq];
    double* k5 = new double[eq]; 
    double* k6 = new double[eq];
    int i; 
    for(i=1; i<inter; i++){
        
        if(i%t_store==0)
            saveData(y,t,h,t_store,i);
        
        //Calculando K1=f(t[i],y[i])
        if(i%t_store==0){
            Sistema(t[t_store-1],y[t_store-1],k1);
            cout<<"Tempo "<<i*h<<" dias"<<endl;
            cerr<<"Interacao "<<i<<endl;
        }else{
            Sistema(t[i%t_store-1],y[i%t_store-1],k1);
        }

        ofstream outputFile1("k1.csv",ios_base::app);
        for(int k=0; k<eq;k++){
            outputFile1<<"\t"<<k1[k];
        }
        outputFile1<<endl;
        outputFile1.close();
        
        //Calculando K2 = f(t[i]+(1/4)*h,y[i]+(1/4)*k1)
        for(int j=0;j<eq;j++){ 
            y[i%t_store][j] = (i%t_store == 0)?y[t_store-1][j]+(1/4)*k1[j]:y[i%t_store-1][j]+(1/4)*k1[j];
        }
        Sistema(t[i%t_store],y[i%t_store],k2);

        ofstream outputFile2("k2.csv",ios_base::app);
        for(int k=0; k<eq;k++){
            outputFile2<<"\t"<<k2[k];
        }
        outputFile2<<endl;
        outputFile2.close();

        //Calculando K3 = f(t[i]+(3/8)*h,y[i]+(3/32)*k1+(9/32)*k2)
        for(int j=0;j<eq;j++){
            if(i%t_store == 0)
                y[i%t_store][j] = y[t_store-1][j]+(3/32)*k1[j]+(9/32)*k2[j];
            else
                y[i%t_store][j] = y[i%t_store-1][j]+(3/32)*k1[j]+(9/32)*k2[j];
        }
        Sistema(t[i%t_store],y[i%t_store],k3);
        
        ofstream outputFile3("k3.csv",ios_base::app);
        for(int k=0; k<eq;k++){
            outputFile3<<"\t"<<k3[k];
        }
        outputFile3<<endl;
        outputFile3.close();

        //Calculando K4 = f(t[i]+(12/13)*h,y[i]+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3) 
        for(int j=0;j<eq;j++){
            if(i%t_store == 0)
                y[i%t_store][j] = y[t_store-1][j]+(1932/2197)*k1[j]-(7200/2197)*k2[j]+(7296/2197)*k3[j];
            else
                y[i%t_store][j] = y[i%t_store-1][j]+(1932/2197)*k1[j]-(7200/2197)*k2[j]+(7296/2197)*k3[j];
        }
        Sistema(t[i%t_store],y[i%t_store],k4);
        
        ofstream outputFile4("k4.csv",ios_base::app);
        for(int k=0; k<eq;k++){
            outputFile4<<"\t"<<k4[k];
        }
        outputFile4<<endl;
        outputFile4.close();

        //Calculando K5 = f(t[i]+h,y[i]+(439/216)*k1-8*k2+(3680/513)*k3-(845/4104)*k4)
        for(int j=0;j<eq;j++){
            if(i%t_store == 0)
                y[i%t_store][j] = y[t_store-1][j];
            else
                y[i%t_store][j] = y[i%t_store-1][j];
            y[i%t_store][j]+=(439/216)*k1[j]-8*k2[j]+(3680/513)*k3[j]-(845/4104)*k4[j];
        }
        Sistema(t[i%t_store],y[i%t_store],k6);
        
        ofstream outputFile5("k5.csv",ios_base::app);
        for(int k=0; k<eq;k++){
            outputFile5<<"\t"<<k5[k];
        }
        outputFile5<<endl;
        outputFile5.close();
        //Calculando K6 = f(t[i]+0,5*h,y[i]-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5))
        for(int j=0;j<eq;j++){
            y[i%t_store][j]= (i%t_store==0) ? y[t_store-1][j] : y[i%t_store-1][j];
            y[i%t_store][j]+=-(8/27)*k1[j]+2*k2[j]-(3544/2565)*k3[j]+(1859/4104)*k4[j]-(11/40)*k5[j];
        }
        Sistema(t[i%t_store],y[i%t_store],k6);

        ofstream outputFile6("k6.csv",ios_base::app);
        for(int k=0; k<eq;k++){
            outputFile6<<"\t"<<k6[k];
        }
        outputFile6<<endl;
        outputFile6.close();

        //Calculando o y para cada eq h*(16*k1[j]/135 + 6656*k3[j]/12825 + 28561*k4[j]/56430 - 9*k5[j]/50 + 2*k6[j]/55);
        for(int j=0;j<eq;j++){
            y[i%t_store][j] = (i%t_store==0)?y[t_store-1][j]:y[i%t_store-1][j];
            y[i%t_store][j] += h*((37/378)*k1[j]+(250/621)*k3[j]+(125/594)*k4[j]+(512/1771)*k6[j]);
        }
        
		//passo de tempo
		t[int(i%t_store)] = (i!=0) ? t[i%t_store-1] + h: t[i%t_store] + h;

	}
    saveData(y,t,h,i%t_store,i);
   
}

int main(){
	
	double t0 = 0.0;
	double t_final = 45.0;
	double h = 0.01;
	int inter = int(abs(t_final-t0)/h);
    cout<<"Numero de interacoes: "<<inter<<endl;
	double t[t_store] = {t0};
	double** y = new double*[t_store];
   	
	for(int i = 0; i<t_store; i++)
		y[i] = new double[eq];
	
	y[0][0] = V0;
	y[0][1] = Ap0;
   	y[0][2] = Apm0;
	y[0][3] = Thn0;
   	y[0][4] = The0;
	y[0][5] = Tkn0;
	y[0][6] = Tke0;
   	y[0][7] = B0;
	y[0][8] = Ps0;
   	y[0][9] = Pl0; 
	y[0][10] = Bm0;
	y[0][11] = A0;

	ofstream outputFile("output.csv");
	outputFile<<"t,V,Ap,Apm,Thn,The,Tkn,Tke,B,Ps,Pl,Bm,A"<<endl;
    outputFile.close();	
    RK5(t,h,y,inter);

    for (int i = 0; i < t_store; ++i)
		delete [] y[i];
	delete [] y;
	
	return 0;
}