#include <stdio.h>
#include <math.h>
#define rho 1.293//諸定数　非圧縮なので圧力一定
#define mu0 0.0000183
#define k 0.4
#define A_plus 26
#define K 0.0168
#define C_CP 1.6
#define C_WK 0.25
#define C_KLEB 0.3
double x;//横方向(平板方向)長さ1くらい
double y;//縦方向長さ0.012くらい
unsigned long n_x;//横方向分割数2000000くらい
int n_y;//縦方向分割数1200くらい
unsigned long i_x;//x方向カウンタ
double y_closs;//μ_oとμ_iの入れ替わる点のμ_t
double U;//一様流速20くらい
double dx;//x方向メッシュ大きさ
double dy;//y方向メッシュ大きさ
double mu_t_last[500];//最終的なmu_t分布を入れる。
double y_plus_last[500];
double div_1[500];
double div_2[500];
double div_3[500];
double div2_1[500];
double div2_2[500];
double div2_3[500];

double u_tau_last;
double tau_w;
double Y[2000];
int ranryu(double*,double*,double*,double*);//速度分布を更新する関数
int main(void){
    printf("x range");//入力
    scanf("%lf",&x);
    printf("y range");
    scanf("%lf",&y);
    printf("x mesh");
    scanf("%lu",&n_x);
    printf("y mesh");
    scanf("%d",&n_y);
    printf("velocity");
    scanf("%lf",&U);
    dx=x/(double)n_x;//メッシュ大きさ計算
    Y[0]=0;
    double y_kyoukai=y/3;
    int y_kyoukai_i=(int)(2*n_y/3);
    for(int i=1;i<y_kyoukai_i;i++){
        Y[i]=i*y_kyoukai/y_kyoukai_i;
    }
    for(int i=y_kyoukai_i;i<n_y-1;i++){
        Y[i]=y_kyoukai+(i-y_kyoukai_i)*(y-y_kyoukai)/(n_y-y_kyoukai_i);
    }
    Y[n_y-1]=y;
    
    
    double x_plus;
    double x_minus;
    double det;
    for(int i=1;i<n_y-1;i++){
        x_plus=Y[i+1]-Y[i];
        x_minus=Y[i-1]-Y[i];
        det=x_plus*x_minus*(x_minus-x_plus);
        div2_1[i]=2*x_plus/det;
        div2_2[i]=2*(x_minus-x_plus)/det;
        div2_3[i]=-2*x_minus/det;
        div_1[i]=-x_plus*x_plus/det;
        div_2[i]=(x_plus*x_plus-x_minus*x_minus)/det;
        div_3[i]=x_minus*x_minus/det;
    }
    double u1[n_y];//x方向速度
    double v1[n_y];//y方向速度
    double u2[n_y];
    double v2[n_y];
    u1[0]=0;//平板表面で速度0。滑りなし
    v1[0]=0;
    u2[0]=0;
    v2[0]=0;
    for(int i=1;i<n_y;i++){//初期条件
        u1[i]=U;
        v1[i]=0;
    }
    for(int i=1;i<n_y;i++){//初期条件
        u2[i]=U;
        v2[i]=0;
    }
    FILE *fpu = fopen("u.csv", "w");
	if (fpu == NULL) {
		perror("opening u failed");
		return 1;
	}
    FILE *fpmu = fopen("mu.csv", "w");
	if (fpmu == NULL) {
		perror("opening mu failed");
		return 1;
	}
    
    FILE *fpv = fopen("v.csv", "w");
	if (fpv == NULL) {
		perror("opening v failed");
		return 1;
	}
    FILE *fpu_plus = fopen("u_plus.csv", "w");
	if (fpu_plus == NULL) {
		perror("opening u_plus failed");
		return 1;
	}
    
    FILE *ftau = fopen("tau.csv", "w");
	if (ftau == NULL) {
		perror("opening tau failed");
		return 1;
	}
    for(i_x=0;i_x<(unsigned long)(n_x/2);i_x++){//速度分布を下流に向かって更新。1と2を交互に入れ替える
        ranryu(u1,v1,u2,v2);
        fprintf(ftau,"%lf\n",tau_w);
        ranryu(u2,v2,u1,v1);
        fprintf(ftau,"%lf\n",tau_w);
        if(i_x%1000==0){
            for(int i=0;i<n_y-1;i++){
                fprintf(fpu, "%lf,",u1[i]);
            }
            fprintf(fpu, "%lf\n",u1[n_y-1]);
            for(int i=0;i<n_y-1;i++){
                fprintf(fpv, "%lf,",v1[i]);
            }
            fprintf(fpv, "%lf\n",v1[n_y-1]);
            for(int i=0;i<n_y-1;i++){
                fprintf(fpmu, "%lf,",mu_t_last[i]);
            }
            fprintf(fpmu, "%lf\n",mu_t_last[n_y-1]);
            for(int i=0;i<n_y-1;i++){
                fprintf(fpu_plus, "%lf,",u1[i]/u_tau_last);
            }
            fprintf(fpu_plus,"%lf\n",u1[n_y-1]/u_tau_last);
        }

    }
	
	fclose(fpu);
	fclose(fpmu);
    fclose(fpv);
    fclose(fpu_plus);
    fclose(ftau);
    
    FILE *fyu_plus = fopen("u_y_plus.csv", "w");
	if (fyu_plus == NULL) {
		perror("opening uy_plus failed");
		return 1;
	}
    for(int i=0;i<n_y;i++){
        fprintf(fyu_plus, "%lf,%lf\n",y_plus_last[i],u1[i]/u_tau_last);
    }
    fclose(fyu_plus);
	printf("writing end");
    return 0;
}
int ranryu(double *u_old,double *v_old,double *u_new,double *v_new){
    
    double y_MAX;
    double F_MAX;
    int flag=0;//内層で0、外層で1のフラグ
    double y_plus;
    double F;
    double mu_t_i;
    double mu_t_o;
    double l;
    double F_WAKE;
    double F_KLEB;
    double mu[n_y];//その点の粘性係数
    double mu_w=mu0;//壁付近の粘性係数
    int flag1=0;//1週目かどうかのフラグ
    double tau_w_new=100;//壁付近の剪断力
    double tau_w_old=1;//比較用の古い壁付近の剪断力
    double u_tau;//粘性速度
    double mu_t;
    double mu_0;
    double mu_1;
    double du_dy_0;
    double du_dy_1;
    double dmudu_dydy;
    while(fabs((tau_w_new-tau_w_old)/tau_w_new)>0.0001&flag1<1){//tau_wを求めるためのイタレーション。
        tau_w_old=tau_w_new;
        y_MAX=0;//定数。最大値を入れるので最初は0
        F_MAX=0;
        flag=0;//外層と内層を分けるためのフラグ
        if(flag1==0){//一週目はひとつ左の速度分布で計算
            tau_w_new=u_old[1]/Y[1]*mu0;//壁付近の速度勾配と剪断力をかけてる。壁は速度0なので(u_old[1]-u_old[0])/dy=(u_old[1]-0)/dy
        }
        else{//tau_wを求めるイタレーションの二週目以降はmu_tを加えたu[1.5]付近の粘性係数で計算
            tau_w_new=u_new[1]/Y[1]*mu_w;
        }
        u_tau=sqrt(tau_w_new/rho);
        for(int i_y=1;i_y<n_y-1;i_y++){//Fmaxとymaxを求める。渦度ω=du/dyは一週目は一つ左の速度分布。二週目以降は今いるxの速度分布
            if(flag1==0){
                y_plus=rho*Y[i_y]*u_tau/mu0;
                F=y*fabs(u_old[i_y-1]*div_1[i_y]+u_old[i_y]*div_2[i_y]+div_3[i_y]*u_old[i_y+1])*(1-exp(-y_plus/A_plus));
            }
            else{
                y_plus=rho*Y[i_y]*u_tau/mu_w;
                F=y*fabs(u_new[i_y-1]*div_1[i_y]+u_new[i_y]*div_2[i_y]+div_3[i_y]*u_new[i_y+1])*(1-exp(-y_plus/A_plus));
            }
            if(F>F_MAX){
                F_MAX=F;
                y_MAX=Y[i_y];
            }
        }
        mu[0]=mu0;
        for(int i_y=1;i_y<n_y-1;i_y++){//内層と外層を求めてu,vを更新
            y_plus=rho*Y[i_y]*u_tau/mu0;
            l=k*Y[i_y]*(1-exp(-y_plus/A_plus));
            if(flag1==0){//渦度ω=du/dyは一週目は一つ左の速度分布。二週目以降は今いるxの速度分布
                mu_t_i=rho*pow(l,2)*fabs(u_old[i_y-1]*div_1[i_y]+u_old[i_y]*div_2[i_y]+div_3[i_y]*u_old[i_y+1]);
            }
            else{
                mu_t_i=rho*pow(l,2)*fabs(u_new[i_y-1]*div_1[i_y]+u_new[i_y]*div_2[i_y]+div_3[i_y]*u_new[i_y+1]);
            }
            if(i_y==1){//壁付近の粘性係数の計算。tau_w用
                mu_w=mu_t_i/2+mu0;
            }
            F_WAKE=fmin(y_MAX*F_MAX,C_WK*y_MAX*pow(U,2)/F_MAX);
            F_KLEB=1/(1+5.5*pow((C_KLEB*Y[i_y]/y_MAX),6));
            mu_t_o=K*C_CP*rho*F_WAKE*F_KLEB;//外層のmu_t
            if(flag==0&mu_t_o<mu_t_i){//外層と内層のmu_tの大きさがひっくり返った瞬間にflag1を立てて以下外層で固定
                flag=1;
            }
            if(flag==1){
                mu_t=mu_t_o;
            }
            else{
                mu_t=mu_t_i;
            }
            if(i_x%1000==0){//最後のメッシュのmu_tを取り出す
                mu_t_last[i_y]=mu_t;
            }
            if(i_x%1000==0){//最後のメッシュのy_plusを取り出す
                y_plus_last[i_y]=y_plus;
                u_tau_last=tau_w_new/rho;
            }
            //mu[i_y]=mu0;
            mu[i_y]=mu0+mu_t;
        }
        v_new[0]=0;
        for(int i_y=1;i_y<n_y-1;i_y++){
            du_dy_0=u_old[i_y-1]*div_1[i_y]+u_old[i_y]*div_2[i_y]+div_3[i_y]*u_old[i_y+1];
            du_dy_1=u_old[i_y-1]*div2_1[i_y]+u_old[i_y]*div2_2[i_y]+div2_3[i_y]*u_old[i_y+1];
            mu_0=mu[i_y];
            mu_1=mu[i_y-1]*div_1[i_y]+mu[i_y]*div_2[i_y]+div_3[i_y]*mu[i_y+1];
            dmudu_dydy=mu_0*du_dy_1+mu_1*du_dy_0;
            u_new[i_y] = u_old[i_y] + (- v_old[i_y] * du_dy_0 + dmudu_dydy / rho ) * dx / u_old[i_y]; 
            v_new[i_y]=v_new[i_y-1]+(Y[i_y]-Y[i_y-1])/2.0/dx*(-u_new[i_y-1]+u_old[i_y-1]-u_new[i_y]+u_old[i_y]);
        }
        flag1=flag1+1;//tau_wのイタレーションの回数のカウンタ
    }
    tau_w=tau_w_new;
    return 0;
    
}
