#ifndef OTOO_WDEOS_H
#define OTOO_WDEOS_H
#include "EOS.hpp"

extern "C" {
    void helm_init_(void);
    void eosx_(double *, double *, double *, double *, double *, PS::F64 *);
    void eosx_return_(double *, double *, double *, double *, double *, PS::F64 *);
}

namespace OTOO {
    class WDEOS : public EOS {
    public:
        WDEOS(double, double, double, double, float);
        ~WDEOS() {};
        
        virtual float GetS(float d, float e) = 0;
        virtual float GetP(float d, float e) = 0;
        virtual float GetE(float d, float e) = 0;
        virtual float GetT(float d, float e) = 0;
        virtual void TestEOS(uint64 n, float *d, float *e= NULL) = 0;
        virtual float GetEmin(float d) = 0; 
        virtual float GetEmin2(float d) = 0; 

    protected:
#define N_TABLE 3999
        void LoadHelmTable();
        void SetupEOSforDumping();
        void SetupEOS();
        
        uint64 get_index_D(float);
        uint64 get_index_E(float, double *);
        double m(double *, double *);
        double CubicHermite(double *, double *, float);
        double getX(uint64 D_index, float eg, double tab[][N_TABLE]);
        
        double sph_dunit, sph_eunit, sph_vunit, sph_punit;
        float initial_temperature;
        
        double D_table[N_TABLE];
        double T_table[N_TABLE][N_TABLE];
        double E_table[N_TABLE][N_TABLE];
        double Emin_table[N_TABLE];
        double P_table[N_TABLE][N_TABLE];
        double C_table[N_TABLE][N_TABLE];
        double D_min, D_max, d_D;
        double T_min, T_max, d_T;
        static const uint64 E_limit_;
        static const double Emin;
        
    private:
        uint64 bisection(float x, uint64 i_min, uint64 i_max, double *tab);
    };
    const uint64 WDEOS::E_limit_ = 20;
    const double WDEOS::Emin = 1.0e8;
    
    WDEOS::WDEOS(double dunit, double eunit, double vunit, double punit, float tt) 
        : sph_dunit(dunit), 
          sph_eunit(eunit), 
          sph_vunit(vunit), 
          sph_punit(punit), 
          initial_temperature(tt)
    {
        LoadHelmTable();
    }

#if 0 // for distributed memories
    void WDEOS::SetupEOSforDumping() {
        D_min = 1.0e2/sph_dunit;
        D_max = 1.0e10/sph_dunit;
        
        d_D = (log10(D_max)-log10(D_min))/N_TABLE;
        
        FILE *fp = fopen("EOS_table.dat", "w");
        fprintf(fp, "#inittemp %e\n", initial_temperature);
        fprintf(fp, "#D T E P C\n");
        for(uint64 j = 0; j < N_TABLE; j++) {
            double D, E;
            double P, C, T;
            
            D = pow(10.0, log10(D_min)+d_D*j)*sph_dunit;
            T = initial_temperature;
            eosx_return_(&T, &D, &P, &E, &C, &CodeUnit::FractionOfCoulombCorrection);
            
            D_table[j]    = D/sph_dunit;
            E_table[0][j] = E/sph_eunit;
            T_table[0][j] = T;
            P_table[0][j] = P/sph_punit;
            C_table[0][j] = C/sph_vunit;
            
            fprintf(fp, "%e %e %e %e %e\n", 
                    D_table[j], 
                    T_table[0][j], 
                    E_table[0][j], 
                    P_table[0][j], 
                    C_table[0][j]);
        }
        fclose(fp);
    }
#else
    void WDEOS::SetupEOSforDumping() {
        D_min = 1.0e2/sph_dunit;
        D_max = 1.0e10/sph_dunit;
        
        d_D = (log10(D_max)-log10(D_min))/N_TABLE;

        if(PS::Comm::getRank() == 0) {
            FILE *fp = fopen("EOS_table.dat", "w");
            fprintf(fp, "#inittemp %e\n", initial_temperature);
            fprintf(fp, "#D T E P C\n");
            for(uint64 j = 0; j < N_TABLE; j++) {
                double D, E;
                double P, C, T;
                
                D = pow(10.0, log10(D_min)+d_D*j)*sph_dunit;
                T = initial_temperature;
                eosx_return_(&T, &D, &P, &E, &C, &CodeUnit::FractionOfCoulombCorrection);
                
                D_table[j]    = D/sph_dunit;
                E_table[0][j] = E/sph_eunit;
                T_table[0][j] = T;
                P_table[0][j] = P/sph_punit;
                C_table[0][j] = C/sph_vunit;
                
                fprintf(fp, "%e %e %e %e %e\n", 
                        D_table[j], 
                        T_table[0][j], 
                        E_table[0][j], 
                        P_table[0][j], 
                        C_table[0][j]);
            }
            fclose(fp);
        } else {
            for(uint64 j = 0; j < N_TABLE; j++) {
                double D, E;
                double P, C, T;
                
                D = pow(10.0, log10(D_min)+d_D*j)*sph_dunit;
                T = initial_temperature;
                eosx_return_(&T, &D, &P, &E, &C, &CodeUnit::FractionOfCoulombCorrection);
                
                D_table[j]    = D/sph_dunit;
                E_table[0][j] = E/sph_eunit;
                T_table[0][j] = T;
                P_table[0][j] = P/sph_punit;
                C_table[0][j] = C/sph_vunit;                
            }
        }

        PS::Comm::barrier();
    }
#endif
    
    void WDEOS::SetupEOS()
    {
        D_min = 1.0e-8;
        D_max = 1.0e10;
        D_min /= sph_dunit;
        D_max /= sph_dunit;
        
        T_min = 1.5e3;
        T_max = 1.0e12;
        
        d_D = (log10(D_max)-log10(D_min))/N_TABLE;
        d_T = (log10(T_max)-log10(T_min))/N_TABLE;
        
        std::stringstream eos_file;
        eos_file << "EOS_" << N_TABLE << ".dat";
        
        int flag = 0;
        int ntab;
        FILE *fp2 = fopen(eos_file.str().c_str(), "r");
        if (fp2 != NULL) {
            fread(&ntab, sizeof(int), 1, fp2);
            double dunit_check;
            fread(&dunit_check, sizeof(double), 1, fp2);
            if (ntab == N_TABLE && dunit_check == sph_dunit) {
                fread(D_table, sizeof(double), ntab, fp2);
                fread(E_table, sizeof(double), ntab*ntab, fp2);
                fread(T_table, sizeof(double), ntab*ntab, fp2);
                fread(P_table, sizeof(double), ntab*ntab, fp2);
                fread(C_table, sizeof(double), ntab*ntab, fp2);
                flag = 1;	
            } 
            fclose(fp2);
        }

#if 0 // for distributed memories
        if (flag == 0) {
            std::cerr << "Generating EOS with N = " << N_TABLE << "\n";
            double D, E, P, C, T;
            for(uint64 j = 0; j < N_TABLE; j++) {
                D = pow(10.0, log10(D_min)+d_D*j);
                D_table[j] = D;
                D *= sph_dunit;
                for(uint64 i = 0; i < N_TABLE; i++) {
                    T = pow(10.0, log10(T_min)+d_T*i);
                    
                    eosx_return_(&T, &D, &P, &E, &C, &CodeUnit::FractionOfCoulombCorrection);
                    
                    E_table[j][i] = E/sph_eunit;
                    T_table[j][i] = T;
                    P_table[j][i] = P/sph_punit;
                    C_table[j][i] = C/sph_vunit;
                }
                if (j % 25 == 0) fputs(".", stderr); fflush(stderr);
            }
            std::cerr << "\n\n";
            std::cerr.flush();
            
            int n = N_TABLE;
            fp2 = fopen(eos_file.str().c_str(), "w");
            fwrite(&n, sizeof(int), 1, fp2);
            fwrite(&sph_dunit, sizeof(double), 1, fp2);
            fwrite(D_table, sizeof(double), n, fp2);
            fwrite(E_table, sizeof(double), n*n, fp2);
            fwrite(T_table, sizeof(double), n*n, fp2);
            fwrite(P_table, sizeof(double), n*n, fp2);
            fwrite(C_table, sizeof(double), n*n, fp2);
            fclose(fp2);
        }
#else
        if (flag == 0) {
            if(PS::Comm::getRank() == 0) {
                std::cerr << "Generating EOS with N = " << N_TABLE << "\n";
            }
            double D, E, P, C, T;
            for(uint64 j = 0; j < N_TABLE; j++) {
                D = pow(10.0, log10(D_min)+d_D*j);
                D_table[j] = D;
                D *= sph_dunit;
                for(uint64 i = 0; i < N_TABLE; i++) {
                    T = pow(10.0, log10(T_min)+d_T*i);
                    
                    eosx_return_(&T, &D, &P, &E, &C, &CodeUnit::FractionOfCoulombCorrection);
                    
                    E_table[j][i] = E/sph_eunit;
                    T_table[j][i] = T;
                    P_table[j][i] = P/sph_punit;
                    C_table[j][i] = C/sph_vunit;
                }
                if (j % 25 == 0) fputs(".", stderr); fflush(stderr);
            }
            if(PS::Comm::getRank() == 0) {
                std::cerr << "\n\n";
                std::cerr.flush();
                
                int n = N_TABLE;
                fp2 = fopen(eos_file.str().c_str(), "w");
                fwrite(&n, sizeof(int), 1, fp2);
                fwrite(&sph_dunit, sizeof(double), 1, fp2);
                fwrite(D_table, sizeof(double), n, fp2);
                fwrite(E_table, sizeof(double), n*n, fp2);
                fwrite(T_table, sizeof(double), n*n, fp2);
                fwrite(P_table, sizeof(double), n*n, fp2);
                fwrite(C_table, sizeof(double), n*n, fp2);
                fclose(fp2);
            }
            PS::Comm::barrier();
        }
#endif
        
        for(uint64 j = 0; j < N_TABLE-2; j++) {
            if (D_table[j]*sph_dunit < 100.0) {
                Emin_table[j] = E_table[j][5];
            } else {
                double dum = E_table[j][E_limit_];
                //	for(uint64 i = -1; i <=2 ; i++) {
                for(int i = -1; i <=2 ; i++) {
                    dum = std::max(dum, E_table[std::min(j+i, (uint64)(N_TABLE-1))][E_limit_]);
                }
                Emin_table[j] = dum;
            }
            Emin_table[j] = std::min(Emin_table[j], Emin);
        }
        
        /*
          FILE *fp0 = fopen("EOS.txt", "w");
          for(int j = 0; j < N_TABLE; j++) {
          for(int i = 0; i < N_TABLE; i++) {
          fprintf(fp0, "%e %e %e %e %e\n", 
          D_table[j]*sph_dunit, 
          E_table[j][i]*sph_eunit, 
          T_table[j][i], 
          P_table[j][i]*sph_punit, 
          C_table[j][i]*sph_vunit);
          }
          }
          fclose(fp0);
        */
    }
    
    void WDEOS::LoadHelmTable()
    {
        double dd, ee, temp0, pres0, cs0;
        helm_init_();
        
        dd = -5.0;
        ee = 12.0;
        
        dd = pow(10.0, dd);
        ee = pow(10.0, ee);
        
        temp0 = 5.0e6;
        eosx_(&temp0, &dd, &pres0, &ee, &cs0, &CodeUnit::FractionOfCoulombCorrection);
    }
    
    uint64 WDEOS::bisection(float x, uint64 i_min, uint64 i_max, double *tab)
    {
        uint64 i_mid;
        double x_mid;
        do {
          i_mid = (i_min+i_max)/2;
          x_mid = tab[i_mid];
          if (x > x_mid) i_min = i_mid;
          else i_max = i_mid;
        } while(i_max - i_min != 1);
        return i_min;
    }
    
    uint64 WDEOS::get_index_D(float x)
    {
//        uint64 i_min, i_max;
//        i_min = 2;
//        i_max = N_TABLE-2;
        
        return (uint64)((log10(x)-log10(D_min))/d_D);
        
        /*
          uint64 ib = bisection(x, i_min, i_max, D_table);
        */
    }
    
    uint64 WDEOS::get_index_E(float x, double *tab)
    {
        uint64 i_min, i_max;
//        double x_min, x_max;
        double x_min;
        
        i_min = 2;
        i_max = N_TABLE-2;
        
        x_min = tab[i_min];
//        x_max = tab[i_max];
        
        /*
          if (x_max < x) {
          std::cerr << "get_index_E: out of range !\n";
          std::cerr << x << "\t" << tab[i_max] << "\n";
          exit(-1);
          }
        */
        
        if (x < x_min) {
      return E_limit_;
        }
        
        return bisection(x, i_min, i_max, tab);
    }
    
    double WDEOS::m(double *x, double *y)
    {
        return 0.5*((y[2]-y[1])/(x[2]-x[1])+(y[1]-y[0])/(x[1]-x[0]));
    }
    
    double WDEOS::CubicHermite(double *x, double *y, float q)
    {
        // 0 1   q   2 3 
        double p0 = y[1];
        double p1 = y[2]; 
        double h00, h10, h01, h11;
        double diff = x[2]-x[1];

        if (diff == 0.0) {
            return 0.5*(y[1]+y[2]);
        } 
        
        
        double t = (x[2]-(double)q)/diff;
        if (t > 1.0) {
            return 0.5*(y[1]+y[2]);
        }
        
        h00 = (1.0+2.0*t)*(1.0-t)*(1.0-t);
        h10 = t*(1.0-t)*(1.0-t);
        h01 = t*t*(3.0-2.0*t);
        h11 = t*t*(t - 1.0);
        
        return h00*p0 + h10*diff*m(&x[0], &y[0]) + h01*p1 + h11*diff*m(&x[1], &y[1]);
    }	
    
    double WDEOS::getX(uint64 D_index, float eg, double tab[][N_TABLE])
    {
        uint64 E_index;
        E_index = get_index_E(eg, E_table[D_index]);
        return CubicHermite(&E_table[D_index][E_index-1], &tab[D_index][E_index-1], eg);
    }
    
    class WDEOSforDumping : public WDEOS {
    public:
        WDEOSforDumping(double dunit, double eunit, double vunit, double punit, float tt) 
            : WDEOS(dunit, eunit, vunit, punit, tt)
            {
                SetupEOSforDumping();
            }
        
        ~WDEOSforDumping() {};
        
        float GetEmin(float d)
            {
                return E_table[0][0];
            }
        
        float GetEmin2(float d)
            {
                return E_table[0][0];
            }
        
        float GetS(float d, float e = 0)
            {
                uint64 D_index = get_index_D(d);
                return CubicHermite(&D_table[D_index-1], &C_table[0][D_index-1], d);
            }
        
        float GetP(float d, float e = 0)
            {
                uint64 D_index = get_index_D(d);
                return CubicHermite(&D_table[D_index-1], &P_table[0][D_index-1], d);
            }
        
        float GetE(float d, float e = 0)
            {
                uint64 D_index = get_index_D(d);
                return CubicHermite(&D_table[D_index-1], &E_table[0][D_index-1], d);
            }
        
        float GetT(float d, float e = 0)
            {
                return initial_temperature;
            }
        
        void TestEOS(uint64 n, float *d, float *e = NULL)
            {
                double p_sum, c_sum, e_sum;
                p_sum = c_sum = e_sum = 0.0;
                std::cerr << "EOS TABLE " << N_TABLE << "\n";
                
                for(uint64 i = 0; i < n; i++) {
                    float eg   = GetE(d[i]);
                    float pres = GetP(d[i]);
                    float cs   = GetS(d[i]);
                    
                    double rho1, temp1, pres1, eg1, cs1;
                    rho1 = d[i]*sph_dunit;
                    temp1 = initial_temperature;
                    eosx_return_(&temp1, &rho1, &pres1, &eg1, &cs1, &CodeUnit::FractionOfCoulombCorrection);
                    
                    eg1   /= sph_eunit;
                    pres1 /= sph_punit;
                    cs1   /= sph_vunit;
                    
                    e_sum += fabs((eg-eg1)/eg1);
                    p_sum += fabs((pres-pres1)/pres1);
                    c_sum += fabs((cs-cs1)/cs1);
                }
                std::cerr << "EOS error E " << e_sum/n << "\n";
                std::cerr << "EOS error P " << p_sum/n << "\n";
                std::cerr << "EOS error C " << c_sum/n << "\n";
            }
    };
    
    class WDEOS_D_E : public WDEOS {
    public:
        WDEOS_D_E(double dunit, double eunit, double vunit, double punit, float tt) 
            : WDEOS(dunit, eunit, vunit, punit, tt)
            {
                SetupEOS();
            }
        
        ~WDEOS_D_E() {};
        
        float GetEmin(float d) {
            uint64 D_index = get_index_D(d);
            return Emin_table[D_index];
        }

        float GetEmin2(float d) {
            uint64 D_index = get_index_D(d);
            return E_table[D_index][E_limit_];
        }
        
        float GetS(float d, float e) {
            float tmp = 0.0;
            tmp = (float)getX_D_E(d, e, C_table);
            if (tmp < 0.0) {
                std::cerr << "\n\n";
                std::cerr << d*sph_dunit << " " << e*sph_eunit << "\n";
                std::cerr << "negative S\n";
                {
                    double X_yy[4];
                    uint64 D_index = get_index_D(d);
                    X_yy[0] = getX(D_index-1, e, C_table);
                    X_yy[1] = getX(D_index-0, e, C_table);
                    X_yy[2] = getX(D_index+1, e, C_table);
                    X_yy[3] = getX(D_index+2, e, C_table);
                    
                    uint64 E_index = get_index_E(e, E_table[D_index]);
                    
                    std::cout << d*sph_dunit << " " << D_index << "\n";
                    std::cout << e*sph_eunit << " " << E_index << "\n";
                    
                    //	  return CubicHermite(&D_table[D_index-1], X_yy, d);
                    std::cout << "tab " << X_yy[0] << "\n";
                    std::cout << "tab " << X_yy[1] << "\n";
                    std::cout << "tab " << X_yy[2] << "\n";
                    std::cout << "tab " << X_yy[3] << "\n";
                }
            }
            return tmp;
        }
        
        float GetP(float d, float e) {
            float tmp = 0.0;
            tmp = (float)getX_D_E(d, e, P_table);
            if (tmp < 0.0) {
                std::cerr << "\n\n";
                std::cerr << d*sph_dunit << " " << e*sph_eunit << "\n";
                std::cerr << "negative P\n";
            }
            return tmp;
        }
        
        float GetE(float d, float e) { 
            std::cerr << "GetE:: Don't call me!!!!!!\n";
            return -1.0;
        }
        
        float GetT(float d, float e) {
            float tmp = 0.0;
            uint64 D_index = get_index_D(d);
            if (e < E_table[D_index][E_limit_]) {
                tmp = T_table[D_index][E_limit_];
            } else if (e > E_table[D_index][N_TABLE-1]) {
                tmp = T_table[D_index][N_TABLE-1];
            } else {
                tmp = (float)getX_D_E(d, e, T_table);
                if (tmp > 1.0e12) {
                    std::cerr << "\n\n";
                    std::cerr << d*sph_dunit << " " << e*sph_eunit << "\n";
                    std::cerr << "too high temp\n";
                    exit(-1);
                }
                if (tmp < 0.0) {
                    uint64 E_index = get_index_E(e, E_table[D_index]);
                    E_index = std::max(E_index, (uint64)0);
                    if (T_table[D_index][E_index] < 0.0) {
                        std::cerr << "\n\n";
                        std::cerr << d*sph_dunit << " " << e*sph_eunit << "\n";
                        std::cerr << "negative temp\n";
                        uint64 D_index = get_index_D(d);
                        uint64 E_index = get_index_E(e, E_table[D_index]);
                        
                        std::cout << d*sph_dunit << " " << D_index << "\n";
                        std::cout << e*sph_eunit << " " << E_index << "\n";
                        
                        std::cout << tmp << "\n";
                        std::cout << "tab " << T_table[D_index-1][E_index] << "\n";
                        std::cout << "tab " << T_table[D_index][E_index] << "\n";
                        std::cout << "tab " << T_table[D_index+1][E_index] << "\n";
                        
                        std::cout << "tab " << T_table[D_index][E_index-1] << "\n";
                        std::cout << "tab " << T_table[D_index][E_index] << "\n";
                        std::cout << "tab " << T_table[D_index][E_index+1] << "\n";
                        exit(-1);
                    } else {
                        return T_table[D_index][E_index];
                    }
                }
            }
            return tmp;
        }
        
        void TestEOS(uint64 n, float *d, float *e);
        
    private:
        double getX_D_E(const float, const float, double tab[][N_TABLE]);
    };
    
    double WDEOS_D_E::getX_D_E(const float d, const float e, double tab[][N_TABLE]) {
        double X_yy[4];
        uint64 D_index = get_index_D(d);
        X_yy[0] = getX(D_index-1, e, tab);
        X_yy[1] = getX(D_index-0, e, tab);
        X_yy[2] = getX(D_index+1, e, tab);
        X_yy[3] = getX(D_index+2, e, tab);
        return CubicHermite(&D_table[D_index-1], X_yy, d);
    }
    
    void WDEOS_D_E::TestEOS(uint64 n, float *d, float *e) {
        std::cerr << "EOS TABLE " << N_TABLE << "\n";
        /*
          {
          double rho1, temp1, pres1, eg1, cs1;
          rho1 = 8.428227e+06;
          eg1  = 1.423587e+17;
          
          temp1 = 1.0e6;
          eosx_(&temp1, &rho1, &pres1, &eg1, &cs1);
          
          rho1 /= sph_dunit;
          eg1  /= sph_eunit;
          
          float pres = GetP(rho1, eg1);
          float cs   = GetS(rho1, eg1);
          float tmp  = GetT(rho1, eg1);
          
          std::cout << pres*sph_punit << "\t" << pres1 << "\n";
          std::cout << cs*sph_vunit << "\t" << cs1 << "\n";
          std::cout << tmp << "\t" << temp1 << "\n";
          }
        */
        double p_sum, c_sum, t_sum; 
        p_sum = c_sum = t_sum = 0.0;
        for(uint64 i = 0; i < n; i++) {
            /*
              std::cout << d[i] << "\n";
              std::cout << e[i] << "\n";
              int D_index = get_index_D(d[i], D_table);
              int E_index = get_index_E(e[i], E_table[D_index]);
              std::cout << d[i]*sph_dunit << " " << D_index << "\n";
              std::cout << e[i]*sph_eunit << " " << E_index << "\n";
              std::cout << E_table[D_index-1][0] << "\n";
              std::cout << E_table[D_index-0][0] << "\n";
              std::cout << E_table[D_index+1][0] << "\n";
              std::cout << E_table[D_index+2][0] << "\n";
            */
            
            double rho1, temp1, pres1, eg1, cs1;
            rho1 = d[i]*sph_dunit;
            eg1  = e[i]*sph_eunit;
            temp1 = 1.0e6;
            eosx_(&temp1, &rho1, &pres1, &eg1, &cs1, &CodeUnit::FractionOfCoulombCorrection);
            
            pres1 /= sph_punit;
            cs1   /= sph_vunit;
            /*
              std::cout << P_table[D_index][E_index]*sph_punit 
              << " " << pres1*sph_punit  << " " 
              << P_table[D_index+1][E_index]*sph_punit  
              << "\n";
            */
            float pres = GetP(d[i], e[i]);
            float cs   = GetS(d[i], e[i]);
            float tmp  = GetT(d[i], e[i]);
            
            //      if (i % 14000 == 0) {
            //	std::cout << rho1 << " " << eg1 << " " << temp1 << " " << tmp << " " << pres1 << " " << pres << "\n";
            //      }
            
            /*
              std::cout << pres*sph_punit << " " 
              << pres1*sph_punit  << " " 
              << cs*sph_punit << " " 
              << cs1*sph_punit << " " 
              << "\n";
            */
            p_sum += fabs((pres-pres1)/pres1);
            c_sum += fabs((cs-cs1)/cs1);
            t_sum += fabs((tmp-temp1)/temp1);
        }
        std::cerr << "EOS error P " << p_sum/n << "\n";
        std::cerr << "EOS error C " << c_sum/n << "\n";
        std::cerr << "EOS error T " << t_sum/n << "\n";
        return;
    }
}
#endif
