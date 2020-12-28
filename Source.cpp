#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Random.h"

using namespace std;

Random rand_gen;
auto fit_dist = rand_gen.real_unif_dist(-0.3, 0.3);

class Species
{
public:
    int nmin;
    int id;
    bool focal;
    double fit;
    long long int birth;
    Species* father;
    Species()
    {
        nmin = 0;
        birth = 0;
        focal = false;
        fit = fit_dist();
        id = 0;
    }
    Species(int n,  Species* orig, int index)
    {
        id = index;
        nmin = n;
        father = orig;
        if (father)
            focal = father->focal;
        else
            focal = 0;
         fit = fit_dist();
    }
};


bool ContinueOld = false;
double N = 1000;
long long int runs = 1000000;
long long int print_each = 1;
long long int print_next = 1;
double gama0 = 0.3;
double delta = 0.2;
double nu = 0.001;


void ParseArgs(int argc, char* argv[])
{
    for (int i = 1; i < argc; i++)
    {
        char* str_i = argv[i];
        char* pch = strtok(str_i, "=");
        if (strcmp(str_i, "-continue") == 0)
        {
            ContinueOld = true;
        }
        else {
            if (strcmp(pch, "N") == 0)
            {
                pch = strtok(NULL, "=");
                int value = (int) strtod(pch, NULL);
                N = value;
            }
            else if (strcmp(str_i, "gamma") == 0 || strcmp(str_i, "gama") == 0)
            {
                pch = strtok(NULL, "=");
                double value = strtod(pch, NULL);
                gama0 = value;
            }
            else if (strcmp(str_i, "runs") == 0)
            {
                pch = strtok(NULL, "=");
                long long int value = (long long int) strtod(pch, NULL);
                runs = value;
            }
            else if (strcmp(str_i, "print_each") == 0)
            {
                pch = strtok(NULL, "=");
                int value = (int)strtod(pch, NULL);
                print_each = value;
                print_next = print_each;
            }
            if (strcmp(pch, "delta") == 0)
            {
                pch = strtok(NULL, "=");
                double value = strtod(pch, NULL);
                delta = value;
            }
            if (strcmp(pch, "nu") == 0)
            {
                pch = strtok(NULL, "=");
                double value = strtod(pch, NULL);
                nu = value;
            }

        }
    }
};

template<typename T>
void LoadPrevious(string file_name, vector<T> vec)
{
    ifstream File(file_name);
    string str;
    for (int x = 0; x < N; x++)
    {
        getline(File, str);
        char* endptr = NULL;
        vec[x] = strtoll(str.c_str(), &endptr, 10);
    }
    File.close();
};
template<typename T>
void print_vec(string file_name, vector<T> vec) 
{
    ofstream File(file_name);
    for (int x = 0; x < N; x++)
        File << vec[x] << "\n";
    File.close();
};

string ftostr(double x, int precision)
{
    double eps = pow(0.1,precision)/2;
    double y = x - int(x);
    bool z = true;
    if (int(x) > 0)
        z = false;
    string str = to_string(int(x+eps)) + ".";
    for (int i = 0; i < precision; i++)
    {
        y *= 10;
        int digit = int(y + eps);
        str += to_string(digit);
        y -= digit;
        if (digit > 0)
            z = false;
    }
    if (x == 0)
        return (str);
    while (z)
    {
        y *= 10;
        int digit = int(y + eps);
        str += to_string(digit);
        y -= digit;
        if (digit > 0)
            z = false;
    }
    return (str);
};

int main(int argc, char* argv[])
{
    ParseArgs(argc, argv);
    fit_dist = rand_gen.real_unif_dist(-gama0/4, gama0/4);
    auto rand_dist = rand_gen.real_unif_dist(0, 1);
    double var = 0;    double g = pow(gama0, 2) * delta;
    string DIR = ".\\";
    string suffix = "_N" + to_string(int(N)) + "g" + ftostr(g, 5) + "nu" + 
        ftostr(nu, 5) + ".dat";


    vector<Species*> SP;

    vector<long long int> Wins(N, 0);
    vector<long long int> Loses(N, 0);
    vector<long long int> Wins1(N, 0);
    vector<long long int> Loses1(N, 0);
    vector<long long int> TWins(N, 0);
    vector<long long int> TLoses(N, 0);
    vector<long long int> TWins1(N, 0);
    vector<long long int> TLoses1(N, 0);
    long long int *T_temp = new long long int[N]();
    long long int* T_temp1 = new long long int[N]();
    long long int* V_temp = new long long int[N]();
    long long int* T_temp_ptr = T_temp;
    long long int* T_temp1_ptr = T_temp1;
    long long int* V_temp_ptr = V_temp;
    /*   double* mm1 = new double[N]();
    double* mm2 = new double[N]();
    double* visits = new double[N]();*/
    vector<double> mm1(N, 0);
    vector<double> mm2(N, 0);
    vector<double> visits(N, 0);
    vector<double> mm1_mes(N, 0);
    vector<double> mm2_mes(N, 0);
    vector<double> visits_mes(N, 0);
    
    vector<long long int> abun(N, 0);
    double* fit = new double[N];
    double tt = 1 / (N * delta);
    

    if (ContinueOld)
    {
        LoadPrevious(DIR+"Wins"+suffix, Wins);
        LoadPrevious(DIR + "TWins" + suffix, TWins);
        LoadPrevious(DIR + "Loses" + suffix, Loses);
        LoadPrevious(DIR + "TLoses" + suffix, TLoses);
        LoadPrevious(DIR + "mm1" + suffix, mm1);
        LoadPrevious(DIR + "mm2" + suffix, mm2);
        LoadPrevious(DIR + "visits" + suffix, visits);
    }
    
    vector<Species* > grid;
    Species MRCA;
    SP.push_back(&MRCA);
    MRCA.nmin = N - 1;
    for (int i = 1; i < N; i++)
        grid.push_back(&MRCA);
    Species Mutant(1, &MRCA, 1);
    Mutant.focal = 1;
    grid.push_back(&Mutant);
    SP.push_back(&Mutant);
    int sr = 2;
    double x = 1;
    T_temp[0] = 1;
    vector<int> death_list;
    int nred = 1;
    long long int t = 0;
    int n1 = rand_gen.random(N);
    int n2 = rand_gen.random(N);
    int w, l;
    double df = grid[n1]->fit - grid[n2]->fit;
    Species* sw, * sl;
    bool Mutate = false;
    double fit_focal = Mutant.fit * Mutant.nmin;
    double fit_env = MRCA.fit * MRCA.nmin;  
    auto mm1_mes_ptr = mm1_mes.begin();
    auto mm2_mes_ptr = mm2_mes.begin();
    auto visits_mes_ptr = visits_mes.begin();
    auto mm1_ptr = mm1.begin();
    auto mm2_ptr = mm2.begin();
    auto visits_ptr = visits.begin();
    int old_x = x;

    for (long long int run = 0; run < runs; run++)
    {
        for (int step = 0; step < N; step++)
        {
            t++;
            if (rand_dist() < tt)
                //if (rand_gen.rand() < tt)
                {
                fit_focal = 0;
                fit_env = 0;
                for (auto s = SP.begin(); s != SP.end(); s++)
                {
                    (*s)->fit = fit_dist();
                    if ((*s)->focal && (*s)->nmin > 0)
                        fit_focal += (*s)->fit * (*s)->nmin;
                    else if ((*s)->focal == false && (*s)->nmin > 0)
                        fit_env += (*s)->fit * (*s)->nmin;
                }
            }
            n1 = rand_gen.random(N);
            n2 = rand_gen.random(N);
            df = grid[n1]->fit - grid[n2]->fit;
            (*mm1_mes_ptr) += df;
            (*mm2_mes_ptr) += pow(df,2);
            (*visits_mes_ptr)++;

           // if (rand_gen.rand() < 0.5 + df)
            if (rand_dist() < 0.5 + df)
                {
                w = n1;
                l = n2;
            }
            else
            {
                w = n2;
                l = n1;
            }

            sl = grid[l];
            sw = grid[w];
            //if (rand_gen.rand() < nu)
            if (rand_dist() < nu)
                {
                unsigned int new_id = SP.size();
                Species* mut;
                if (death_list.size() > 0)
                {
                    new_id = death_list.back();
                    death_list.pop_back();
                    mut = *(SP.begin() + new_id);
                    mut->nmin = 1;
                    mut->focal = sw->focal;
                    mut->father = sw;
                    mut->birth = t;
                    mut->fit = fit_dist();
                }
                else
                {
                    mut = new Species(1, sw, new_id);
                    mut->birth = t;
                    SP.push_back(mut);
                }
                grid[l] = mut;
                sr = sr + 1;
                if (mut->focal)
                {
                    x++;
                    fit_focal += mut->fit;
                }
                else
                    fit_env += mut->fit;

            }
            else
            {
                grid[l] = grid[w];
                grid[w]->nmin++;
                if (sw->focal)
                {
                    x++;
                    fit_focal += sw->fit;
                }
                else
                    fit_env += sw->fit;
            }
            sl->nmin--;
            if (sl->focal)
            {
                x--;
                fit_focal -= sl->fit;
            }
            else
                fit_env -= sl->fit;

            if (sl->nmin == 0)
            {
                sr = sr - 1;
                if (sl->id < SP.size())
                    death_list.push_back(sl->id);
                else
                {
                    Species* s = SP.back();
                    delete s;
                    SP.pop_back();
                }
            }
            // mutation to preserve color

            if (x > 0 && x < N)
            {
                if (x > old_x)
                {
                    mm1_ptr++;
                    mm2_ptr++;
                    visits_ptr++;
                    mm1_mes_ptr++;
                    mm2_mes_ptr++;
                    visits_mes_ptr++;
                    V_temp_ptr++;
                    T_temp_ptr++;
                    T_temp1_ptr++;
                }
                else if (x < old_x)
                {
                    mm1_ptr--;
                    mm2_ptr--;
                    visits_ptr--;
                    mm1_mes_ptr--;
                    mm2_mes_ptr--;
                    visits_mes_ptr--;
                    V_temp_ptr--;
                    T_temp_ptr--;
                    T_temp1_ptr--;
                }
                old_x = x;

                (*V_temp_ptr)++;
                if (*T_temp1_ptr == 0)
                    *T_temp1_ptr = t;
                 *T_temp_ptr += t;
            }
            else if (x == N)
            {
                auto w = Wins.begin();
                auto w1 = Wins1.begin();
                auto tw = TWins.begin();
                auto tw1 = TWins1.begin();
                long long int* vt = V_temp;
                long long int* t_t = T_temp;
                long long int* t_t1 = T_temp1;
                *w += *vt;
                if ((*t_t) > 0)
                {
                    (*w1)++;
                    *tw1 +=  t - (*t_t1);
                    *tw += (*vt) * t - (*t_t);
                }
                (*vt) = 0;
                (*t_t) = 0;
                (*t_t1) = 0;
                for (int ii = 1; ii < N; ii++)
                {
                    ++vt;
                    ++t_t1;
                    ++t_t;
                    ++w;
                    ++w1;
                    ++tw;
                    ++tw1;

                    if ((*t_t) > 0)
                    {
                        *w += *vt;
                        (*w1) ++;
                        *tw1 += t - (*t_t1);
                        *tw += (*vt) * t - (*t_t);
                    }
                    (*vt) = 0;
                    (*t_t) = 0;
                    (*t_t1) = 0;
                }

                Mutate = true;
            }
            else if (x == 0)
            {
                auto l = Loses.begin();
                auto l1 = Loses1.begin();
                auto tl = TLoses.begin();
                auto tl1 = TLoses1.begin();
                long long int* vt = V_temp;
                long long int* t_t = T_temp;
                long long int* t_t1 = T_temp1;
                *l += *vt;
                if (*t_t > 0)
                {
                    (*l1)++;
                    *tl1 += t - *t_t1;
                    *tl += (*vt)*t - *t_t;
                }
                (*vt) = 0;
                (*t_t) = 0;
                (*t_t1) = 0;
                for (int ii = 1; ii < N; ii++)
                {
                    ++l;
                    ++tl;
                    ++l1;
                    ++tl1;
                    ++vt;
                    ++t_t;
                    ++t_t1;
                    (*l) += *vt;
                    if (*t_t > 0)
                    {
                        (*l1)++;
                        (*tl1) += t - *t_t1;
                        (*tl) += (*vt)*t - *t_t;
                    }
                    (*vt) = 0;
                    (*t_t) = 0;
                    (*t_t1) = 0;
                }

                Mutate = true;
            }

            if (Mutate)
            {
                Mutate = false;
                int n_mut = rand_gen.random(N);
                if (grid[n_mut]->nmin == 1)
                {
                    if (grid[n_mut]->focal)
                        fit_focal -= grid[n_mut]->fit;
                    else
                        fit_env -= grid[n_mut]->fit;
                    grid[n_mut]->focal = !grid[n_mut]->focal;
                    grid[n_mut]->birth = t;
                    grid[n_mut]->father = nullptr;
                    grid[n_mut]->fit = fit_dist();
                    if (grid[n_mut]->focal)
                        fit_focal += grid[n_mut]->fit;
                    else
                        fit_env += grid[n_mut]->fit;
                    if (x == 0)
                        x = 1;
                    else if (x == N)
                        x = N - 1;
                    else
                        cout << "Error - Wrong Forced Mutation (singelton)";

                    }
                else
                {
                    grid[n_mut]->nmin--;
                    if (grid[n_mut]->focal)
                        fit_focal -= grid[n_mut]->fit;
                    else
                        fit_env -= grid[n_mut]->fit;

                    unsigned int new_id = SP.size();
                    if (death_list.size() > 0)
                    {
                        new_id = death_list.back();
                        death_list.pop_back();
                    }
                    Species* mut = new Species(1, nullptr, new_id);
                    if (x == 0)
                    {
                        mut->focal = true;
                        x = 1;
                    }
                    else if (x == N)
                    {
                        mut->focal = false;
                        x = N - 1;
                    }
                    else
                        cout << "Error - Wrong Forced Mutation";
                    if (SP.size() > new_id)
                    {
                        mut->id = new_id;
                        SP[new_id] = mut;
                    }
                    else
                    {
                        mut->id = SP.size();
                        SP.push_back(mut);
                    }
                    grid[n_mut] = mut;
                    if (mut->focal)
                        fit_focal = mut->fit;
                    else
                        fit_env = mut->fit;

                    sr = sr + 1;
                    Mutate = false;
                }
            }

          /*  int sum_n = 0;
            double fit1_test = 0, fit2_test = 0;
            for (unsigned int i = 0; i < SP.size(); i++)
            {
                auto s = SP[i];
                sum_n += s->nmin;
                if (s->focal)
                    fit1_test += s->nmin * s->fit;
                else
                    fit2_test += s->nmin * s->fit;
                if (s->nmin == 0)
                {
                    bool found = false;
                    for (int j = 0; j < death_list.size(); j++)
                        if (death_list[j] == i)
                            found = true;
                    if (!found)
                        cout << "Death list Error";
                }

            }
            if(abs(fit1_test - fit_focal) > 0.0000001)
                cout << "Error: df1 - " << (fit1_test - fit_focal) << "\n";
            if (abs(fit2_test - fit_env) > 0.0000001)
                cout << "Error: df2 - " << (fit2_test - fit_env) << "\n";
            if (sum_n != N)
                cout << "sum_n Error";*/
        }
        

        //if (*mm1_ptr != mm1[x - 1])
        //    cout << "mm1_ptr Error";
        //if (*mm2_ptr != mm2[x - 1])
        //    cout << "mm2_ptr Error";
        //if (*visits_ptr != visits[x - 1])
        //    cout << "visits_ptr Error";
        //if (*V_temp_ptr != V_temp[int(x - 1)])
        //    cout << "V_temp_ptr Error";
        //if (*T_temp_ptr != T_temp[int(x - 1)])
        //    cout << "T_temp_ptr Error";
        //if (*T_temp1_ptr != T_temp1[int(x - 1)])
        //    cout << "T_temp1_ptr Error";
        //if (*mm1_mes_ptr != mm1_mes[x - 1])
        //    cout << "mm1_mes_ptr error";
        //if (*mm2_mes_ptr != mm2_mes[x - 1])
        //    cout << "mm2_mes_ptr error";
        //if (*visits_mes_ptr != visits_mes[x - 1])
        //    cout << "visits_mes_ptr error";

        double deltafit = (fit_focal / x) - (fit_env / ((long double)N - (long double)x));
        (*visits_ptr) ++;
        *mm1_ptr += deltafit;
        *mm2_ptr += pow(deltafit , 2);
 
      

        if(run > print_next)
        {
            cout << "Wins after " << run << " steps: ";
            cout << *Wins1.begin() << "\n";
            print_vec(DIR + "Wins" + suffix, Wins);
            print_vec(DIR + "TWins" + suffix, TWins);
            print_vec(DIR + "Loses" + suffix, Loses);
            print_vec(DIR + "TLoses" + suffix, TLoses);
            print_vec(DIR + "Wins1" + suffix, Wins1);
            print_vec(DIR + "TWins1" + suffix, TWins1);
            print_vec(DIR + "Loses1" + suffix, Loses1);
            print_vec(DIR + "TLoses1" + suffix, TLoses1);
            print_vec(DIR + "mm1" + suffix, mm1);
            print_vec(DIR + "mm2" + suffix, mm2);
            print_vec(DIR + "visits" + suffix, visits);
            print_vec(DIR + "mm1_mes" + suffix, mm1_mes);
            print_vec(DIR + "mm2_mes" + suffix, mm2_mes);
            print_vec(DIR + "visits_mes" + suffix, visits_mes);
            print_next += print_each;
        }
        
}
   



 delete[] fit;
 delete[] T_temp ;
 delete[] T_temp1;
 delete[] V_temp ;
}