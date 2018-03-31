#pragma once
#include <cmath>
#include <iomanip> // std::setprecision
#include <iostream>
#include <vector>

#define MIN -32768
#define MAX 32768
#define TRIES 100
#define TAIL 0
#define SHUFFLE 256
#define MULT 25173
#define INCR 13849
#define MOD ((long int)65536)

#define print_out false
#define print_evolution false

#define flip() ((int)((newrand() * (long)2) / (long)65535))

using namespace std;

class Genocop
{
public:
    typedef vector<vector<float> > t_matrix;
    typedef vector<float> t_vector;
    typedef vector<vector<int> > t_imatrix;
    typedef vector<int> t_ivector;

    // random
    long rseed;
    unsigned int rand_array[SHUFFLE];
    // fin random

    int total_variables;
    int total_eq;
    int total_ineq;
    int total_domains;

    int pop_size;
    int generations;
    int P1, P2, P3, P4, P5, P6, P7;
    float Q;
    int isMinimization;
    int init_val;
    int B;
    int STEP;
    int test_num;

    t_matrix final_mat;
    t_matrix equalities;
    t_ivector eq_co;
    t_vector eq_rhs;
    t_matrix a1;
    t_matrix a2;
    t_matrix inv_a1;
    t_matrix inva1_a2;
    t_vector inva1_b;
    t_matrix new_in_eq;

    t_matrix inequalities;
    t_ivector ineq_co;
    t_vector ineq_rhs;
    t_matrix c1;
    t_matrix c2;
    t_matrix org_ineq;

    t_matrix domains;
    t_vector ldomain;
    t_vector udomain;

    t_vector l1;
    t_vector u1;
    t_ivector x1;

    t_vector l2;
    t_vector u2;
    t_ivector x2;

    t_vector X;
    t_imatrix var_order;
    t_ivector cart;

    float Ackley(t_vector& X, float _d)
    {
        float _a, _b, _c, sum, sum2;
        _a = 20.0;
        _b = 0.2;
        _c = 2 * 3.1415926536f;

        sum = 0.0;
        sum2 = 0.0;
        for(int i = 1; i <= (int)_d; i++) {
            sum += X[i] * X[i];
            sum2 += cos(_c * X[i]);
        }

        return (-_a * exp(-_b * sqrt((1.0 / _d) * sum)) - exp((1.0f / _d) * sum2) + _a + exp(1));
    }

    float Schwefel(t_vector& X, int _d)
    {
        float sum = 0.0f;
        for(int i = 1; i <= (int)_d; i++) {
            sum += X[i] * sin(sqrt(abs(X[i])));
        }

        return (418.9829f * _d - sum);
    }

    float function3(t_vector& X, int _d)
    {
        float sum = 0.0;
        for(int i = 1; i <= (int)_d; i++) {
            sum += X[i] * X[i];
        }

        return (0.5f - (pow(sin(sqrt(sum)), 2) - 0.5f) / (pow(1.0f + 0.001f * sum, 2)));
    }

    float evaluate(t_vector& X)
    {

        float _a, _b, _c, _d, sum, sum2;

        switch(test_num) {

        // 2 dimensiones
        case 1:
            return Ackley(X, 2.0f);
        case 2:
            return Schwefel(X, 2.0f);
        case 3:
            return function3(X, 2.0f);
        // 8 dimensiones
        case 4:
            return Ackley(X, 8.0f);
        case 5:
            return Schwefel(X, 8.0f);
        case 6:
            return function3(X, 8.0f);

        default:
            if(print_out)
                cout << "Invalid test case " << test_num << endl;
            return 0.0;
        }

        return 0.0f;
    }

    void resize_matrix(int rows, int cols, t_matrix& matrix)
    {
        matrix.resize(rows + 1);
        for(int i = 0; i < rows + 1; i++) {
            matrix[i].resize(cols + 1);
        }
    }
    void resize_matrix(int rows, int cols, t_imatrix& matrix)
    {
        matrix.resize(rows + 1);
        for(int i = 0; i < rows + 1; i++) {
            matrix[i].resize(cols + 1);
        }
    }

    void resize_vector(int size, t_vector& _vector)
    {
        _vector.resize(size + 1);
    }

    void resize_vector(int size, t_ivector& _vector)
    {
        _vector.resize(size + 1);
    }

    void create_from_file(const string& path)
    {
        float t1, t3;
        int t2;
        FILE* input;
        input = fopen(path.c_str(), "r");
        if(input == NULL) {
            if(print_out)
                cout << "No se puede leer el archivo" << endl;
            return;
        }

        seed();

        fscanf(input, " %d", &total_variables);
        fscanf(input, " %d", &total_eq);
        fscanf(input, " %d", &total_ineq);
        fscanf(input, " %d", &total_domains);

        int fin_r = total_ineq + total_variables;
        int fin_c = total_variables - total_eq + 2;
        int org_col = total_variables - total_eq + 1;
        int x2_vari = total_variables - total_eq;
        int newin_r = total_eq;
        int newin_c = fin_c;
        int a1a2_r = total_eq;
        int a1a2_c = org_col;

        // reservamos memoria para los vectores y matrices
        resize_matrix(fin_r, fin_c, final_mat);
        resize_matrix(total_eq, total_variables + 1, equalities);
        resize_vector(total_variables, eq_co);
        resize_vector(total_eq, eq_rhs);
        resize_matrix(total_eq, total_eq, a1);
        resize_matrix(total_eq, x2_vari, a2);
        resize_matrix(total_eq, total_eq, inv_a1);
        resize_matrix(total_eq, x2_vari, inva1_a2);
        resize_vector(total_eq, inva1_b);
        resize_matrix(total_eq, fin_c, new_in_eq);

        resize_matrix(total_ineq, total_variables + 1, inequalities);

        resize_vector(total_variables, ineq_co);
        resize_vector(total_ineq, ineq_rhs);

        resize_matrix(total_ineq, total_eq, c1);
        resize_matrix(total_ineq, x2_vari, c2);
        resize_matrix(total_ineq, org_col, org_ineq);
        resize_matrix(total_variables, 3, domains);

        ldomain.resize(total_variables);
        udomain.resize(total_variables);

        resize_vector(total_eq, l1);
        resize_vector(total_eq, u1);
        resize_vector(total_eq, x1);

        resize_vector(x2_vari, l2);
        resize_vector(x2_vari, u2);
        resize_vector(x2_vari, x2);

        resize_vector(total_variables, X);
        resize_matrix(total_variables, 2, var_order);
        resize_vector(total_eq, cart);

        // leemos y llenamos las ecuaciones
        if(total_eq != 0) {
            for(int i = 1; i <= total_eq; i++) {
                for(int j = 1; j <= total_variables + 1; j++) {
                    fscanf(input, " %f", &equalities[i][j]);
                }
            }
        }
        // leemos y llenamos las inecuaciones
        if(total_ineq != 0) {
            for(int i = 1; i <= total_ineq; i++) {
                for(int j = 1; j <= total_variables + 1; j++) {
                    fscanf(input, " %f", &inequalities[i][j]);
                }
            }
        }
        // leemos y llenamos las restricciones de dominio
        for(int i = 1; i <= total_variables; i++) {
            domains[i][1] = MIN;
            domains[i][2] = (float)i;
            domains[i][3] = MAX;
        }
        if(total_domains != 0) {
            for(int i = 1; i <= total_domains; i++) {

                fscanf(input, " %f %d %f", &t1, &t2, &t3);
                domains[t2][1] = t1;
                domains[t2][3] = t3;
            }
        }

        // inicializacion
        for(int i = 1; i <= total_variables; i++) {
            eq_co[i] = i;
            ineq_co[i] = i;
        }
        for(int i = 1; i <= total_eq; i++) {
            eq_rhs[i] = equalities[i][total_variables + 1];
        }
        for(int i = 1; i <= total_ineq; i++) {
            ineq_rhs[i] = inequalities[i][total_variables + 1];
        }

        if(total_eq != 0) {

            get_var_order(total_variables, total_eq, cart, var_order);
            find_x1_x2(total_variables, var_order, x1, x2);
            find_ac1_ac2(total_eq, total_eq, x2_vari, x1, x2, equalities, a1, a2);
            inverse(a1, inv_a1, total_eq);
            mmprod(total_eq, total_eq, x2_vari, inva1_a2, inv_a1, a2);
            mvprod(total_eq, total_eq, inva1_b, inv_a1, eq_rhs);
            find_ac1_ac2(total_eq, total_ineq, x2_vari, x1, x2, inequalities, c1, c2);
            find_limits(total_variables, domains, ldomain, udomain);
            find_lu1_lu2(total_variables, total_eq, x1, x2, ldomain, l1, l2);
            find_lu1_lu2(total_variables, total_eq, x1, x2, udomain, u1, u2);
            find_new_in_eq(inva1_b, inva1_a2, l1, u1, newin_r, newin_c, new_in_eq);
            find_org_in_eq(inva1_b, inva1_a2, ineq_rhs, c1, c2, total_ineq, a1a2_r, a1a2_c, org_ineq);

            initialize(final_mat, fin_r, fin_c);

            find_final_mat1(l2, u2, final_mat, x2_vari, fin_c);
            find_final_mat2(new_in_eq, total_eq, fin_c, org_col, final_mat);
            find_final_mat3(org_ineq, total_ineq, org_col, total_variables + 1, final_mat);

        } else {
            for(int i = 1; i <= total_variables; i++) {
                l2[i] = domains[i][1];
                x2[i] = domains[i][2];
                u2[i] = domains[i][3];
            }

            initialize(final_mat, fin_r, fin_c);
            find_final_mat1(l2, u2, final_mat, total_variables, fin_c);

            if(total_ineq != 0) {
                find_final_mat3(inequalities, total_ineq, org_col, total_variables + 1, final_mat);
            }
        }

        // terminamos de leer los demás parametros
        fscanf(input, "%d %d %d %d %d %d %d %d %d", &pop_size, &generations, &P1, &P2, &P3, &P4, &P5, &P6, &P7);
        fscanf(input, "%f ", &Q);
        fscanf(input, "%d ", &isMinimization);
        fscanf(input, "%d ", &init_val);
        fscanf(input, "%d ", &B);
        fscanf(input, "%d ", &STEP);
        fscanf(input, "%d ", &test_num);
        fclose(input);

        if(total_eq != 0) {
            optimization(X, x1, x2, final_mat, fin_r, fin_c, total_eq, inva1_b);
        } else {
            optimization(X, x2, x2, final_mat, fin_r, fin_c, total_eq, inva1_b);
        }
    }

    // optimization

    void optimization(t_vector& X,
        t_ivector& x1,
        t_ivector& x2,
        t_matrix& fin_mat,
        int rows,
        int cols,
        int tot_eq,
        t_vector& a1_b)
    {
        t_matrix new_genera;
        t_matrix population;
        t_matrix temp_m;

        t_matrix new_gen;

        t_vector probab;
        t_vector cum_probab;
        t_vector temp_v;

        t_ivector live;
        t_ivector die;

        unsigned long cout_generation = 1;
        unsigned long peak_count;
        float peak_val;
        char response_char;
        int _PROGEND, same;
        int dup_count;
        float Teval;

        int j1, j2, j3, j4, j5, j6, j7;
        int oper;

        int first, first_live, second_live, first_die, second_die;

        if(print_out) {

            cout << "Test case number        : " << test_num << endl;
            cout << "Numbers of operators    : " << P1 << " " << P2 << " " << P3 << " " << P4 << " " << P5 << " " << P6
                 << " " << P7 << endl;
            cout << "Number of generations   : " << generations << endl;
            cout << "Population size         : " << pop_size << endl;
            cout << "Parameter B             : " << B << endl;
            cout << "Parameter Q             : " << Q << endl;
        }

        int P = P1 + P2 + P3 + P4 + P5 + P6 + P7;

        if(P > pop_size) {
            if(print_out)
                cout << "The total number of operators greater than population" << endl;
            return;
        }

        peak_val = 0.0f;
        peak_count = 0;
        int x2_vari = cols - 2;

        resize_matrix(pop_size, x2_vari + 1, population);
        resize_matrix(pop_size, x2_vari + 1, new_genera);
        resize_matrix(2, x2_vari, temp_m);
        resize_vector(pop_size, probab);
        resize_vector(x2_vari, temp_v);
        resize_vector(pop_size, cum_probab);
        resize_vector(pop_size, live);
        resize_vector(pop_size, die);

        if(init_val == 1) {
            if(print_out)
                cout << "USING SINGLE POINT INITIAL POPULATION" << endl;
            _PROGEND = initialize_x2(fin_mat, rows, cols, x1, x2, tot_eq, X, a1_b);
            for(int j = 1; j <= pop_size; j++) {
                for(int i = 1; i <= x2_vari; i++) {
                    population[j][i] = X[x2[i]];
                    population[j][x2_vari + 1] = 0.0f;
                }
                population[j][0] = evaluate(X);
            }
            if(print_out)
                cout << "The initial point of the population is " << endl;
            print_vector(X, 1, tot_eq + x2_vari);
        } else {
            if(print_out)
                cout << "USING MULTIPLE POINT INITIAL POPULATION..." << endl;
            int j = 1;
            while(j <= pop_size) {
                _PROGEND = initialize_x2(fin_mat, rows, cols, x1, x2, total_eq, X, a1_b);
                if(_PROGEND == 1) {
                    for(int i = 1; i <= x2_vari; i++) {
                        population[j][i] = X[x2[i]];
                        population[j][x2_vari + 1] = 0.0f;
                    }
                    j++;
                } else {
                    if(print_out)
                        cout << "Do you wish to include/replicate this vector in the population? (y/n): ";
                    cin >> response_char;
                    if(response_char == 'Y' || response_char == 'y') {
                        do {
                            if(print_out)
                                cout << "How many copies (min. 1, max." << pop_size - j + 1 << ") :" << endl;
                            cin >> dup_count;
                            if((dup_count < 1) || (dup_count > (pop_size - j + 1))) {
                                if(print_out)
                                    cout << "Invalid entry. Must be in the range 1 to " << (pop_size - j + 1) << endl;
                            }
                        } while((dup_count < 1) || (dup_count > (pop_size - j + 1)));

                        for(int k = 1; k <= dup_count; k++) {
                            for(int i = 1; i <= x2_vari; i++) {
                                population[j][i] = X[x2[i]];
                                population[j][x2_vari + 1] = 0.0f;
                            }
                            j++;
                        }
                    }
                }
            }
        }

        // eval. and sort first generation -------------------
        for(int i = 1; i <= pop_size; i++) {
            if(total_eq != 0) {
                for(int j = 1; j <= x2_vari + tot_eq; j++) {
                    X[j] = 0.0;
                }
                find_X(fin_mat, population[i], X, x2_vari, tot_eq, x1, x2, a1_b);
            } else {
                for(int j = 1; j <= x2_vari; j++) {
                    X[j] = population[i][j];
                }
                population[i][0] = evaluate(X);
            }
        }

        /* Evaluate the best in population and assign it to Teval */

        sort_gen(isMinimization, population, pop_size);

        if(print_out)
            cout << "Generation# "
                 << "Solution Value" << endl;

        assign_probab(probab, pop_size, Q);
        find_cum_probab(cum_probab, probab, pop_size);

        Teval = population[1][0];
        peak_val = population[1][0];

        for(int j = 1; j <= pop_size; j++) {
            new_genera[j][x2_vari + 1] = 0.0f;
        }

        /*Reproducing and evaluating for the total number of generations times*/

        // cout<<population[1].size() <<endl;

        do {

            /*
            for(int i = 1; i<pop_size; i++)
            {
                population
            }*/

            /*Initializing the live and die vectors*/
            for(int j = 1; j <= pop_size; j++) {
                live[j] = die[j] = 0.0f;
                for(int i = 0; i <= x2_vari + 1; i++) {
                    new_genera[j][i] = population[j][i];
                }
            }
            /*Finding the agents that will die and the agents that will reproduce*/
            find_live(cum_probab, live, pop_size, P4 + P5 + P7);
            j1 = j2 = j3 = j4 = j5 = j6 = j7 = 0.0;

            while(j1 + j2 + j3 + j4 + j4 + j5 + j5 + j6 + j7 + j7 < P) {
                oper = irange_ran(1, 7);
                switch(oper) {
                case 1:
                    // cout<<"c1"<<endl;
                    /*Applying the first operator, uniform mutation, for the number of times specified*/
                    if(j1 != P1) {
                        do {
                            first = irange_ran(2, pop_size);
                        } while(die[first] == 1);
                        die[first] = 1;

                        new_genera[first][x2_vari + 1] = 1.0f;

                        for(int i = 1; i <= x2_vari; i++) {
                            temp_v[i] = population[first][i];
                        }
                        oper1(temp_v, fin_mat, rows, cols);
                        for(int i = 1; i <= x2_vari; i++) {
                            new_genera[first][i] = temp_v[i];
                        }
                        j1++;
                    }
                    break;
                case 2:
                    // cout<<"c2"<<endl;
                    /*Applying the second operator, boundary mutation, for the number of times specified*/
                    if(j2 != P2) {
                        do {
                            first = irange_ran(2, pop_size);
                        } while(die[first] == 1);

                        die[first] = 1;
                        new_genera[first][x2_vari + 1] = 2.0f;
                        for(int i = 1; i <= x2_vari; i++) {
                            temp_v[i] = population[first][i];
                        }
                        oper2(temp_v, fin_mat, rows, cols);

                        for(int i = 1; i <= x2_vari; i++) {
                            new_genera[first][i] = temp_v[i];
                        }
                        j2++;
                    }

                    break;
                case 3:
                    // cout<<"c3"<<endl;
                    if(j3 != P3) {
                        do {
                            first = irange_ran(2, pop_size);
                        } while(die[first] == 1);
                        die[first] = 1;

                        new_genera[first][x2_vari + 1] = 3.0f;

                        for(int i = 1; i <= x2_vari; i++) {
                            temp_v[i] = population[first][i];
                        }
                        // joao oper3
                        oper3(temp_v, fin_mat, rows, cols, generations, cout_generation, B);
                        for(int i = 1; i <= x2_vari; i++) {
                            new_genera[first][i] = temp_v[i];
                        }
                        j3++;
                    }
                    break;
                case 4:
                    // cout<<"c4"<<endl;
                    /*Applying the fourth operator, whole arithmetical crossover*/
                    if(j4 != (int)P4 / 2) {

                        /*Find two distinct parents for crossover operator 4*/
                        first_live = find_parent(live, pop_size);
                        second_live = find_parent(live, pop_size);
                        same = 1;

                        for(int i = 1; i <= x2_vari; i++) {
                            if(population[first_live][i] != population[second_live][i]) {
                                same = 0;
                            }
                        }
                        if(same == 0) {
                            first_die = find_die(cum_probab, die, pop_size);
                            second_die = find_die(cum_probab, die, pop_size);
                            die[first_die] = 1;
                            die[second_die] = 1;

                            new_genera[first_die][x2_vari + 1] = 4.0f;
                            new_genera[second_die][x2_vari + 1] = 4.0f;

                            for(int i = 1; i <= x2_vari; i++) {
                                temp_m[1][i] = population[first_live][i];
                                temp_m[2][i] = population[second_live][i];
                            }

                            oper4(temp_m[1], temp_m[2], x2_vari);
                            for(int i = 1; i <= x2_vari; i++) {
                                new_genera[first_die][i] = temp_m[1][i];
                                new_genera[second_die][i] = temp_m[2][i];
                            }
                        }
                        j4++;
                    }
                    break;
                case 5:
                    /*Applying the fifth operator, simple arithmetical crossover*/
                    if(j5 != (int)P5 / 2) {
                        first_live = find_parent(live, pop_size);
                        second_live = find_parent(live, pop_size);
                        same = 1;
                        for(int i = 1; i <= x2_vari; i++) {
                            if(population[first_live][i] != population[second_live][i]) {
                                same = 0;
                            }
                        }

                        if(same == 0) {
                            first_die = find_die(cum_probab, die, pop_size);
                            second_die = find_die(cum_probab, die, pop_size);
                            die[first_die] = 1;
                            die[second_die] = 1;
                            new_genera[first_die][x2_vari + 1] = 5.0;
                            new_genera[second_die][x2_vari + 1] = 5.0;
                            for(int i = 1; i <= x2_vari; i++) {
                                temp_m[1][i] = population[first_live][i];
                                temp_m[2][i] = population[second_live][i];
                            }
                            oper5(temp_m[1], temp_m[2], STEP, rows, cols, fin_mat);
                            for(int i = 1; i <= x2_vari; i++) {
                                new_genera[first_die][i] = temp_m[1][i];
                                new_genera[second_die][i] = temp_m[2][i];
                            }
                        }
                        j5++;
                    }
                    break;
                case 6:
                    // cout<<"c6"<<endl;
                    /*Applying the sixth operator, whole non-uniform mutation, for the number of times specified*/
                    if(j6 != P6) {
                        do {
                            first = irange_ran(2, pop_size);
                        } while(die[first] == 1);
                        die[first] = 1;

                        new_genera[first][x2_vari + 1] = 6.0;
                        for(int i = 1; i <= x2_vari; i++) {
                            temp_v[i] = population[first][i];
                        }
                        oper6(temp_v, fin_mat, rows, cols, generations, cout_generation, B);
                        for(int i = 1; i <= x2_vari; i++) {
                            new_genera[first][i] = temp_v[i];
                        }
                        j6++;
                    }
                    break;
                case 7:
                    // cout<<"c7"<<endl;
                    /*Applying the seventh operator*/
                    if(j7 != (int)P7 / 2) {
                        /*Find two distinct parents for operator 7*/
                        first_live = find_parent(live, pop_size);
                        second_live = find_parent(live, pop_size);
                        same = 1;
                        for(int i = 1; i <= x2_vari; i++) {
                            if(population[first_live][i] != population[second_live][i]) {
                                same = 1;
                            }
                        }
                        if(same == 0) {
                            first_die = find_die(cum_probab, die, pop_size);
                            die[first_die] = 1;
                            new_genera[first_die][x2_vari + 1] = 7.0;
                            for(int i = 1; i <= x2_vari; i++) {
                                if(first_live < second_live) {
                                    temp_m[2][i] = population[first_live][i];
                                    temp_m[1][i] = population[second_live][i];
                                } else {
                                    temp_m[2][i] = population[second_live][i];
                                    temp_m[1][i] = population[first_live][i];
                                }
                            }
                            oper7(temp_m[1], temp_m[2], rows, cols, fin_mat);
                            for(int i = 1; i <= x2_vari; i++) {
                                new_genera[first_die][i] = temp_m[1][i];
                            }
                        }
                        j7++;
                    }
                    break;
                }
            }

            /*Replace the population with the new generation */
            new_gen = new_genera;
            new_genera = population;
            population = new_gen;

            // para mostrar los datos de evolucion:--------------------------------------
            if(print_evolution) {
                cout << "pop_generation " << pop_size << " " << cout_generation << endl;
                for(int i = 1; i < population.size(); i++) {
                    find_X(fin_mat, population[i], X, x2_vari, tot_eq, x1, x2, a1_b);
                    cout << "values " << i << " ";
                    for(int j = 1; j <= x2_vari + tot_eq; j++) {
                        // cout<< X[j]<<endl;
                        cout << setprecision(9) << X[j] << " ";
                        // printf("%18.8f\n", j, X[j]);
                    }
                    cout << endl;
                }
            }
            //---------------------------------------------------------------------------

            /*Evaluate each of the agents in the new population*/
            for(int i = 1; i <= pop_size; i++) {
                if(tot_eq != 0) {
                    for(int j = 1; j <= x2_vari + tot_eq; j++) {
                        X[j] = 0.0;
                    }
                    find_X(fin_mat, population[i], X, x2_vari, tot_eq, x1, x2, a1_b);
                } else
                    for(int j = 1; j <= x2_vari; j++)
                        X[j] = population[i][j];

                population[i][0] = evaluate(X);
            }

            sort_gen(isMinimization, population, pop_size);

            switch(isMinimization) {
            case 0:
                if(Teval > population[1][0]) {
                    Teval = population[1][0];

                    if(print_out)
                        cout << cout_generation << " " << population[1][0] << endl;

                    peak_count = cout_generation;
                    peak_val = population[1][0];
                }
                break;
            case 1:
                if(Teval < population[1][0]) {
                    Teval = population[1][0];
                    if(print_out)
                        cout << cout_generation << " " << population[1][0] << endl;
                    peak_count = cout_generation;
                    peak_val = population[1][0];
                }
                break;
            }
            cout_generation++;
        } while(cout_generation <= generations);

        if(print_out)
            printf("\nBest solution was found at generation %lu (solution value = %.8f)\n", peak_count, peak_val);
        if(print_out)
            printf("\n\nBest solution found:\n\n");

        if(!print_evolution) {
            printf("%.32f\n", peak_val);
        } else {
        }

        /* print best solution */
        find_X(fin_mat, population[1], X, x2_vari, tot_eq, x1, x2, a1_b);
        for(int j = 1; j <= x2_vari + tot_eq; j++) {
            if(print_out)
                printf(" X[%2d] :\t%18.8f\n", j, X[j]);
        }

        // para mostrar los datos de evolucion:--------------------------------------
        if(print_evolution) {
            cout << "optimal ";
            for(int j = 1; j <= x2_vari + tot_eq; j++) {
                // cout<< X[j]<<endl;
                cout << setprecision(9) << X[j] << " ";
                // printf("%18.8f\n", j, X[j]);
            }
        }
    }

    // operators

    void oper7(t_vector& p1, t_vector& p2, int rows, int cols, t_matrix& fin_mat)
    {
        t_vector child;
        int _CHECK = 0; /*Check to see if the newly created vector satisfies the*/
                        /*set of constraints*/
        int i, n = 2, tries = 10;
        float A;

        resize_vector(cols - 2, child);

        do {
            A = frange_ran(0.0, 1.0);
            for(i = 1; i <= cols - 2; i++)
                child[i] = A * (p2[i] - p1[i]) + p2[i];

            /*Check to see if it satisfies the constraints*/
            _CHECK = satis_con(child, fin_mat, rows, cols);
            n++;
            /*If the constraints not satisfied, then try again */
        } while((n <= tries) && (_CHECK == 0));

        if(_CHECK == 1)
            for(i = 1; i <= cols - 2; i++)
                p1[i] = child[i];
        child.clear();
    }

    void oper6(t_vector& parent, t_matrix& fin_mat, int rows, int cols, unsigned long T, unsigned long t, int B)
    {
        int comp, i, num;
        t_ivector next;
        float llim, ulim;

        resize_vector(cols - 2, next);

        for(i = 1; i <= cols - 2; i++)
            next[i] = 0;

        for(i = 1; i <= cols - 2; i++) {
            do
                comp = irange_ran(1, cols - 2);
            while(next[comp] == 1);
            next[comp] = 1;

            find_range(llim, ulim, comp, fin_mat, rows, cols, parent);

            /*From the current value of the component to be mutated, chooose at random*/
            /*whether to mutate with a lesser value or a greater value*/
            /*Then find a value lesser or greater than the original value from the*/
            /*function get_f()*/
            parent[comp] = (flip() == TAIL) ? parent[comp] - get_F(T, t, parent[comp] - llim, B) :
                                              parent[comp] + get_F(T, t, ulim - parent[comp], B);
        }
        next.clear();
    }

    void oper5(t_vector& p1, t_vector& p2, int STEP, int rows, int cols, t_matrix& fin_mat)
    {

        t_matrix child;
        int _CHECK1 = 0, /*Check to see if the newly created vectors satisfies the*/
            _CHECK2 = 0; /*set of constraints*/
        int i, n = 1, cut;

        resize_matrix(2, cols - 2, child);

        /*Get a random spot on the vector for crossover*/
        cut = irange_ran(1, cols - 2);
        /*Copy the parent vectors on to the child vectors*/
        for(i = 1; i <= cut; i++) {
            child[1][i] = p1[i];
            child[2][i] = p2[i];
        }
        do {
            /*Cross the two vectors*/
            for(i = cut + 1; i <= cols - 2; i++) {
                child[1][i] = p1[i] * (float)n / (float)STEP + p2[i] * (1.0 - (float)n / (float)STEP);
                child[2][i] = p2[i] * (float)n / (float)STEP + p1[i] * (1.0 - (float)n / (float)STEP);
            }

            /*Check to see if they satisfy the constraints*/
            _CHECK1 = satis_con(child[1], fin_mat, rows, cols);
            _CHECK2 = satis_con(child[2], fin_mat, rows, cols);
            n++;
            /*If the constraints not satisfied, then generate another*/
            /*set of crossed over values*/
        } while((n <= STEP) && ((_CHECK1 == 0) || (_CHECK2 == 0)));

        for(i = 1; i <= cols - 2; i++) {
            p1[i] = child[1][i];
            p2[i] = child[2][i];
        }

        child.clear();
    }

    int satis_con(t_vector& child, t_matrix& fin_mat, int rows, int cols)
    {
        int i, j;
        float tot;

        for(j = 1; j <= cols - 2; j++)
            if((child[j] > fin_mat[j][cols]) || (child[j] < fin_mat[j][1]))
                return 0;

        /*Substitute the values for all the variables, with the new values of*/
        /*the vector passed and check if the limits are crossed*/
        for(j = 1; j <= rows; j++) {
            tot = 0.0;
            for(i = 2; i <= cols - 1; i++)
                tot = tot + fin_mat[j][i] * child[i - 1];

            /*If the limits are crossed, return FALSE*/
            if((tot < fin_mat[j][1]) || (tot > fin_mat[j][cols]))
                return 0;
        }

        /*If the limits are not crossed, return TRUE*/
        return 1;
    }

    void oper4(t_vector& p1, t_vector& p2, int x2_vari)
    {
        t_matrix child;
        int i;
        float A;

        resize_matrix(2, x2_vari, child);

        do {
            A = frange_ran(0.0, 1.0);
        } while(A == 0); /* insure A is above 0 */

        for(i = 1; i <= x2_vari; i++) {
            child[1][i] = p1[i] * A + p2[i] * (1.0 - A);
            child[2][i] = p2[i] * A + p1[i] * (1.0 - A);
        }
        for(i = 1; i <= x2_vari; i++) {
            p1[i] = child[1][i];
            p2[i] = child[2][i];
        }

        child.clear();
    }

    int find_die(t_vector& cum_probab, t_ivector& die, int pop_size)
    {

        float random;
        int i;
        bool done = false;

        do {

            /*Choosing a random cumulative probability*/
            random = frange_ran(0.0, 1.0);
            i = 0;
            /*Finding the agent with the chosen cumulative probability*/
            do {
                // cout<<"do 2"<<endl;
                i++;
            } while((random > cum_probab[i]) && (i < pop_size));

            /*Chosing the agent to die*/
            if((die[pop_size + 1 - i] == 0) && (i < pop_size)) {
                done = true;
            }
        } while(!done);
        return (pop_size + 1 - i);
    }

    int find_parent(t_ivector& live, int pop_size)
    {
        int i, temp, t1, tot = 0;

        /*Finding the total number of parents to reproduce*/
        for(i = 1; i <= pop_size; i++)
            tot = tot + live[i];
        if(tot == 0) {
            if(print_out)
                cout << "No agents to select" << endl;
            return 0;
        }

        /*Choosing one of them randomly*/
        temp = irange_ran(1, tot);

        tot = 0;
        i = 1;
        do {
            if(live[i] != 0)
                t1 = i;
            tot = tot + live[i++];
        } while(tot < temp);

        /*Decrementing the number of times the parent chosen is going to reproduce*/
        live[t1]--;
        return (t1);
    }

    void oper3(t_vector& parent, t_matrix& fin_mat, int rows, int cols, unsigned long T, unsigned long t, int B)
    {
        int comp, i;
        float llim, ulim;

        comp = irange_ran(1, cols - 2);
        find_range(llim, ulim, comp, fin_mat, rows, cols, parent);

        /*From the current value of the component to be mutated, chooose at random*/
        /*whether to mutate with a lesser value or a greater value*/
        /*Then find a value lesser or greater than the original value from the*/
        /*function get_f()*/
        parent[comp] = (flip() == TAIL) ? parent[comp] - get_F(T, t, parent[comp] - llim, B) :
                                          parent[comp] + get_F(T, t, ulim - parent[comp], B);
    }

    float get_F(unsigned long T, unsigned long t, float y, int B)
    {
        float factor;

        factor = (float)pow(1.0 - (float)t / (float)T, (float)B);
        factor = factor * frange_ran(0.0, 1.0);
        if(factor < 0.00001)
            factor = 0.00001;
        return (y * factor);
    }

    void oper2(t_vector& parent, t_matrix& fin_mat, int rows, int cols)
    {
        int comp, i;
        float llim, ulim; /*Lower and Upper limits of the value to be mutated*/

        comp = irange_ran(1, cols - 2);

        /*Finding the lower and upper limits between which the values are to be mutated*/
        find_range(llim, ulim, comp, fin_mat, rows, cols, parent);

        /*Replace either the lower limit or the upper limit at random,*/
        /*for the old value*/
        parent[comp] = (flip() == TAIL) ? llim : ulim;
    }

    void oper1(t_vector& parent, t_matrix& fin_mat, int rows, int cols)
    {
        int comp, i;
        float llim, ulim; /*Lower and Upper limits of the value to be mutated*/

        comp = irange_ran(1, cols - 2);

        /*Finding the lower and upper limits between which the values are to be mutated*/
        find_range(llim, ulim, comp, fin_mat, rows, cols, parent);

        /*Find a random value between the lower and the upper limits, to substitute*/
        /*for the old value*/
        parent[comp] = frange_ran(llim, ulim);
    }

    void find_range(float& llim, float& ulim, int comp, t_matrix& fin_mat, int rows, int cols, t_vector& parent)
    {
        int i, j;
        float tot, templ, tempu, temp;
        int _CHANGE = 0;

        /*Replace the values of the variables, and for the values of the vectors*/
        /*except the one which is to be mutated*/
        /*Find the lower and upper limits within which the mutation component's*/
        /*value should be found*/
        for(j = 1; j <= rows; j++)
            if(fin_mat[j][comp + 1] != 0.0) {
                tot = 0.0;
                for(i = 2; i <= cols - 1; i++)
                    if(i != comp + 1)
                        tot = tot + fin_mat[j][i] * parent[i - 1];

                templ = (fin_mat[j][1] - tot) / fin_mat[j][comp + 1];
                tempu = (fin_mat[j][cols] - tot) / fin_mat[j][comp + 1];

                if(fin_mat[j][comp + 1] < 0)
                    swap(templ, tempu);

                if(!_CHANGE) {
                    llim = templ;
                    ulim = tempu;
                    _CHANGE = 1;
                } else {
                    if(llim < templ)
                        llim = templ;
                    if(ulim > tempu)
                        ulim = tempu;
                }
            }
    }

    // genocop functions
    void find_live(t_vector& cum_probab, t_ivector& live, int pop_size, int P)
    {
        float random;
        int count = 0, /*Count of the number of agents chosen to live*/
            i;

        do {
            /*Choosing a random cumulative probability*/
            random = frange_ran(0.0, 1.0);
            i = 0;
            /*Finding the agent with the chosen cumulative probability*/
            do {
                i++;
            } while((random > cum_probab[i]) && (i < pop_size));

            /*Chosing the parent with that probability to reproduce*/
            if(count < P) {
                live[i]++;
                count++;
            }
        } while(count < P);
    }

    void find_cum_probab(t_vector& cum_probab, t_vector& probab, int pop_size)
    {
        int i;
        cum_probab[1] = probab[1];
        for(i = 2; i <= pop_size; i++)
            cum_probab[i] = cum_probab[i - 1] + probab[i];
    }

    void assign_probab(t_vector& probab, int pop_size, double Q)
    {
        int i;
        double Q1;

        /* Q, Q(1-Q)^1, Q(1-Q)^2 ... Q(1-Q)^n */
        Q1 = Q / (1 - pow(1 - Q, pop_size));
        for(i = 1; i <= pop_size; i++)
            probab[i] = Q1 * pow(1 - Q, i - 1);
    }

    void sort_gen(int MinMax, t_matrix& population, int pop_size)
    {
        int i, j, k;
        switch(MinMax) {
        case 0:
            for(i = 1; i <= pop_size; i++)
                for(j = i + 1; j <= pop_size; j++)
                    if(population[i][0] > population[j][0])
                        swap_gen(&population[i], &population[j]);
            break;

        case 1:
            for(i = 1; i <= pop_size; i++)
                for(j = i + 1; j <= pop_size; j++)
                    if(population[i][0] < population[j][0])
                        swap_gen(&population[i], &population[j]);
            break;
        default:
            if(print_out)
                cout << "Incorrect data: Must be a 0 or 1" << endl;
            return;
        }
    }

    void swap_gen(t_vector* x, t_vector* y)
    {
        t_vector temp;
        temp = *x;
        *x = *y;
        *y = temp;
    }

    void swap_gen(float& x, float& y)
    {
        float temp;
        temp = x;
        x = y;
        y = temp;
    }

    // print
    void print_vector(t_vector& arr, int l, int u)
    {
        for(int i = l; i <= u; i++) {
            if(print_out)
                cout << arr[i] << endl;
        }
    }
    // Fitness function

    void find_X(t_matrix& final_mat,
        t_vector& agent,
        t_vector& X,
        int x2_vari,
        int tot_eq,
        t_ivector& x1,
        t_ivector& x2,
        t_vector& a1_b)
    {
        int i, j;

        for(j = 1; j <= x2_vari; j++)
            X[x2[j]] = agent[j];

        for(j = 1; j <= tot_eq; j++) {
            X[x1[j]] = a1_b[j];
            for(i = 1; i <= x2_vari; i++)
                X[x1[j]] = X[x1[j]] + agent[i] * final_mat[j + x2_vari][i + 1];
        }
    }

    // Random----------------

    int irange_ran(int llim, int ulim)
    {
        int num;

        do
            num = llim + ((int)((newrand() * (long)(ulim - llim + 1)) / (long)65535));
        while((num < llim) || (num > ulim));
        return (num);
    }

    void seed()
    {
        long l_time;
        int n;
        unsigned int x;

        l_time = time(NULL);
        l_time = l_time % 65536;
        rseed = l_time;

        for(n = 0; n < SHUFFLE; n++) /* initialize random number array */
        {
            rand_array[n] = randint();
        }
        for(n = 0; n < 1000; n++) /* warm up the generator */
        {
            x = newrand();
        }
    }
    unsigned int randint()
    {
        int num;

        rseed = (MULT * rseed + INCR) % MOD;
        num = rseed % 65536;
        return (num);
    }

    unsigned int newrand()
    {
        int r_offset;
        unsigned int ran_int;
        r_offset = (int)(randint() % SHUFFLE);
        ran_int = rand_array[r_offset];
        rand_array[r_offset] = randint();
        return (ran_int);
    }

    float frange_ran(float llim, float ulim)
    {
        float diff, num1;

        diff = ulim - llim;
        if(diff == 0) {
            return llim;
        } else if(diff < 0.0001) {
            return ((flip() == TAIL) ? llim : ulim);
        }
        do {
            num1 = llim + (float)((newrand() * (ulim - llim)) / 65535);
        } while((num1 < llim) || (num1 > ulim));
        return num1;
    }

    int initialize_x2(t_matrix& final_mat,
        int rows,
        int cols,
        t_ivector& x1,
        t_ivector& x2,
        int x1_vari,
        t_vector& X,
        t_vector& a1_b)
    {
        int x2_vari = cols - 2;
        int itry = 0;
        int _LOW;
        int _HIGH;
        t_vector temp;
        t_matrix trymat;
        double sum = 0.0;
        int _CHECK = 1;
        int LOCAL_CHECK = 1;

        resize_vector(x2_vari, temp);
        resize_matrix(rows - x2_vari, x2_vari, trymat);
        LOCAL_CHECK = 1;
        _CHECK = 1;

        do {
            if(itry < TRIES) {
                ++itry;
            }
            for(int i = 1; i <= x2_vari; i++) {
                temp[i] = frange_ran(final_mat[i][1], final_mat[i][cols]);
            }
            if(x2_vari != rows) {
                for(int j = x2_vari + 1; j <= rows; j++) {
                    sum = 0.0f;
                    for(int i = 1; i <= x2_vari; i++) {
                        sum = sum + temp[i] * final_mat[j][i + 1];
                    }
                    if((!_LOW) || (!_HIGH)) {
                        break;
                    }
                }
            } else {
                _LOW = _HIGH = 1;
            }
        } while((itry < TRIES) && ((!_LOW) || (!_HIGH)));

        if(itry >= TRIES) {
            _CHECK = 0;
            if(print_out)
                cout << "No feasible random solution found after " << TRIES << endl;
            do {
                LOCAL_CHECK = 1;
                if(print_out)
                    cout << "Please input initial values" << endl;
                for(int i = 0; i <= x1_vari + x2_vari; i++) {
                    cout << i << endl;
                    cin >> X[i];
                }

                for(int j = 1; j <= rows; j++) {
                    sum = 0.0f;
                    for(int i = 1; i <= x2_vari; i++) {
                        sum = sum + X[x2[i]] * final_mat[j][i + 1];
                    }

                    _LOW = (sum >= final_mat[j][1]) ? 1 : 0;
                    _HIGH = (sum <= final_mat[j][cols]) ? 1 : 0;

                    if((!_LOW) || (!_HIGH)) {
                        if(print_out)
                            cout << "The input values do not satisfy the constraint #" << j << endl;
                        LOCAL_CHECK = 0;
                    }
                }
            } while(LOCAL_CHECK == 0);
        } else {
            if(x1_vari != 0) {
                for(int i = 1; i <= x2_vari; i++) {
                    X[x2[i]] = temp[i];
                }
                for(int j = 1; j <= x1_vari; j++) {
                    X[x1[j]] = a1_b[j];
                    for(int i = 1; i <= x2_vari; i++) {
                        X[x1[j]] = X[x1[j]] + temp[i] * final_mat[j + x2_vari][i + 1];
                    }
                }
            } else {
                for(int i = 1; i <= x2_vari; i++) {
                    X[i] = temp[i];
                }
            }
        }

        temp.clear();
        trymat.clear();

        return _CHECK;
    }

    // Change order----------
    void find_final_mat1(t_vector& l2, t_vector& u2, t_matrix& finmat, int row, int col)
    {
        int j = 2;
        for(int i = 1; i <= row; i++) {
            finmat[i][1] = l2[i];
            finmat[i][col] = u2[i];
            finmat[i][j++] = 1.0f;
        }
    }

    void find_final_mat2(t_matrix& newin, int r, int c, int finr, t_matrix& finmat)
    {
        for(int i = 1; i <= r; i++) {
            for(int j = 1; j <= c; j++) {
                finmat[finr][j] = newin[i][j];
            }
            finr++;
        }
    }

    void find_final_mat3(t_matrix& orgin, int r, int c, int finr, t_matrix& finmat)
    {
        for(int i = 1; i <= r; i++) {
            finmat[finr][1] = MIN;
            for(int j = 1; j <= c; j++) {
                finmat[finr][j + 1] = orgin[i][j];
            }
            finr++;
        }
    }

    void initialize(t_matrix& mat, int rows, int cols)
    {
        for(int i = 1; i <= rows; i++) {
            for(int j = 1; j <= cols; j++) {
                mat[i][j] = 0.0f;
            }
        }
    }

    void find_org_in_eq(t_vector& a1_b,
        t_matrix& a1_a2,
        t_vector& vec_d,
        t_matrix& c1,
        t_matrix& c2,
        int c1row,
        int a1a2_row,
        int a1a2_col,
        t_matrix& org_ineq)
    {
        t_vector temp;
        t_matrix mat;

        resize_vector(c1row, temp);
        resize_matrix(c1row, a1a2_col - 1, mat);

        mvprod(c1row, a1a2_row, temp, c1, a1_b);
        mmprod(c1row, a1a2_row, a1a2_col - 1, mat, c1, a1_a2);

        for(int i = 1; i <= c1row; i++) {
            for(int j = 1; j <= a1a2_col; j++) {
                if(j == a1a2_col) {
                    org_ineq[i][j] = vec_d[i] - temp[i];
                } else {
                    org_ineq[i][j] = c2[i][j] - mat[i][j];
                }
            }
        }

        temp.clear();
        mat.clear();
    }

    void find_new_in_eq(t_vector& a1b, t_matrix& a1a2, t_vector& ll, t_vector& ul, int rows, int cols, t_matrix& newin)
    {
        for(int i = 1; i <= rows; i++) {
            for(int j = 1; j <= cols; j++) {
                if(j == 1) {
                    newin[i][j] = ll[i] - a1b[i];
                } else if(j == cols) {
                    newin[i][j] = ul[i] - a1b[i];
                } else {
                    newin[i][j] = 0 - a1a2[i][j - 1];
                }
            }
        }
    }

    void find_lu1_lu2(int total_variables,
        int total_eq,
        t_ivector& x1,
        t_ivector& x2,
        t_vector& dom,
        t_vector& dom1,
        t_vector& dom2)
    {
        for(int i = 1; i <= total_eq; i++) {
            dom1[i] = dom[x1[i]];
        }

        for(int i = 1; i <= total_variables - total_eq; i++) {
            dom2[i] = dom[x2[i]];
        }
    }

    void find_limits(int tot, t_matrix& domains, t_vector& llim, t_vector& ulim)
    {
        for(int i = 1; i <= tot; i++) {
            llim[i] = domains[i][1];
            ulim[i] = domains[i][3];
        }
    }

    void find_ac1_ac2(int t1, int t2, int t3, t_ivector& x1, t_ivector& x2, t_matrix& mat, t_matrix& ac1, t_matrix& ac2)
    {
        for(int i = 1; i <= t1; i++) {
            for(int j = 1; j <= t2; j++) {
                ac1[j][i] = mat[j][x1[i]];
            }
        }

        for(int i = 1; i <= t3; i++) {
            for(int j = 1; j <= t2; j++) {
                ac2[j][i] = mat[j][x2[i]];
            }
        }
    }

    void find_x1_x2(int tot, t_imatrix& var_order, t_ivector& x1, t_ivector& x2)
    {
        int j = 1;
        int k = 1;
        for(int i = 1; i <= tot; i++) {
            if(var_order[i][2] == 1) {
                x1[j++] = var_order[i][1];
            } else {
                x2[k++] = var_order[i][1];
            }
        }
    }

    void get_var_order(int total_variables, int total_eq, t_ivector& cart, t_imatrix& var_order)
    {
        int in_val;
        for(int i = 1; i <= total_variables; i++) {
            var_order[i][1] = i;
            var_order[i][2] = 0;
        }
        if(print_out)
            cout << total_eq << " variables can be eliminated out of " << total_variables << endl;
        if(print_out)
            cout << "Please enter the chosen ones below. " << endl << endl;
        for(int i = 1; i <= total_eq; i++) {
            do {
                if(print_out)
                    cout << "Enter index (subscript) of variable" << i << " to be eliminated : ";
                cin >> in_val;
                if(in_val < 1 || in_val > total_variables) {
                    if(print_out)
                        cout << "Invalid subscript.  Must be in the range 1 to " << total_variables << endl;
                }
            } while(in_val < 1 || in_val > total_variables);
            if(print_out)
                cout << "Value " << in_val << " accepted" << endl;
            var_order[in_val][2] = 1;
        }
    }

    // MATRIX AND VECTORS OPERATIONS

    void mmprod(int m, int nm, int n, t_matrix& mul_cm, t_matrix& mul_am, t_matrix& mul_bm)
    {
        for(int i = 1; i <= m; i++) {
            for(int j = 1; j <= n; j++) {
                mul_cm[i][j] = 0.0f;
                for(int k = 1; k < nm + 1; k++) {
                    mul_cm[i][j] = mul_am[i][k] * mul_bm[k][j] + mul_cm[i][j];
                }
            }
        }
    }

    void mvprod(int m, int nm, t_vector& cm, t_matrix& am, t_vector& bm)
    {
        for(int i = 1; i <= m; i++) {
            cm[i] = 0.0f;
            for(int k = 1; k <= nm; k++) {
                cm[i] = cm[i] + am[i][k] * bm[k];
            }
        }
    }

    void inverse(t_matrix& a, t_matrix& a_inverse, int n)
    {
        int n1, k1, l1;
        float d1, v, w;
        t_matrix b;
        resize_matrix(n, n, b);

        for(int i = 1; i <= n; i++) {
            for(int j = 1; j <= n; j++) {
                b[i][j] = 0.0f;
            }
        }

        d1 = det(a, n);

        if(d1 == 0.0f) {
            if(print_out)
                cout << "Variables are dependent : No inverse exists." << endl;
            return;
        } else {
            v = -1.0;
            for(int i = 1; i <= n; i++) {
                v = -v;
                w = -1.0f;
                for(int j = 1; j <= n; j++) {
                    w = -w;
                    n1 = n - 1;
                    for(int k = 1; k <= n1; k++) {
                        k1 = k;
                        if(k >= i) {
                            k1 = k + 1;
                        }

                        for(int l = 1; l <= n1; l++) {
                            l1 = l;
                            if(l >= j) {
                                l1 = l + 1;
                            }
                            b[k][l] = a[k1][l1];
                        }
                    }
                    a_inverse[j][i] = det(b, n1) / d1 * v * w;
                }
            }
        }
        b.clear();
    }

    float det(t_matrix& input_matrix, int nl)
    {
        int ia, ib, i1, j1, n2, nz, n3, l, l1, m, m1;
        float d, x;
        t_matrix temp_input_matrix, b1;

        resize_matrix(nl, nl, temp_input_matrix);
        resize_matrix(nl, nl, b1);

        for(i1 = 1; i1 <= nl; i1++) {
            for(j1 = 1; j1 <= nl; j1++) {
                temp_input_matrix[i1][j1] = input_matrix[i1][j1];
            }
        }

        d = 1.0f;

        for(ia = 1; ia <= nl; ia++) {
            x = -1.0;
            n2 = nl - ia + 1;
            n3 = n2 - 1;
            nz = 0;
            for(ib = 1; ib <= n2; ib++) {
                if(temp_input_matrix[ib][1] != 0.0f) {
                    nz = ib;
                }
            }
            if(nz == 0.0f) {
                d = 0.0f;
            } else {
                for(i1 = 1; i1 <= nz; i1++) {
                    x = -x;
                }

                d *= temp_input_matrix[nz][1] * x;

                for(l = 1; l <= n3; l++) {
                    l1 = l;
                    if(l >= nz) {
                        l1 = l + 1;
                    }
                    for(m = 1; m <= n3; m++) {
                        m1 = m + 1;
                        b1[l][m] = temp_input_matrix[l1][m1] -
                            temp_input_matrix[l1][1] / temp_input_matrix[nz][1] * temp_input_matrix[nz][m1];
                    }
                }
            }
            for(i1 = 1; i1 <= n3; i1++) {
                for(j1 = 1; j1 <= n3; j1++) {
                    temp_input_matrix[i1][j1] = b1[i1][j1];
                }
            }
        }

        b1.clear();
        temp_input_matrix.clear();

        return d;
    }

    //------------------

    Genocop(const string& path)
    {
        create_from_file(path);
    }
};