#pragma once
#include <iostream>
#include <vector>

#define MIN -32768
#define MAX 32768

using namespace std;

class Genocop
{
public:
    typedef vector<vector<float> > t_matrix;
    typedef vector<float> t_vector;
    typedef vector<vector<int> > t_imatrix;
    typedef vector<int> t_ivector;

    int total_variables;
    int total_eq;
    int total_ineq;
    int total_domains;

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
            cout << "No se puede leer el archivo" << endl;
            return;
        }
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
        }
    }

    //-----------------//
    
    void find_org_in_eq(t_vector& a1_b, )
    {
        
    }
    
    void find_new_in_eq( t_vector &a1b, t_matrix &a1a2, t_vector &ll, t_vector &ul, int rows, int cols, t_matrix &newin )
    {
        for( int i = 1 ; i<=rows; i++)
        {
            for( int j = 1; j<=cols; j++ )
            {
                if( j==1 )
                {
                    newin[i][j] = ll[i] - a1b[i];
                }
                else if( j == cols )
                {
                    newin[i][j] = ul[i] - a1b[i];
                }
                else
                {
                    newin[i][j] = 0 - a1a2[i][j-1];
                }
            }
        }
    }

    void find_lu1_lu2(int total_variables,
        int total_eq,
        t_ivector &x1,
        t_ivector &x2,
        t_vector &dom,
        t_vector &dom1,
        t_vector &dom2)
    {
        for( int i = 1; i<=total_eq; i++ )
        {
            dom1[i] = dom[ x1[i] ];
        }
        
        for( int i = 1; i<=total_variables - total_eq; i++ )
        {
            dom2[i] = dom[ x2[i] ];
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
        cout << total_eq << " variables can be eliminated out of " << total_variables << endl;
        cout << "Please enter the chosen ones below. " << endl << endl;
        for(int i = 1; i <= total_eq; i++) {
            do {
                cout << "Enter index (subscript) of variable" << i << " to be eliminated : ";
                cin >> in_val;
                if(in_val < 1 || in_val > total_variables) {
                    cout << "Invalid subscript.  Must be in the range 1 to " << total_variables << endl;
                }
            } while(in_val < 1 || in_val > total_variables);

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