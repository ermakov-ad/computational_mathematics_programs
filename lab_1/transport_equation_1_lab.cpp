// solution of transport equation
// U'_t + lambda * U'_x = 0
// scheme:
// n+1: -   -   - * -   -
// n  : - * - * - * -   -
// n-1: -   -   -   - * -
//       m-2 m-1  m  m+1
// A*U(n, m-2) + B*U(n, m-1) + C*U(n, m) + D*U(n-1, m+1) = U(n+1, m)

#include <iostream>
#include <fstream>
#include <math.h>
#include <exception>

const double cfl = 0.2;                 // Current number
const double lambda = 1.0;
const double h = 0.01;                  // coordinate step
const double tau = cfl * h / lambda;    // time step
const double left_bound = 0.0;
const double right_bound = 2.0;
const double start_left_bound = 0.4;    // place with difference in start conditions
const double start_right_bound = 0.6;
std::string default_file_name = "transport_equation.txt";

struct scheme_coeffs {
    double A;
    double B;
    double C;
    double D;
};

double initial_conditions (double x) {
    double ans;
    if (x >= start_left_bound && x <= start_right_bound)
    {
        ans = sqrt(1.0 - 100.0*(x - 0.5)*(x - 0.5));
    }
    else
    {
        ans = 0.0;
    }
    return ans;
}

double analytical_solution (double x, double t) {
    return initial_conditions(x - lambda * t);
}

double left_boundary_condition (double t) {
    return 0.0;     // phi(0, t) = 0
}

double right_boundary_condition (double t) {
    return 0.0;     // phi(L, t) = 0
}

class transport_equation {

    public:
        transport_equation(double left_boundary, double right_boundary, unsigned number_of_cells, double time_step, bool save_txt_file, unsigned number_of_schemes);
        ~transport_equation();
        void set_initial_layers(double* first_layer, double* second_layer);
        void calculate_process(double time, scheme_coeffs c);       //calculate all process with fixed coefficients from structure c
        void print_results();

        void set_scheme_coeffs(scheme_coeffs* c);
        void calculate_process_hybrid(double time);                 //calculate process with a hybrid scheme - combine some schemes, which takes from previous method
        void enter_file_name_for_saving(std::string str);

    private:
        void change_layer();
        double calc_new_value(scheme_coeffs c, unsigned ind);
        bool check_monotony(double check_value, double present_value, double previous_value);

        bool save_txt_layers;
        unsigned count_of_schemes;
        unsigned count_of_cells;
        double l_bound;
        double r_bound;
        double t_step;
        scheme_coeffs* coeffs;
        bool entered_coeffs = false;
        std::string save_file_name = default_file_name;

        // arrays with state in 3 times level
        double* previous_layer;
        double* present_layer;
        double* next_layer;
};

transport_equation::transport_equation(double left_boundary, double right_boundary, unsigned number_of_cells, double time_step, bool save_txt_file, unsigned number_of_schemes) {
    l_bound = left_boundary;
    r_bound = right_boundary;
    t_step = time_step;
    save_txt_layers = save_txt_file;
    count_of_schemes = number_of_schemes;
    count_of_cells = number_of_cells;

    previous_layer = new double[count_of_cells];
    present_layer = new double[count_of_cells];
    next_layer = new double[count_of_cells];
    coeffs = new scheme_coeffs[number_of_schemes];
    std::cout << "class exemplar was created" << std::endl;
}

transport_equation::~transport_equation() {
    delete[] previous_layer;
    delete[] present_layer;
    delete[] next_layer;
    delete[] coeffs;
    std::cout << "class exemplar was deleted" << std::endl;
}

void transport_equation::set_initial_layers(double* first_layer, double* second_layer) {
    for (unsigned i = 0; i < count_of_cells; ++i) {
        previous_layer[i] = first_layer[i];
        present_layer[i] = second_layer[i];
    }
    std::cout << "first and second layers was inicialised" << std::endl;
}

void transport_equation::set_scheme_coeffs(scheme_coeffs* c) {
    for (unsigned i = 0; i < count_of_schemes; ++i) {
        coeffs[i] = c[i];
    }
    entered_coeffs = true;
    std::cout << "scheme coefficients was inicialised, hybrid scheme unlocked" << std::endl;
}

double transport_equation::calc_new_value(scheme_coeffs c, unsigned ind) {
    return c.A*present_layer[ind - 2] + c.B*present_layer[ind - 1] + c.C*present_layer[ind] + c.D*previous_layer[ind + 1];
}

bool transport_equation::check_monotony(double check_value, double present_value, double previous_value) {
    bool ans;
    if ((std::min(present_value, previous_value) <= check_value) && (std::max(present_value, previous_value) >= check_value)) {
        ans = true;
    }
    else {
        ans = false;
    }
    return ans;
}

void transport_equation::change_layer() {
    for (unsigned i = 0; i < count_of_cells; ++i) {
        previous_layer[i] = present_layer[i];
        present_layer[i] = next_layer[i];
    }
}

void transport_equation::calculate_process(double time, scheme_coeffs c) {
    unsigned count_of_steps = time / t_step;

    if (save_txt_layers == true) {
        std::ofstream outfile;
        outfile.open(save_file_name, std::ios_base::app);

        //save initial conditions
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            outfile << previous_layer[i] << ", ";
        }
        outfile << previous_layer[count_of_cells - 1] << std::endl;
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            outfile << present_layer[i] << ", ";
        }
        outfile << present_layer[count_of_cells - 1] << std::endl;

        //time cicle 
        for (unsigned time = 1; time < count_of_steps; ++time) {
            next_layer[0] = left_boundary_condition(time * t_step);
            next_layer[1] = left_boundary_condition(time * t_step);
            next_layer[count_of_cells - 1] = right_boundary_condition(time * t_step);
            outfile << next_layer[0] << ", " << next_layer[1] << ", ";
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                next_layer[coord] = calc_new_value(c, coord);
                outfile << next_layer[coord] << ", ";
            }
            if (time < count_of_steps - 1) {
                outfile << next_layer[count_of_cells - 1] << std::endl;
            }
            else {
                outfile << next_layer[count_of_cells - 1];
            }
            change_layer();
        }
        outfile.close();
    }
    else {
        for (unsigned time = 1; time < count_of_steps; ++time) {
            next_layer[0] = left_boundary_condition(time * t_step);
            next_layer[1] = left_boundary_condition(time * t_step);
            next_layer[count_of_cells - 1] = right_boundary_condition(time * t_step);
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                next_layer[coord] = calc_new_value(c, coord);
            }
            change_layer();
        }
    }
    print_results();
}

void transport_equation::calculate_process_hybrid(double time) {
    if (entered_coeffs == false) {
        throw;
    }
    unsigned count_of_steps = time / t_step;
    double variable;

    if (save_txt_layers == true) {
        std::ofstream outfile;
        outfile.open(save_file_name, std::ios_base::app);

        //save initial conditions
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            outfile << previous_layer[i] << ", ";
        }
        outfile << previous_layer[count_of_cells - 1] << std::endl;
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            outfile << present_layer[i] << ", ";
        }
        outfile << present_layer[count_of_cells - 1] << std::endl;

        //time cicle 
        for (unsigned time = 1; time < count_of_steps; ++time) {
            next_layer[0] = left_boundary_condition(time * t_step);
            next_layer[1] = left_boundary_condition(time * t_step);
            next_layer[count_of_cells - 1] = right_boundary_condition(time * t_step);
            outfile << next_layer[0] << ", " << next_layer[1] << ", ";
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                for (unsigned k = 0; k < count_of_schemes; ++k) {
                    variable = calc_new_value(coeffs[k], coord);
                    if (check_monotony(variable, present_layer[coord], previous_layer[coord-1]) == true) {
                        k = count_of_schemes;
                    }
                }
                next_layer[coord] = variable;
                outfile << variable << ", ";
            }
            if (time < count_of_steps - 1) {
                outfile << next_layer[count_of_cells - 1] << std::endl;
            }
            else {
                outfile << next_layer[count_of_cells - 1];
            }
            change_layer();
        }
        outfile.close();
    }
    
    else {
        for (unsigned time = 1; time < count_of_steps; ++time) {
            next_layer[0] = left_boundary_condition(time * t_step);
            next_layer[1] = left_boundary_condition(time * t_step);
            next_layer[count_of_cells - 1] = right_boundary_condition(time * t_step);
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                for (unsigned k = 0; k < count_of_schemes; ++k) {
                    variable = calc_new_value(coeffs[k], coord);
                    if (check_monotony(variable, present_layer[coord], previous_layer[coord-1]) == true) {
                        k = count_of_schemes;
                    }
                }
                next_layer[coord] = variable;
            }
            change_layer();
        }
    }
    print_results();
}

void transport_equation::print_results() {
    for (unsigned i = 0; i < count_of_cells; ++i) {
        std::cout << present_layer[i] << " ";
    }
    std::cout << std::endl;
}

void transport_equation::enter_file_name_for_saving(std::string str) {
    save_file_name = str;
}

void analitical_solution_save(double t, double l_bound, unsigned count_of_cells, double coord_step) {
    std::ofstream outfile2;
    outfile2.open("transport_equation_analitical_solution.txt", std::ios_base::app);
    for (unsigned i = 0; i < count_of_cells - 1; ++i) {
        outfile2 << analytical_solution(l_bound + i*coord_step, t) << ", ";
    }
    outfile2 << analytical_solution(l_bound + (count_of_cells - 1)*coord_step, t);
    outfile2.close();
}

int main()
{
    try {
        unsigned length = (right_bound - left_bound) / h + 1;
        double* init_1 = new double[length];
        for (unsigned i = 0; i < length; ++i) {
            init_1[i] = initial_conditions(left_bound + i*h);
        }
        double* init_2 = new double[length];
        for (unsigned i = 0; i < length; ++i) {
            init_2[i] = analytical_solution(left_bound + i*h, tau);
        }

        scheme_coeffs point_O, point_X, point_Y, point_Z, point_H, point_T;
        point_O.A = 0.0;                //first approximation order
        point_O.B = cfl;
        point_O.C = 1.0 - cfl;
        point_O.D = 0.0;
        std::string point_O_name = "transport_eq_point_O.txt";

        point_X.A = 0.0;                //first approximation order
        point_X.B = 7.0 / 11.0;
        point_X.C = 0.0;
        point_X.D = 4.0 / 11.0;
        std::string point_X_name = "transport_eq_point_X.txt";

        point_Y.A = 0.1;                //first approximation order
        point_Y.B = 0.0;
        point_Y.C = 0.9;
        point_Y.D = 0.0;
        std::string point_Y_name = "transport_eq_point_Y.txt";

        point_Z.A = 7.0 / 16.0;        //first approximation order
        point_Z.B = 0.0;
        point_Z.C = 0.0;
        point_Z.D = 9.0 / 16.0;
        std::string point_Z_name = "transport_eq_point_Z.txt";

        point_H.A = - 25.0 / 857.0;     //second approximation order
        point_H.B = 909.0 / 4285.0;
        point_H.C = 3666.0 / 4285.0;
        point_H.D = - 33.0 / 857.0;
        std::string point_H_name = "transport_eq_point_H_2_order.txt";

        point_T.A = -0.035;             //third approximtion order
        point_T.B = 63.0/275.0;
        point_T.C = 0.84;
        point_T.D = -3.0/88.0;
        std::string point_T_name = "transport_eq_point_T_3_order.txt";

        scheme_coeffs* approximation_c = new scheme_coeffs[2];
        approximation_c[0] = point_H;
        approximation_c[1] = point_O;
        std::string hybrid_scheme_name = "transport_eq_hybrid_scheme.txt";

        //analitical_solution_save(100.0*tau, left_bound, length, h);

        transport_equation solution_hybrid(left_bound, right_bound, length, tau, true, 2);
        solution_hybrid.set_initial_layers(init_1, init_2);
        solution_hybrid.set_scheme_coeffs(approximation_c);
        solution_hybrid.enter_file_name_for_saving(hybrid_scheme_name);
        solution_hybrid.calculate_process_hybrid(100.0*tau);
        
        transport_equation solution_1_1(left_bound, right_bound, length, tau, true, 1);
        solution_1_1.set_initial_layers(init_1, init_2);
        solution_1_1.enter_file_name_for_saving(point_O_name);
        solution_1_1.calculate_process(100.0*tau, point_O);

        transport_equation solution_1_2(left_bound, right_bound, length, tau, true, 1);
        solution_1_2.set_initial_layers(init_1, init_2);
        solution_1_2.enter_file_name_for_saving(point_X_name);
        solution_1_2.calculate_process(100.0*tau, point_X);

        transport_equation solution_1_3(left_bound, right_bound, length, tau, true, 1);
        solution_1_3.set_initial_layers(init_1, init_2);
        solution_1_3.enter_file_name_for_saving(point_Y_name);
        solution_1_3.calculate_process(100.0*tau, point_Y);

        transport_equation solution_1_4(left_bound, right_bound, length, tau, true, 1);
        solution_1_4.set_initial_layers(init_1, init_2);
        solution_1_4.enter_file_name_for_saving(point_Z_name);
        solution_1_4.calculate_process(100.0*tau, point_Z);

        transport_equation solution_2(left_bound, right_bound, length, tau, true, 1);
        solution_2.set_initial_layers(init_1, init_2);
        solution_2.enter_file_name_for_saving(point_H_name);
        solution_2.calculate_process(100.0*tau, point_H);

        transport_equation solution_3(left_bound, right_bound, length, tau, true, 1);
        solution_3.set_initial_layers(init_1, init_2);
        solution_3.enter_file_name_for_saving(point_T_name);
        solution_3.calculate_process(100.0*tau, point_T);
    }
    
    catch (...) {
        std::cout << "firstly, entered coefficients for any schemes with 'set_scheme_coeffs' method" << std::endl;
    }
return 0;
}
