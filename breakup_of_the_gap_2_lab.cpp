// solving the Riemann problem on the decay of the gap by reducing the problem to 2 transfer equations
// u'_t + 1/rho * p'_x = 0
// p'_t + rho * c0^2 * u'_x = 0
// convert to:
// y'_t + c0 * y'_x = 0
// z'_t - c0 * z'_x = 0
//
// direct scheme:
// n+1: -   -   - * -   -   -
// n  : - * - * - * -   -   -
// n-1: -   -   -   - * -   -
//       m-2 m-1  m  m+1 m+2
// U(n+1, m) = A*U(n, m-2) + B*U(n, m-1) + C*U(n, m) + D*U(n-1, m+1)
//
// reverse scheme:
// n+1: -   -   - * -   -   -
// n  : -   -   - * - * - * -
// n-1: -   - * -   -   -   -
//       m-2 m-1  m  m+1 m+2
// U(n+1, m) = A*U(n, m+2) + B*U(n, m+1) + C*U(n, m) + D*U(n-1, m-1)

#include <iostream>
#include <fstream>
#include <math.h>
#include <exception>

const double cfl = 0.2;                 // Current number
const double h = 0.01;                  // coordinate step
const double left_bound = -1.0;
const double right_bound = 1.0;
const double rho_0 = 0.25;
const double c_0 = 2.0;
const double u_left = 1.0;
const double u_right = 0.0;
const double p_left = 5.0;
const double p_right = 2.0;
const double start_bound_of_gap = 0.0;
const double tau = cfl * h / c_0;    // time step
std::string default_file_name ("breakup_of_the_gap");

struct scheme_coeffs {
    double A;
    double B;
    double C;
    double D;
};

struct conservative_variables {
    double u_plus;
    double u_minus;
};

struct start_variables {
    double u;
    double p;
};

start_variables initial_conditions (double x) {
    start_variables ans;
    if (x < start_bound_of_gap)
    {
        ans.p = p_left;
        ans.u = u_left;
    }
    else
    {
        ans.p = p_right;
        ans.u = u_right;
    }
    return ans;
}

conservative_variables start_var_to_conservative_var (start_variables var) {
    conservative_variables ans;
    ans.u_plus = var.u + var.p/rho_0/c_0;
    ans.u_minus = var.u - var.p/rho_0/c_0;
    return ans;
}

start_variables conservative_var_to_start_var (double u_plus, double u_minus) {
    start_variables ans;
    ans.u = 0.5 * (u_plus + u_minus);
    ans.p = 0.5 * (u_plus - u_minus) * rho_0 * c_0;
    return ans;
}

start_variables analytical_solution (double x, double t) {
    start_variables ans;
    if (x <= start_bound_of_gap - c_0*t) 
    {
        ans = initial_conditions(-0.5);
    }
    else if (x >= start_bound_of_gap + c_0*t)
    {
        ans = initial_conditions(0.5);
    }
    else
    {
        ans.u = 0.5 * (u_left + u_right) + 0.5 * (p_left - p_right) / rho_0 / c_0;
        ans.p = (u_left - u_right) * 0.5 * rho_0 * c_0 + 0.5 * (p_left + p_right);
    }
    return ans;
}

conservative_variables left_boundary_condition (double t) {
    return start_var_to_conservative_var(initial_conditions(start_bound_of_gap-0.5));     // p(-L, t) = pL, u(-L, t) = uL
}

conservative_variables right_boundary_condition (double t) {
    return start_var_to_conservative_var(initial_conditions(start_bound_of_gap + 0.5));      // p(L, t) = pR, u(L, t) = uR
}

class breakup_of_the_gap {

    public:
        breakup_of_the_gap(double left_boundary, double right_boundary, unsigned number_of_cells, double time_step, bool save_txt_file, unsigned number_of_schemes);
        ~breakup_of_the_gap();
        void set_initial_layers(double* first_layer_u_plus, double* second_layer_u_plus, double* first_layer_u_minus, double* second_layer_u_minus);
        void calculate_process(double time, scheme_coeffs c);       //calculate all process with fixed coefficients from structure c
        void print_results();

        void set_scheme_coeffs(scheme_coeffs* c);
        void calculate_process_hybrid(double time);                 //calculate process with a hybrid scheme - combine some schemes, which takes from previous method
        void enter_file_name_for_saving(std::string str_p, std::string str_u);

    private:
        void change_layer();
        double calc_new_value_plus(scheme_coeffs c, unsigned ind);
        double calc_new_value_minus(scheme_coeffs c, unsigned ind);
        bool check_monotony(double check_value, double present_value, double previous_value);

        bool save_txt_layers;
        unsigned count_of_schemes;
        unsigned count_of_cells;
        double l_bound;
        double r_bound;
        double t_step;
        scheme_coeffs* coeffs;
        bool entered_coeffs = false;
        std::string save_file_name_p = default_file_name + "_p.txt";
        std::string save_file_name_u = default_file_name + "_u.txt";

        // arrays with state in 3 times level for 2 variables
        double* previous_layer_u_plus;
        double* present_layer_u_plus;
        double* next_layer_u_plus;
        double* previous_layer_u_minus;
        double* present_layer_u_minus;
        double* next_layer_u_minus;
};

breakup_of_the_gap::breakup_of_the_gap(double left_boundary, double right_boundary, unsigned number_of_cells, double time_step, bool save_txt_file, unsigned number_of_schemes) {
    l_bound = left_boundary;
    r_bound = right_boundary;
    t_step = time_step;
    save_txt_layers = save_txt_file;
    count_of_schemes = number_of_schemes;
    count_of_cells = number_of_cells;

    previous_layer_u_plus = new double[count_of_cells];
    present_layer_u_plus = new double[count_of_cells];
    next_layer_u_plus = new double[count_of_cells];
    previous_layer_u_minus = new double[count_of_cells];
    present_layer_u_minus = new double[count_of_cells];
    next_layer_u_minus = new double[count_of_cells];
    coeffs = new scheme_coeffs[number_of_schemes];
    std::cout << "class exemplar was created" << std::endl;
}

breakup_of_the_gap::~breakup_of_the_gap() {
    delete[] previous_layer_u_plus;
    delete[] present_layer_u_plus;
    delete[] next_layer_u_plus;
    delete[] previous_layer_u_minus;
    delete[] present_layer_u_minus;
    delete[] next_layer_u_minus;
    delete[] coeffs;
    std::cout << "class exemplar was deleted" << std::endl;
}

void breakup_of_the_gap::set_initial_layers(double* first_layer_u_plus, double* second_layer_u_plus, double* first_layer_u_minus, double* second_layer_u_minus) {
    for (unsigned i = 0; i < count_of_cells; ++i) {
        previous_layer_u_plus[i] = first_layer_u_plus[i];
        present_layer_u_plus[i] = second_layer_u_plus[i];
        previous_layer_u_minus[i] = first_layer_u_minus[i];
        present_layer_u_minus[i] = second_layer_u_minus[i];
    }
    std::cout << "first and second layers was inicialised" << std::endl;
}

void breakup_of_the_gap::set_scheme_coeffs(scheme_coeffs* c) {
    for (unsigned i = 0; i < count_of_schemes; ++i) {
        coeffs[i] = c[i];
    }
    entered_coeffs = true;
    std::cout << "scheme coefficients was inicialised, hybrid scheme unlocked" << std::endl;
}

double breakup_of_the_gap::calc_new_value_plus(scheme_coeffs c, unsigned ind) {
    return c.A*present_layer_u_plus[ind - 2] + c.B*present_layer_u_plus[ind - 1] + c.C*present_layer_u_plus[ind] + c.D*previous_layer_u_plus[ind + 1];
}

double breakup_of_the_gap::calc_new_value_minus(scheme_coeffs c, unsigned ind) {
    return c.A*present_layer_u_minus[ind + 2] + c.B*present_layer_u_minus[ind + 1] + c.C*present_layer_u_minus[ind] + c.D*previous_layer_u_minus[ind - 1];
}

bool breakup_of_the_gap::check_monotony(double check_value, double present_value, double previous_value) {
    bool ans;
    if ((std::min(present_value, previous_value) <= check_value) && (std::max(present_value, previous_value) >= check_value)) {
        ans = true;
    }
    else {
        ans = false;
    }
    return ans;
}

void breakup_of_the_gap::change_layer() {
    for (unsigned i = 0; i < count_of_cells; ++i) {
        previous_layer_u_plus[i] = present_layer_u_plus[i];
        present_layer_u_plus[i] = next_layer_u_plus[i];

        previous_layer_u_minus[i] = present_layer_u_minus[i];
        present_layer_u_minus[i] = next_layer_u_minus[i];
    }
}

void breakup_of_the_gap::calculate_process(double time, scheme_coeffs c) {
    unsigned count_of_steps = time / t_step;
    start_variables var;

    if (save_txt_layers == true) {
        std::ofstream outfile_u;
        outfile_u.open(save_file_name_u, std::ios_base::app);
        std::ofstream outfile_p;
        outfile_p.open(save_file_name_p, std::ios_base::app);

        //save initial conditions
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            var = conservative_var_to_start_var(previous_layer_u_plus[i], previous_layer_u_minus[i]);
            outfile_u << var.u << ", ";
            outfile_p << var.p << ", ";
        }
        var = conservative_var_to_start_var(previous_layer_u_plus[count_of_cells - 1], previous_layer_u_minus[count_of_cells - 1]);
        outfile_u << var.u << std::endl;
        outfile_p << var.p << std::endl;
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            var = conservative_var_to_start_var(present_layer_u_plus[i], present_layer_u_minus[i]);
            outfile_u << var.u << ", ";
            outfile_p << var.p << ", ";
        }
        var = conservative_var_to_start_var(present_layer_u_plus[count_of_cells - 1], present_layer_u_minus[count_of_cells - 1]);
        outfile_u << var.u << std::endl;
        outfile_p << var.p << std::endl;

        //time cicle 
        for (unsigned time = 1; time < count_of_steps; ++time) {
            // calculate new layer, set boundary condition
            next_layer_u_plus[0] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[1] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_plus;

            next_layer_u_minus[0] = left_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 2] = right_boundary_condition(time * t_step).u_minus;

            // cicles for coordinate
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                next_layer_u_plus[coord] = calc_new_value_plus(c, coord);
            }
            for (unsigned coord = count_of_cells - 3; coord > 0; coord--) {
                next_layer_u_minus[coord] = calc_new_value_minus(c, coord);
            }
            
            // save data in files
            for (unsigned coord = 0; coord < count_of_cells-1; ++coord) {
                var = conservative_var_to_start_var(next_layer_u_plus[coord], next_layer_u_minus[coord]);
                outfile_u << var.u << ", ";
                outfile_p << var.p << ", ";
            }
            if (time < count_of_steps - 1) {
                var = conservative_var_to_start_var(next_layer_u_plus[count_of_cells - 1], next_layer_u_minus[count_of_cells - 1]);
                outfile_u << var.u << std::endl;
                outfile_p << var.p << std::endl;
            }
            else {
                var = conservative_var_to_start_var(next_layer_u_plus[count_of_cells - 1], next_layer_u_minus[count_of_cells - 1]);
                outfile_u << var.u;
                outfile_p << var.p;
            }
            change_layer();
        }
        outfile_u.close();
        outfile_p.close();
    }
    else {
        for (unsigned time = 1; time < count_of_steps; ++time) {
            // calculate new layer, set boundary condition
            next_layer_u_plus[0] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[1] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_plus;

            next_layer_u_minus[0] = left_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 2] = right_boundary_condition(time * t_step).u_minus;

            // cicles for coordinate
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                next_layer_u_plus[coord] = calc_new_value_plus(c, coord);
            }
            for (unsigned coord = count_of_cells - 3; coord > 0; coord--) {
                next_layer_u_minus[coord] = calc_new_value_minus(c, coord);
            }
            change_layer();
        }
    }
    print_results();
}

void breakup_of_the_gap::calculate_process_hybrid(double time) {
    if (entered_coeffs == false) {
        throw;
    }
    unsigned count_of_steps = time / t_step;
    start_variables var;
    double temporary_variable;

    if (save_txt_layers == true) {
        std::ofstream outfile_u;
        outfile_u.open(save_file_name_u, std::ios_base::app);
        std::ofstream outfile_p;
        outfile_p.open(save_file_name_p, std::ios_base::app);

        //save initial conditions
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            var = conservative_var_to_start_var(previous_layer_u_plus[i], previous_layer_u_minus[i]);
            outfile_u << var.u << ", ";
            outfile_p << var.p << ", ";
        }
        var = conservative_var_to_start_var(previous_layer_u_plus[count_of_cells - 1], previous_layer_u_minus[count_of_cells - 1]);
        outfile_u << var.u << std::endl;
        outfile_p << var.p << std::endl;
        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            var = conservative_var_to_start_var(present_layer_u_plus[i], present_layer_u_minus[i]);
            outfile_u << var.u << ", ";
            outfile_p << var.p << ", ";
        }
        var = conservative_var_to_start_var(present_layer_u_plus[count_of_cells - 1], present_layer_u_minus[count_of_cells - 1]);
        outfile_u << var.u << std::endl;
        outfile_p << var.p << std::endl;

        //time cicle 
        for (unsigned time = 1; time < count_of_steps; ++time) {
            // calculate new layer, set boundary condition
            next_layer_u_plus[0] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[1] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_plus;

            next_layer_u_minus[0] = left_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 2] = right_boundary_condition(time * t_step).u_minus;

            // cicles for coordinate
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                for (unsigned k = 0; k < count_of_schemes; ++k) {
                    temporary_variable = calc_new_value_plus(coeffs[k], coord);
                    if (check_monotony(temporary_variable, present_layer_u_plus[coord], previous_layer_u_plus[coord-1]) == true) {
                        k = count_of_schemes;
                    }
                }
                next_layer_u_plus[coord] = temporary_variable;
            }
            for (unsigned coord = count_of_cells - 3; coord > 0; coord--) {
                for (unsigned k = 0; k < count_of_schemes; ++k) {
                    temporary_variable = calc_new_value_minus(coeffs[k], coord);
                    if (check_monotony(temporary_variable, present_layer_u_minus[coord], previous_layer_u_minus[coord+1]) == true) {
                        k = count_of_schemes;
                    }
                }
                next_layer_u_minus[coord] = temporary_variable;
            }

            // save data in files
            for (unsigned coord = 0; coord < count_of_cells-1; ++coord) {
                var = conservative_var_to_start_var(next_layer_u_plus[coord], next_layer_u_minus[coord]);
                outfile_u << var.u << ", ";
                outfile_p << var.p << ", ";
            }
            if (time < count_of_steps - 1) {
                var = conservative_var_to_start_var(next_layer_u_plus[count_of_cells - 1], next_layer_u_minus[count_of_cells - 1]);
                outfile_u << var.u << std::endl;
                outfile_p << var.p << std::endl;
            }
            else {
                var = conservative_var_to_start_var(next_layer_u_plus[count_of_cells - 1], next_layer_u_minus[count_of_cells - 1]);
                outfile_u << var.u;
                outfile_p << var.p;
            }
            change_layer();
        }
        outfile_u.close();
        outfile_p.close();
    }
    
    else {
        for (unsigned time = 1; time < count_of_steps; ++time) {
            // calculate new layer, set boundary condition
            next_layer_u_plus[0] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[1] = left_boundary_condition(time * t_step).u_plus;
            next_layer_u_plus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_plus;

            next_layer_u_minus[0] = left_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 1] = right_boundary_condition(time * t_step).u_minus;
            next_layer_u_minus[count_of_cells - 2] = right_boundary_condition(time * t_step).u_minus;

            // cicles for coordinate
            for (unsigned coord = 2; coord < count_of_cells-1; ++coord) {
                for (unsigned k = 0; k < count_of_schemes; ++k) {
                    temporary_variable = calc_new_value_plus(coeffs[k], coord);
                    if (check_monotony(temporary_variable, present_layer_u_plus[coord], previous_layer_u_plus[coord-1]) == true) {
                        k = count_of_schemes;
                    }
                }
                next_layer_u_plus[coord] = temporary_variable;
            }
            for (unsigned coord = count_of_cells - 3; coord > 0; coord--) {
                for (unsigned k = 0; k < count_of_schemes; ++k) {
                    temporary_variable = calc_new_value_minus(coeffs[k], coord);
                    if (check_monotony(temporary_variable, present_layer_u_minus[coord], previous_layer_u_minus[coord+1]) == true) {
                        k = count_of_schemes;
                    }
                }
                next_layer_u_minus[coord] = temporary_variable;
            }
            change_layer();
        }
    }
    print_results();
}

void breakup_of_the_gap::print_results() {
    start_variables var;
    std::cout << "u in layer: ";
    for (unsigned i = 0; i < count_of_cells; ++i) {
        var = conservative_var_to_start_var(present_layer_u_plus[i], present_layer_u_minus[i]);
        std::cout << var.u << " ";
    }
    std::cout << std::endl;
    std::cout << "p in layer: ";
    for (unsigned i = 0; i < count_of_cells; ++i) {
        var = conservative_var_to_start_var(present_layer_u_plus[i], present_layer_u_minus[i]);
        std::cout << var.p << " ";
    }
    std::cout << std::endl;
}

void breakup_of_the_gap::enter_file_name_for_saving(std::string str_p, std::string str_u) {
    save_file_name_p = str_p;
    save_file_name_u = str_u;
}

void analitical_solution_save(double t, double l_bound, unsigned count_of_cells, double coord_step) {
    std::ofstream outfile_u;
    outfile_u.open("breakup_of_the_gap_analitical_solution_u.txt", std::ios_base::app);
    std::ofstream outfile_p;
    outfile_p.open("breakup_of_the_gap_analitical_solution_p.txt", std::ios_base::app);
    
    start_variables var;
    for (unsigned i = 0; i < count_of_cells - 1; ++i) {
        var = analytical_solution(l_bound + i*coord_step, t);
        outfile_u << var.u << ", ";
        outfile_p << var.p << ", ";
    }
    var = analytical_solution(l_bound + (count_of_cells - 1)*coord_step, t);
    outfile_u << var.u;
    outfile_p << var.p;
    outfile_u.close();
    outfile_p.close();
}

int main()
{
    start_variables var;
    try {
        unsigned length = (right_bound - left_bound) / h + 1;
        double* init_1_u_plus = new double[length];
        double* init_1_u_minus = new double[length];
        double* init_2_u_plus = new double[length];
        double* init_2_u_minus = new double[length];

        for (unsigned i = 0; i < length; ++i) {
            var = initial_conditions(left_bound + i*h);
            init_1_u_plus[i] = start_var_to_conservative_var(var).u_plus;
            init_1_u_minus[i] = start_var_to_conservative_var(var).u_minus;
        }
        for (unsigned i = 0; i < length; ++i) {
            var = analytical_solution(left_bound + i*h, tau);
            init_2_u_plus[i] = start_var_to_conservative_var(var).u_plus;
            init_2_u_minus[i] = start_var_to_conservative_var(var).u_minus;
        }

        // setting schemes coeeficients
        scheme_coeffs point_O;
        point_O.A = 0.0;                //first approximation order
        point_O.B = cfl;
        point_O.C = 1.0 - cfl;
        point_O.D = 0.0;

        scheme_coeffs point_H;
        point_H.A = - 25.0 / 857.0;     //second approximation order
        point_H.B = 909.0 / 4285.0;
        point_H.C = 3666.0 / 4285.0;
        point_H.D = - 33.0 / 857.0;

        scheme_coeffs point_T;
        point_T.A = -0.035;             //third approximtion order
        point_T.B = 63.0/275.0;
        point_T.C = 0.84;
        point_T.D = -3.0/88.0;

        scheme_coeffs* approximation_c = new scheme_coeffs[3];
        approximation_c[0] = point_T;
        approximation_c[1] = point_H;
        approximation_c[2] = point_O;

        breakup_of_the_gap solution_1(left_bound, right_bound, length, tau, true, 1);
        solution_1.set_initial_layers(init_1_u_plus, init_2_u_plus, init_1_u_minus, init_2_u_minus);
        solution_1.enter_file_name_for_saving("breakup_of_the_gap_point_O_p.txt", "breakup_of_the_gap_point_O_u.txt");
        solution_1.calculate_process(100.0*tau, point_O);

        breakup_of_the_gap solution_2(left_bound, right_bound, length, tau, true, 1);
        solution_2.set_initial_layers(init_1_u_plus, init_2_u_plus, init_1_u_minus, init_2_u_minus);
        solution_2.enter_file_name_for_saving("breakup_of_the_gap_point_H_p.txt", "breakup_of_the_gap_point_H_u.txt");
        solution_2.calculate_process(100.0*tau, point_H);

        breakup_of_the_gap solution_hybrid(left_bound, right_bound, length, tau, true, 3);
        solution_hybrid.set_initial_layers(init_1_u_plus, init_2_u_plus, init_1_u_minus, init_2_u_minus);
        solution_hybrid.enter_file_name_for_saving("breakup_of_the_gap_hybrid_p.txt", "breakup_of_the_gap_hybrid_u.txt");
        solution_hybrid.set_scheme_coeffs(approximation_c);
        solution_hybrid.calculate_process_hybrid(100.0*tau);

        // analitical solution
        analitical_solution_save(100.0*tau, left_bound, length, h);
    }
    
    catch (...) {
        std::cout << "firstly, entered coefficients for any schemes with 'set_scheme_coeffs' method" << std::endl;
    }
return 0;
}
