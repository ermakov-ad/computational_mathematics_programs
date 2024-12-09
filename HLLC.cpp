#include <iostream>
#include <math.h>
#include <fstream>

const double left_bound = 0.0;
const double right_bound = 1.0;
const unsigned cells_count = 100;
const double h = (right_bound - left_bound) / cells_count;

// first test:
// const double start_bound_of_gap = 0.3;
// const double rhoL = 1.0;
// const double rhoR = 0.125;
// const double uL = 0.75;
// const double uR = 0.0;
// const double pL = 1.0;
// const double pR = 0.1;
// const double time_period = 0.2;
// const double save_time = 0.01;

// second test:
// const double start_bound_of_gap = 0.4;
// const double rhoL = 5.99924;
// const double rhoR = 5.99242;
// const double uL = 19.5975;
// const double uR = -6.19633;
// const double pL = 460.894;
// const double pR = 46.0950;
// const double time_period = 0.035;
// const double save_time = 0.001;

// break with reaction:
const double start_bound_of_gap = 0.4;
const double rhoL = 2.0;
const double rhoR = 1.5;
const double uL = 200.0;
const double uR = -180.0;
const double pL = 8000.0;
const double pR = 10000.0;
const double time_period = 0.005;
const double save_time = 0.001;

const double Courant_number = 0.9;

// T = 1200
const double R = 8.314;
const double gamma_ = 1.4;

struct start_variables {
    double rho;
    double u;
    double p;
};

struct conservative_variables {
    double rho;
    double rho_u;
    double E;
};

struct flow_variables {
    double rho_u;
    double rho_u_u_p;
    double u_E_p;
};

conservative_variables operator+ (conservative_variables u1, conservative_variables u2) {
    conservative_variables ans;
    ans.rho = u1.rho + u2.rho;
    ans.rho_u = u1.rho_u + u2.rho_u;
    ans.E = u1.E + u2.E;
    return ans;
}

conservative_variables operator- (conservative_variables u1, conservative_variables u2) {
    conservative_variables ans;
    ans.rho = u1.rho - u2.rho;
    ans.rho_u = u1.rho_u - u2.rho_u;
    ans.E = u1.E - u2.E;
    return ans;
}

flow_variables operator+ (flow_variables f1, flow_variables f2) {
    flow_variables ans;
    ans.rho_u = f1.rho_u + f2.rho_u;
    ans.rho_u_u_p = f1.rho_u_u_p + f2.rho_u_u_p;
    ans.u_E_p = f1.u_E_p + f2.u_E_p;
    return ans;
}

flow_variables operator- (flow_variables f1, flow_variables f2) {
    flow_variables ans;
    ans.rho_u = f1.rho_u - f2.rho_u;
    ans.rho_u_u_p = f1.rho_u_u_p - f2.rho_u_u_p;
    ans.u_E_p = f1.u_E_p - f2.u_E_p;
    return ans;
}

// E = rho * (u^2/2 + e)
// e = E/rho - u^2/2 = (E - 0.5*(rho*u)^2/rho) / rho
double internal_energy(conservative_variables state) {
    return (state.E - 0.5*state.rho_u*state.rho_u/state.rho) / state.rho;
}

// e = p / rho / (gamma - 1)
double internal_energy(start_variables state) {
    return state.p / state.rho / (gamma_ - 1.0);
}

// e = Cv * T; Cv = R / (gamma - 1)
// T = Cv / e = R / (gamma - 1) / e
double calculate_temperature(conservative_variables state) {
    double e = internal_energy(state);
    return e * (gamma_ - 1.0) / R; 
}

// e = E/rho - u^2/2; e = p / rho / (gamma - 1) =>
// p = e * rho * (gamma - 1)
double pressure_from_conservative(conservative_variables state) {
    double e = internal_energy(state);
    return e * state.rho * (gamma_ - 1.0);
}

conservative_variables start_to_conservative(start_variables state) {
    double e = internal_energy(state);
    conservative_variables ans;
    ans.rho = state.rho;
    ans.rho_u = state.rho * state.u;
    ans.E = state.rho * (0.5 * state.u * state.u + e);
    return ans;
}

start_variables conservative_to_start(conservative_variables state) {
    start_variables ans;
    ans.rho = state.rho;
    ans.u = state.rho_u / state.rho;
    ans.p = pressure_from_conservative(state);
    return ans;
}

double find_sound_speed(conservative_variables state) {
    double p = pressure_from_conservative(state);
    return sqrt(gamma_ * p / state.rho);
}

flow_variables find_flow(conservative_variables state) {
    double p = pressure_from_conservative(state);
    flow_variables ans;
    ans.rho_u = state.rho_u;
    ans.rho_u_u_p = state.rho_u * state.rho_u / state.rho + p;
    ans.u_E_p = state.rho_u * (state.E + p) / state.rho;
    return ans;
}

flow_variables conservative_to_flow(double wave_speed, conservative_variables U) {
    flow_variables ans;
    ans.rho_u = wave_speed * U.rho;
    ans.rho_u_u_p = wave_speed * U.rho_u;
    ans.u_E_p = wave_speed * U.E;
    return ans;
}

conservative_variables flow_to_conservative(double reverse_wave_speed, flow_variables F) {
    conservative_variables ans;
    ans.rho = reverse_wave_speed * F.rho_u;
    ans.rho_u = reverse_wave_speed * F.rho_u_u_p;
    ans.E = reverse_wave_speed * F.u_E_p;
    return ans;
}

conservative_variables left_boundary_condition(double t) {
    start_variables left_state;
    left_state.p = pL;
    left_state.rho = rhoL;
    left_state.u = uL;
    conservative_variables ans = start_to_conservative(left_state);
    return ans;
}

conservative_variables right_boundary_condition(double t) {
    start_variables right_state;
    right_state.p = pR;
    right_state.rho = rhoR;
    right_state.u = uR;
    conservative_variables ans = start_to_conservative(right_state);
    return ans;
}

start_variables initial_conditions(double coord) {
    start_variables ans;
    if (coord < start_bound_of_gap) {
        ans.rho = rhoL;
        ans.u = uL;
        ans.p = pL;
    }
    else {
        ans.rho = rhoR;
        ans.u = uR;
        ans.p = pR;
    }
    return ans;
}

class HLLC_solver {

    public:

        HLLC_solver(unsigned cells_count, double coordinate_step);
        ~HLLC_solver();

        void set_initial_distribution(start_variables *distribution);
        void calculate_process(double process_time, double write_time, double CFL);

    private:

        double find_time_step();
        flow_variables find_flow_between_cells(unsigned index_left);        //F(i+1/2) = flow between U[i] and U[i+1]
        void save_present_layer();

        unsigned count_of_cells;
        double h;

        double end_time;
        double calculation_time = 0.0;
        double write_time_step;
        double previous_save_time = 0.0;
        double cfl;

        conservative_variables* present_layer;
        conservative_variables* next_layer;
        flow_variables* present_flow;
};

HLLC_solver::HLLC_solver(unsigned N_cells, double coordinate_step) {
    count_of_cells = N_cells;
    h = coordinate_step;
    present_layer = new conservative_variables[N_cells];
    next_layer = new conservative_variables[N_cells];
    present_flow = new flow_variables[N_cells];
    std::cout << "class exemplar was created" << std::endl;
}

HLLC_solver::~HLLC_solver() {
    delete[] present_layer;
    delete[] next_layer;
    delete[] present_flow;
    std::cout << "class exemplar was deleted" << std::endl;
}

void HLLC_solver::set_initial_distribution(start_variables *distribution) {
    for (unsigned i = 0; i < count_of_cells; ++i) {
        present_layer[i] = start_to_conservative(distribution[i]);
    }
}

double HLLC_solver::find_time_step() {
    double tau = end_time;
    for (unsigned i = 0; i < count_of_cells; ++i) {
        double u = abs(present_layer[i].rho_u / present_layer[i].rho);
        double c = find_sound_speed(present_layer[i]);
        double dt = h / (u + c);
        if (dt < tau) {
            tau = dt;
        }
    }
    tau *= cfl;
    if (calculation_time + tau > previous_save_time + write_time_step) {
        tau = previous_save_time + write_time_step - calculation_time;
    } 
    return tau;
}

flow_variables HLLC_solver::find_flow_between_cells(unsigned index_left) {
    double a_l = find_sound_speed(present_layer[index_left]);
    double a_r = find_sound_speed(present_layer[index_left + 1]);
    double u_l = present_layer[index_left].rho_u / present_layer[index_left].rho;
    double u_r = present_layer[index_left + 1].rho_u / present_layer[index_left + 1].rho;
    double SL = std::min(u_l - a_l, u_r - a_r);
    double SR = std::max(u_l + a_l, u_r + a_r);
    double p_l = pressure_from_conservative(present_layer[index_left]);
    double p_r = pressure_from_conservative(present_layer[index_left + 1]);
    double rho_l = present_layer[index_left].rho;
    double rho_r = present_layer[index_left + 1].rho;
    double S_star = (p_r - p_l + rho_l * u_l * (SL - u_l) - rho_r * u_r * (SR - u_r)) / (rho_l * (SL - u_l) - rho_r * (SR - u_r));

    flow_variables ans;

    if (SL >= 0.0) {
        ans = find_flow(present_layer[index_left]);     //F(i+1/2) = FL
    }
    else if (S_star >= 0.0) {   //SL < 0.0
        //F(i+1/2) = F*L = FL + SL (U*L − UL)
        conservative_variables U_star_L;
        double multiplier = rho_l * (SL - u_l) / (SL - S_star);
        U_star_L.rho = multiplier;
        U_star_L.rho_u = multiplier * S_star;
        U_star_L.E = multiplier * (present_layer[index_left].E / rho_l + (S_star - u_l) * (S_star + p_l / (rho_l*(SL - u_l))));
        ans = find_flow(present_layer[index_left]) + conservative_to_flow(SL, U_star_L - present_layer[index_left]);
    }
    else if (SR <= 0.0) {
        ans = find_flow(present_layer[index_left + 1]);     //F(i+1/2) = FR
    }
    else {          //SR > 0.0, S* < 0.0
        //F(i+1/2) = F*R = FR + SR (U*R − UR)
        conservative_variables U_star_R;
        double multiplier = rho_r * (SR - u_r) / (SR - S_star);
        U_star_R.rho = multiplier;
        U_star_R.rho_u = multiplier * S_star;
        U_star_R.E = multiplier * (present_layer[index_left + 1].E / rho_r + (S_star - u_r) * (S_star + p_r / (rho_r*(SR - u_r))));
        ans = find_flow(present_layer[index_left + 1]) + conservative_to_flow(SR, U_star_R - present_layer[index_left + 1]);
    }
    return ans;
}

void HLLC_solver::save_present_layer() {
    std::ofstream outfile_rho;
    outfile_rho.open("rho_field.txt", std::ios_base::app);
    std::ofstream outfile_u;
    outfile_u.open("u_field.txt", std::ios_base::app);
    std::ofstream outfile_p;
    outfile_p.open("p_field.txt", std::ios_base::app);
    std::ofstream outfile_T;
    outfile_T.open("T_field.txt", std::ios_base::app);

    outfile_rho << "time = " << calculation_time << std::endl;
    outfile_u << "time = " << calculation_time << std::endl;
    outfile_p << "time = " << calculation_time << std::endl;
    outfile_T << "time = " << calculation_time << std::endl;

    start_variables var;
    for (unsigned i = 0; i < count_of_cells - 1; ++i) {
        var = conservative_to_start(present_layer[i]);
        outfile_rho << var.rho << " ";
        outfile_u << var.u << " ";
        outfile_p << var.p << " ";
        outfile_T << calculate_temperature(present_layer[i]) << " ";
    }
    var = conservative_to_start(present_layer[count_of_cells - 1]);
    outfile_rho << var.rho << std::endl;
    outfile_u << var.u << std::endl;
    outfile_p << var.p << std::endl;
    outfile_T << calculate_temperature(present_layer[count_of_cells - 1]) << std::endl;

    outfile_rho.close();
    outfile_u.close();
    outfile_p.close();
    outfile_T.close();
}

void HLLC_solver::calculate_process(double process_time, double write_time, double CFL) {
    end_time = process_time;
    write_time_step = write_time;
    cfl = CFL;

    save_present_layer();
    while (calculation_time < end_time) {
        double dt = find_time_step();
        // time step
        next_layer[0] = left_boundary_condition(calculation_time);
        next_layer[count_of_cells - 1] = right_boundary_condition(calculation_time);

        for (unsigned i = 0; i < count_of_cells - 1; ++i) {
            present_flow[i] = find_flow_between_cells(i);   //F[i + 1/2]
        }
        for (unsigned i = 1; i < count_of_cells - 1; ++i) {
            // U(n+1, i) = U(n, i) - dt/dx * (F[i+1/2] - F[i-1/2])
            next_layer[i] = present_layer[i] - flow_to_conservative(dt / h, present_flow[i] - present_flow[i - 1]);
        }
        calculation_time += dt;
        // changing layer
        for (unsigned i = 0; i < count_of_cells; ++i) {
            present_layer[i] = next_layer[i];
        }
        //save layer if needed
        if (previous_save_time + write_time_step == calculation_time) {
            save_present_layer();
            previous_save_time += write_time_step;
        }
    }
}

int main()
{
    start_variables* init = new start_variables[cells_count];
    for (unsigned i = 0; i < cells_count; ++i) {
        init[i] = initial_conditions(left_bound + i * h + 0.5 * h);
    }
    HLLC_solver task(cells_count, h);
    task.set_initial_distribution(init);
    task.calculate_process(time_period, save_time, Courant_number);
return 0;
}