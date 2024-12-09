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

// 2 CH4 = C2H2 + 3 H2
const double reaction_temperature = 1673.0;
const double reacting_constant = 150.0;
const double q1L = 1.0;
const double q2L = 0.0;
const double q3L = 0.0;
const double q1R = 1.0;
const double q2R = 0.0;
const double q3R = 0.0;

const double time_period = 0.005;
const double save_time = 0.0001;

const double Courant_number = 0.5;

// T = 1200
const double R = 8.314;
const double gamma1 = 4.0 / 3.0;
const double gamma2 = 4.0 / 3.0;
const double gamma3 = 1.4;

struct start_variables {
    double rho;
    double u;
    double p;
    double q1;      // CH4
    double q2;      // C2H2
    double q3;      // H2; q1 + q2 + q3 = 1
};

struct conservative_variables {
    double rho;
    double rho_u;
    double E;
    double rho_q1;
    double rho_q2;
    double rho_q3;
};

struct flow_variables {
    double rho_u;
    double rho_u_u_p;
    double u_E_p;
    double rho_u_q1;
    double rho_u_q2;
    double rho_u_q3;
};

conservative_variables operator+ (conservative_variables u1, conservative_variables u2) {
    conservative_variables ans;
    ans.rho = u1.rho + u2.rho;
    ans.rho_u = u1.rho_u + u2.rho_u;
    ans.E = u1.E + u2.E;
    ans.rho_q1 = u1.rho_q1 + u2.rho_q1;
    ans.rho_q2 = u1.rho_q2 + u2.rho_q2;
    ans.rho_q3 = u1.rho_q3 + u2.rho_q3;
    return ans;
}

conservative_variables operator- (conservative_variables u1, conservative_variables u2) {
    conservative_variables ans;
    ans.rho = u1.rho - u2.rho;
    ans.rho_u = u1.rho_u - u2.rho_u;
    ans.E = u1.E - u2.E;
    ans.rho_q1 = u1.rho_q1 - u2.rho_q1;
    ans.rho_q2 = u1.rho_q2 - u2.rho_q2;
    ans.rho_q3 = u1.rho_q3 - u2.rho_q3;
    return ans;
}

flow_variables operator+ (flow_variables f1, flow_variables f2) {
    flow_variables ans;
    ans.rho_u = f1.rho_u + f2.rho_u;
    ans.rho_u_u_p = f1.rho_u_u_p + f2.rho_u_u_p;
    ans.u_E_p = f1.u_E_p + f2.u_E_p;
    ans.rho_u_q1 = f1.rho_u_q1 + f2.rho_u_q1;
    ans.rho_u_q2 = f1.rho_u_q2 + f2.rho_u_q2;
    ans.rho_u_q3 = f1.rho_u_q3 + f2.rho_u_q3;
    return ans;
}

flow_variables operator- (flow_variables f1, flow_variables f2) {
    flow_variables ans;
    ans.rho_u = f1.rho_u - f2.rho_u;
    ans.rho_u_u_p = f1.rho_u_u_p - f2.rho_u_u_p;
    ans.u_E_p = f1.u_E_p - f2.u_E_p;
    ans.rho_u_q1 = f1.rho_u_q1 - f2.rho_u_q1;
    ans.rho_u_q2 = f1.rho_u_q2 - f2.rho_u_q2;
    ans.rho_u_q3 = f1.rho_u_q3 - f2.rho_u_q3;
    return ans;
}

// E = rho * (u^2/2 + e)
// e = E/rho - u^2/2 = (E - 0.5*(rho*u)^2/rho) / rho
double internal_energy(conservative_variables state) {
    return (state.E - 0.5*state.rho_u*state.rho_u/state.rho) / state.rho;
}

// e = E/rho - u^2/2; e = p / rho / (gamma - 1) =>
// p = e * rho * (gamma - 1)
double pressure_from_conservative(conservative_variables state) {
    double e = internal_energy(state);
    return e * state.rho * state.rho / (state.rho_q1 / (gamma1 - 1.0) + state.rho_q2 / (gamma2 - 1.0) + state.rho_q3 / (gamma3 - 1.0));
}

// e = p / rho / (gamma - 1)
double internal_energy(start_variables state) {
    return (state.q1 / (gamma1 - 1.0) + state.q2 / (gamma2 - 1.0) + state.q3 / (gamma3 - 1.0) ) * state.p / state.rho;
}

// e = Cv * T; Cv = R / (gamma - 1)
// T = Cv / e = R / (gamma - 1) / e
double calculate_temperature(conservative_variables state) {
    double e = internal_energy(state);
    return e * state.rho / R / (state.rho_q1 / (gamma1 - 1.0) + state.rho_q2 / (gamma2 - 1.0) + state.rho_q3 / (gamma3 - 1.0)); 
}

conservative_variables start_to_conservative(start_variables state) {
    double e = internal_energy(state);
    conservative_variables ans;
    ans.rho = state.rho;
    ans.rho_u = state.rho * state.u;
    ans.E = state.rho * (0.5 * state.u * state.u + e);
    ans.rho_q1 = state.rho * state.q1;
    ans.rho_q2 = state.rho * state.q2;
    ans.rho_q3 = state.rho * state.q3;
    return ans;
}

start_variables conservative_to_start(conservative_variables state) {
    start_variables ans;
    ans.rho = state.rho;
    ans.u = state.rho_u / state.rho;
    ans.p = pressure_from_conservative(state);
    ans.q1 = state.rho_q1 / state.rho;
    ans.q2 = state.rho_q2 / state.rho;
    ans.q3 = state.rho_q3 / state.rho;
    return ans;
}

double find_sound_speed(conservative_variables state) {
    double p = pressure_from_conservative(state);
    return std::max(std::max(sqrt(gamma1 * p * state.rho_q1 / state.rho / state.rho), sqrt(gamma2 * p * state.rho_q2 / state.rho / state.rho)), sqrt(gamma3 * p * state.rho_q3 / state.rho / state.rho));
}

flow_variables find_flow(conservative_variables state) {
    double p = pressure_from_conservative(state);
    flow_variables ans;
    ans.rho_u = state.rho_u;
    ans.rho_u_u_p = state.rho_u * state.rho_u / state.rho + p;
    ans.u_E_p = state.rho_u * (state.E + p) / state.rho;
    ans.rho_u_q1 = state.rho_u * state.rho_q1 / state.rho;
    ans.rho_u_q2 = state.rho_u * state.rho_q2 / state.rho;
    ans.rho_u_q3 = state.rho_u * state.rho_q3 / state.rho;
    return ans;
}

flow_variables conservative_to_flow(double wave_speed, conservative_variables U) {
    flow_variables ans;
    ans.rho_u = wave_speed * U.rho;
    ans.rho_u_u_p = wave_speed * U.rho_u;
    ans.u_E_p = wave_speed * U.E;
    ans.rho_u_q1 = wave_speed * U.rho_q1;
    ans.rho_u_q2 = wave_speed * U.rho_q2;
    ans.rho_u_q3 = wave_speed * U.rho_q3;
    return ans;
}

conservative_variables flow_to_conservative(double reverse_wave_speed, flow_variables F) {
    conservative_variables ans;
    ans.rho = reverse_wave_speed * F.rho_u;
    ans.rho_u = reverse_wave_speed * F.rho_u_u_p;
    ans.E = reverse_wave_speed * F.u_E_p;
    ans.rho_q1 = reverse_wave_speed * F.rho_u_q1;
    ans.rho_q2 = reverse_wave_speed * F.rho_u_q2;
    ans.rho_q3 = reverse_wave_speed * F.rho_u_q3;
    return ans;
}

conservative_variables left_boundary_condition(double t) {
    start_variables left_state;
    left_state.p = pL;
    left_state.rho = rhoL;
    left_state.u = uL;
    left_state.q1 = q1L;
    left_state.q2 = q2L;
    left_state.q3 = q3L;
    conservative_variables ans = start_to_conservative(left_state);
    return ans;
}

conservative_variables right_boundary_condition(double t) {
    start_variables right_state;
    right_state.p = pR;
    right_state.rho = rhoR;
    right_state.u = uR;
    right_state.q1 = q1R;
    right_state.q2 = q2R;
    right_state.q3 = q3R;
    conservative_variables ans = start_to_conservative(right_state);
    return ans;
}

start_variables initial_conditions(double coord) {
    start_variables ans;
    if (coord < start_bound_of_gap) {
        ans.rho = rhoL;
        ans.u = uL;
        ans.p = pL;
        ans.q1 = q1L;
        ans.q2 = q2L;
        ans.q3 = q3L;
    }
    else {
        ans.rho = rhoR;
        ans.u = uR;
        ans.p = pR;
        ans.q1 = q1R;
        ans.q2 = q2R;
        ans.q3 = q3R;
    }
    return ans;
}

double find_SL(conservative_variables UL, conservative_variables UR) {
    double a_l = find_sound_speed(UL);
    double a_r = find_sound_speed(UR);
    double u_l = UL.rho_u / UL.rho;
    double u_r = UR.rho_u / UR.rho;
    double p_l = pressure_from_conservative(UL);
    double p_r = pressure_from_conservative(UR);
    double rho_l = UL.rho;
    double rho_r = UR.rho;

    double p0 = std::max(0.0, 0.5*(p_l + p_r) - 0.125*(u_r - u_l)*(rho_l + rho_r)*(a_l + a_r));
    double ql = 1.0;
    double qr = 1.0;
    if (p0 > p_l) {
        ql = sqrt(1.0 + 0.5 * (gamma3 + 1.0) / gamma3 * (p0 / p_l - 1.0));
    }
    if (p0 > p_r) {
        qr = sqrt(1.0 + 0.5 * (gamma3 + 1.0) / gamma3 * (p0 / p_r - 1.0));
    }
    return std::min(u_l - a_l*ql, u_r - a_r*qr);
}

double find_SR(conservative_variables UL, conservative_variables UR) {
    double a_l = find_sound_speed(UL);
    double a_r = find_sound_speed(UR);
    double u_l = UL.rho_u / UL.rho;
    double u_r = UR.rho_u / UR.rho;
    double p_l = pressure_from_conservative(UL);
    double p_r = pressure_from_conservative(UR);
    double rho_l = UL.rho;
    double rho_r = UR.rho;

    double p0 = std::max(0.0, 0.5*(p_l + p_r) - 0.125*(u_r - u_l)*(rho_l + rho_r)*(a_l + a_r));
    double ql = 1.0;
    double qr = 1.0;
    if (p0 > p_l) {
        ql = sqrt(1.0 + 0.5 * (gamma3 + 1.0) / gamma3 * (p0 / p_l - 1.0));
    }
    if (p0 > p_r) {
        qr = sqrt(1.0 + 0.5 * (gamma3 + 1.0) / gamma3 * (p0 / p_r - 1.0));
    }
    return std::max(u_l + a_l*ql, u_r + a_r*qr);
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
        void conducting_a_reaction(double time_step);

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
    double u_l = present_layer[index_left].rho_u / present_layer[index_left].rho;
    double u_r = present_layer[index_left + 1].rho_u / present_layer[index_left + 1].rho;
    double p_l = pressure_from_conservative(present_layer[index_left]);
    double p_r = pressure_from_conservative(present_layer[index_left + 1]);
    double rho_l = present_layer[index_left].rho;
    double rho_r = present_layer[index_left + 1].rho;

    double SL = find_SL(present_layer[index_left], present_layer[index_left + 1]);
    double SR = find_SR(present_layer[index_left], present_layer[index_left + 1]);
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
        U_star_L.rho_q1 = multiplier * present_layer[index_left].rho_q1 / present_layer[index_left].rho;
        U_star_L.rho_q2 = multiplier * present_layer[index_left].rho_q2 / present_layer[index_left].rho;
        U_star_L.rho_q3 = multiplier * present_layer[index_left].rho_q3 / present_layer[index_left].rho;
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
        U_star_R.rho_q1 = multiplier * present_layer[index_left + 1].rho_q1 / present_layer[index_left + 1].rho;
        U_star_R.rho_q2 = multiplier * present_layer[index_left + 1].rho_q2 / present_layer[index_left + 1].rho;
        U_star_R.rho_q3 = multiplier * present_layer[index_left + 1].rho_q3 / present_layer[index_left + 1].rho;
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
    std::ofstream outfile_q1;
    outfile_q1.open("q1_field.txt", std::ios_base::app);
    std::ofstream outfile_q2;
    outfile_q2.open("q2_field.txt", std::ios_base::app);
    std::ofstream outfile_q3;
    outfile_q3.open("q3_field.txt", std::ios_base::app);

    // outfile_rho << "time = " << calculation_time << std::endl;
    // outfile_u << "time = " << calculation_time << std::endl;
    // outfile_p << "time = " << calculation_time << std::endl;
    // outfile_T << "time = " << calculation_time << std::endl;

    start_variables var;
    for (unsigned i = 0; i < count_of_cells - 1; ++i) {
        var = conservative_to_start(present_layer[i]);
        outfile_rho << var.rho << " ";
        outfile_u << var.u << " ";
        outfile_p << var.p << " ";
        outfile_T << calculate_temperature(present_layer[i]) << " ";
        outfile_q1 << present_layer[i].rho_q1 / present_layer[i].rho << " ";
        outfile_q2 << present_layer[i].rho_q2 / present_layer[i].rho << " ";
        outfile_q3 << present_layer[i].rho_q3 / present_layer[i].rho << " ";
    }
    var = conservative_to_start(present_layer[count_of_cells - 1]);
    outfile_rho << var.rho << std::endl;
    outfile_u << var.u << std::endl;
    outfile_p << var.p << std::endl;
    outfile_T << calculate_temperature(present_layer[count_of_cells - 1]) << std::endl;
    outfile_q1 << present_layer[count_of_cells - 1].rho_q1 / present_layer[count_of_cells - 1].rho << std::endl;
    outfile_q2 << present_layer[count_of_cells - 1].rho_q2 / present_layer[count_of_cells - 1].rho << std::endl;
    outfile_q3 << present_layer[count_of_cells - 1].rho_q3 / present_layer[count_of_cells - 1].rho << std::endl;

    outfile_rho.close();
    outfile_u.close();
    outfile_p.close();
    outfile_T.close();
    outfile_q1.close();
    outfile_q2.close();
    outfile_q3.close();
}

void HLLC_solver::conducting_a_reaction(double time_step) {
    // time_step - заготовка для добавления скорости реакции
    // пока просто считаю, что половина CH4 разложилась на компоненты
    for (unsigned i = 0; i < count_of_cells; ++i) {
        double T = calculate_temperature(present_layer[i]);
        if (T > reaction_temperature) {
            double velocity = present_layer[i].rho_u / present_layer[i].rho;
            double d_CH4 = reacting_constant * time_step * present_layer[i].rho_q1;
            double d_C2H2 = 0.5 * d_CH4;
            double d_H2 = 1.5 * d_CH4;

            present_layer[i].rho_q1 -= d_CH4;
            present_layer[i].rho_q2 += d_C2H2;
            present_layer[i].rho_q3 += d_H2;
            // обновляем плотность и скорость
            present_layer[i].rho = present_layer[i].rho_q1 + present_layer[i].rho_q2 + present_layer[i].rho_q3;
            present_layer[i].rho_u = present_layer[i].rho * velocity;
        }
    }
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
        // accounting for possible decay
        conducting_a_reaction(dt);
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