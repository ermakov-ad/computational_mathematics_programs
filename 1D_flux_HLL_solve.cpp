#include <iostream>
#include <fstream>
#include <cmath>

//constant physical parameters
double gamma = 1.4;
//constant environment parameters
double p_l = 1.0;
double p_r = 0.1;
double rho_l = 1.0;
double rho_r = 0.125;
double u_l = 0.75;
double u_r = 0.0;
double CFL = 0.1;
//a = sqrt(gamma*p / rho)
double a_l = sqrt(gamma * p_l / rho_l);
double a_r = sqrt(gamma * p_r / rho_r);
//p = <p>/2 + <u>/2 * <rho>/2 * <a>/2
double p_stable = 0.5*(p_l + p_r) - 0.125*(u_r - u_l)*(rho_l + rho_r)*(a_l + a_r);

//e = p/(rho*(gamma - 1))
double e_l = p_l / rho_l / (gamma - 1.0);    //specific internal energy
double e_r = p_r / rho_r / (gamma - 1.0);

double SL = std::min(u_l - a_l, u_r - a_r);
double SR = std::max(u_r + a_r, u_r + a_r);

//grid parameters
double LR = 0.7;    //length of right tube
double LL = 0.3;    //length of left tube
double h = 0.01;    //characteristic distance step
double T = 0.2;    //limit of time
double tau = CFL * h / std::max(SL, SR);    //characteristic time step <- T/random coefficient

struct U {
    double rho;
    double rho_u;
    double E;
};

//E = rho*(0.5*u^2 + e); e = p/(rho*(gamma - 1)) => p = (E/rho - u^2/2)*rho*(gamma - 1)
double pressure (U state) {
    return (state.E - 0.5*state.rho_u*state.rho_u/state.rho)*(gamma - 1.0);
}

//a = sqrt(gamma*p / rho)
double find_a (U state) {
    return sqrt(gamma * pressure(state)/ state.rho);
}

double find_SL (U left, U right) {
    return std::min(left.rho_u/left.rho - find_a(left), right.rho_u/right.rho - find_a(right));
}
double find_SR (U left, U right) {
    return std::max(left.rho_u/left.rho + find_a(left), right.rho_u/right.rho + find_a(right));
}

U flow(U start) {
    U ans;
    double p = pressure(start);
    ans.rho = start.rho_u;
    ans.rho_u = start.rho_u * start.rho_u / start.rho + p;
    ans.E = start.rho_u * (start.E + p) / start.rho;
    return ans;
}

U operator+ (U u1, U u2) {
    U ans;
    ans.rho = u1.rho + u2.rho;
    ans.rho_u = u1.rho_u + u2.rho_u;
    ans.E = u1.E + u2.E;
    return ans;
}

U operator- (U u1, U u2) {
    U ans;
    ans.rho = u1.rho - u2.rho;
    ans.rho_u = u1.rho_u - u2.rho_u;
    ans.E = u1.E - u2.E;
    return ans;
}

U operator* (double x, U var) {
    U ans;
    ans.rho = x * var.rho;
    ans.rho_u = x * var.rho_u;
    ans.E = x * var.E;
    return ans;
}

U operator/ (U var, double x) {
    U ans;
    ans.rho = var.rho / x;
    ans.rho_u = var.rho_u / x;
    ans.E = var.E / x;
    return ans;
}

U F_hll_calc(U last, U next) {
    U ans;
    double S_L = find_SL(last, next);
    double S_R = find_SR(last, next);
    ans = (S_R*flow(last) - S_L*flow(next) + S_R*S_L*(next - last)) / (S_R - S_L);
    return ans;
}

double HLL_solver(){
    //initial and boundary conditionals
    //E = rho*(0.5*u^2 + e)
    U uL, uR;
    uL.rho = rho_l;
    uL.rho_u = rho_l * u_l;
    uL.E = rho_l*e_l + 0.5*rho_l*u_l*u_l;
    uR.rho = rho_r;
    uR.rho_u = rho_r * u_r;
    uR.E = rho_r*e_r + 0.5*rho_r*u_r*u_r;

    //approximate solution of the Riemann problem
    int nt = T / tau, nl = LL / h, nr = LR / h;
    U* approximate_tube = new U[nl + nr];
    U* tube = new U[nl + nr];
    for (int i = 0; i < nl; ++i)
    {
        tube[i] = uL;
        approximate_tube[i] = uL;
    }
    for (int i = nl; i < nl + nr; ++i)
    {
        tube[i] = uR;
        approximate_tube[i] = uR;
    }

    //save state of tube for animate
    std::ofstream outfile;
    outfile.open("tube.txt", std::ios_base::app);
    for (int i = 0; i < nl + nr - 1; ++i)
    {
        if (i < nl)
        {
            outfile << p_l << ", ";
        }
        else
        {
            outfile << p_r << ", ";
        }
    }
    outfile << p_r << std::endl;

    U FL = flow(uL);
    U FR = flow(uR);
    U U_hll = (SR*uR - SL*uL + FL - FR) / (SR - SL);    //approximate state between waves

    U* flux_in_tube = new U[nl + nr];
    flux_in_tube[0] = FL;
    flux_in_tube[nl + nr - 1] = FR;
    for (int i = 1; tau*i < T; ++i)
    {
        int lim_l = tau * i * SL / h;
        int lim_r = tau * i * SR / h;
        for (int j = nl + lim_l; j <= nl + lim_r; ++j)
        {
            approximate_tube[j] = U_hll;
        }

        //recalculation state in tube
        for (int j = 0; j < nr + nl - 1; ++j)
        {
            flux_in_tube[j] = F_hll_calc(tube[j], tube[j+1]);
        }

        for (int j = 1; j < nr + nl; ++j)
        {
            tube[j] = tube[j] - tau * (flux_in_tube[j] - flux_in_tube[j-1]) / h;    //U[t+tau] = U[t] - dt/dx * (F[i+0.5] - F[i-0.5])
        }
        tube[0] = tube[1];

        for (int j = 0; j < nl + nr - 1; ++j)
        {
            outfile << pressure(tube[j]) << ", ";
        }
        outfile << pressure(tube[nl + nr - 1]) << std::endl;
    }
    outfile.close();

    delete[] tube;
    delete[] approximate_tube;
    delete[] flux_in_tube;
return U_hll.rho;
}

int main()
{
    std::cout << "SL = " << SL << "; SR = " << SR << std::endl;
    std::cout << "p in middle zone = " << p_stable << std::endl;
    double ans;
    ans = HLL_solver();
    std::cout << "rho in middle zone = " << ans;
return 0;
}
