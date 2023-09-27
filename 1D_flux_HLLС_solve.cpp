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
//a = sqrt(gamma*p / rho)
double a_l = sqrt(gamma * p_l / rho_l);
double a_r = sqrt(gamma * p_r / rho_r);
//p = <p>/2 + <u>/2 * <rho>/2 * <a>/2
double p_stable = 0.5*(p_l + p_r) - 0.125*(u_r - u_l)*(rho_l + rho_r)*(a_l + a_r);

//e = p/(rho*(gamma - 1))
double e_l = p_l / rho_l / (gamma - 1.0);    //specific internal energy
double e_r = p_r / rho_r / (gamma - 1.0);

//grid parameters
double LR = 0.7;    //length of right tube
double LL = 0.3;    //length of left tube
double h = 0.01;    //characteristic distance step
double T = 0.2;    //limit of time
double CFL = 0.1;

struct U_state {
    double rho;
    double rho_u;
    double E;
};

//E = rho*(0.5*u^2 + e); e = p/(rho*(gamma - 1)) => p = (E/rho - u^2/2)*rho*(gamma - 1)
double pressure (U_state state) {
    return (state.E - 0.5*state.rho_u*state.rho_u/state.rho)*(gamma - 1.0);
}

//a = sqrt(gamma*p / rho)
double find_a (U_state state) {
    return sqrt(gamma * pressure(state)/ state.rho);
}
double find_SL (U_state left, U_state right) {
    return std::min(left.rho_u/left.rho - find_a(left), right.rho_u/right.rho - find_a(right));
}
double find_SR (U_state left, U_state right) {
    return std::max(left.rho_u/left.rho + find_a(left), right.rho_u/right.rho + find_a(right));
}

double find_S_star (U_state left, U_state right) {
    double SL = find_SL(left, right);
    double SR = find_SR(left, right);
    double ul = left.rho_u / left.rho;
    double ur = right.rho_u / right.rho;
    return (pressure(right) - pressure(left) + left.rho_u*(SL - ul) - right.rho_u*(SR - ur)) / (left.rho*(SL - ul) - right.rho*(SR - ur));
}

U_state operator+ (U_state u1, U_state u2) {
    U_state ans;
    ans.rho = u1.rho + u2.rho;
    ans.rho_u = u1.rho_u + u2.rho_u;
    ans.E = u1.E + u2.E;
    return ans;
}

U_state operator- (U_state u1, U_state u2) {
    U_state ans;
    ans.rho = u1.rho - u2.rho;
    ans.rho_u = u1.rho_u - u2.rho_u;
    ans.E = u1.E - u2.E;
    return ans;
}

U_state operator* (double x, U_state var) {
    U_state ans;
    ans.rho = x * var.rho;
    ans.rho_u = x * var.rho_u;
    ans.E = x * var.E;
    return ans;
}

U_state operator/ (U_state var, double x) {
    U_state ans;
    ans.rho = var.rho / x;
    ans.rho_u = var.rho_u / x;
    ans.E = var.E / x;
    return ans;
}

double find_plr (U_state left, U_state right) {
    double SL = find_SL(left, right);
    double SR = find_SR(left, right);
    double s_star = find_S_star(left, right);
    double ul = left.rho_u / left.rho;
    double ur = right.rho_u / right.rho;
    return 0.5*(pressure(right) + pressure(left) + left.rho*(SL - ul)*(s_star - ul) + right.rho*(SR - ur)*(s_star - ur));
}

U_state find_D(U_state start) {
    U_state ans;
    ans.rho = 0.0;
    ans.rho_u = 1.0;
    ans.E = start.rho_u / start.rho;
    return ans;
}

U_state find_D_star(U_state left, U_state right) {
    U_state ans;
    ans.rho = 0.0;
    ans.rho_u = 1.0;
    ans.E = find_S_star(left, right);
    return ans;
}

U_state flow(U_state start) {
    U_state ans;
    ans = start.rho_u*start/start.rho + pressure(start)*find_D(start);
    return ans;
}

U_state flow_l_star(U_state left, U_state right) {
    U_state D_star = find_D_star(left, right);
    double s_star = find_S_star(left, right);
    double sl = find_SL(left, right);
    U_state FL = flow(left);
    double plr = find_plr(left, right);
    return (s_star*(sl*left - FL) + sl*plr*D_star) / (sl - s_star);
}

U_state flow_r_star(U_state left, U_state right) {
    U_state D_star = find_D_star(left, right);
    double s_star = find_S_star(left, right);
    double sr = find_SR(left, right);
    U_state FR = flow(right);
    double plr = find_plr(left, right);
    return (s_star*(sr*right - FR) + sr*plr*D_star) / (sr - s_star);
}

U_state flow_hllc(U_state left, U_state right) {
    U_state ans;
    double SL = find_SL(left, right);
    double SR = find_SR(left, right);
    double S_star = find_S_star(left, right);
    if (SL >= 0)
    {
        ans = flow(left);
    }
    else if (SL <= 0 && S_star >= 0)
    {
        ans = flow_l_star(left, right);
    }
    else if (S_star <= 0 && SR >= 0)
    {
        ans = flow_r_star(left, right);
    }
    else
    {
        ans = flow(right);
    }
    return ans;
}

int HLLC_solver(){
    //initial and boundary conditionals
    //E = rho*(0.5*u^2 + e)
    U_state UL, UR, D;
    UL.rho = rho_l;
    UL.rho_u = rho_l * u_l;
    UL.E = rho_l*e_l + 0.5*rho_l*u_l*u_l;
    UR.rho = rho_r;
    UR.rho_u = rho_r * u_r;
    UR.E = rho_r*e_r + 0.5*rho_r*u_r*u_r;

    double tau = CFL * h / std::max(find_SL(UL, UR), find_SR(UL, UR));    //characteristic time step
    int nt = T / tau, nl = LL / h, nr = LR / h;

    U_state* tube = new U_state[nl + nr];
    U_state* flux = new U_state[nl + nr];
    for (int i = 0; i < nl; ++i)
    {
        tube[i] = UL;
    }
    for (int i = nl; i < nl + nr; ++i)
    {
        tube[i] = UR;
    }

    //save state of tube for animate
    std::ofstream outfile;
    outfile.open("tube_HLLC.txt", std::ios_base::app);
    for (int i = 0; i < nl + nr - 1; ++i)
    {
        if (i < nl)
        {
            outfile << pressure(UL) << ", ";
        }
        else
        {
            outfile << pressure(UR) << ", ";
        }
    }
    outfile << pressure(UR) << std::endl;

    for (int i = 1; tau*i < T; ++i)
    {
        for (int j = 1; j < nl + nr - 1; ++j)
        {
            flux[j] = flow_hllc(tube[j], tube[j+1]);
        }
        flux[0] = flux[1];
        flux[nl + nr - 1] = flux[nl + nr - 2];

        for (int j = 1; j < nl + nr; ++j)
        {
            tube[j] = tube[j] - tau * (flux[j] - flux[j-1]) / h;    //U[t+tau] = U[t] - dt/dx * (F[i+0.5] - F[i-0.5])
        }
        tube[0] = tube[1];

        for (int j = 0; j < nr + nl - 1; ++j)
        {
            outfile << pressure(tube[j]) << ", ";
        }
        outfile << pressure(tube[nl + nr - 1]) << std::endl;
    }
    outfile.close();

    delete[] tube;
    delete[] flux;
return nt;
}

int main()
{
    int ans;
    ans = HLLC_solver();
    std::cout << ans;
return 0;
}
