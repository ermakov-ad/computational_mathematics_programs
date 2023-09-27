#include <iostream>
#include <fstream>

const double L = 500.0;
const double Z = 10.0;
const double B = 50.0;
const unsigned int N = 100;
//pressure in pascals
double start_pressure = 100.0*101000.0;
double p_inj = 150.0*101000.0;
double p_prod = 50.0*101000.0;

//constants in the filtration equation
const double k = 1.0*1e-14;
const double mu = 1.0*1e-3;
const double phi = 0.2;
const double p_0 = 120.0*101000.0; //pressure in control experiments
const double rho_0 = 1000.0;    //density in control experiments
const double cf = 9.9*1e-10; //compressibility of the liquid

class layer_of_ground {
public:

    layer_of_ground(double length, double width, double height, unsigned int n_x);
    ~layer_of_ground();

    void define_boundary_pressure (double p1, double p2);
    void define_pressure_start (double p);
    void define_times (double calculation_time, double time_step);
    void start_solving();
    void save_results();

private:

    double find_rho(double p);
    void find_diagonals();
    double find_rho_plus(unsigned int index);
    double find_rho_minus(unsigned int index);

    double length;
    double width;
    double height;
    unsigned int NX;
    double tau;
    double h;
    double calculation_time;

    double *pressure;
    double *a_diag;
    double *b_diag;
    double *c_diag;
    double *d_column;
    double *alpha_koef;
    double *betta_koef;
};

layer_of_ground::layer_of_ground(double size1, double size2, double size3, unsigned int nx) {
    length = size1;
    width = size2;
    height = size3;
    NX = nx;
    h = size1 / nx;

    pressure = new double[nx];
    a_diag = new double[nx];
    b_diag = new double[nx];
    c_diag = new double[nx];
    d_column = new double[nx];
    alpha_koef = new double[nx + 1];
    betta_koef = new double[nx + 1];

    alpha_koef[0] = 0.0;
    betta_koef[0] = 0.0;
    a_diag[0] = 1.0;
    a_diag[NX-1] = 1.0;
    b_diag[0] = 0.0;
    b_diag[NX-1] = 0.0;
    c_diag[0] = 0.0;
    c_diag[NX-1] = 0.0;
    std::cout << "object created" << std::endl;
}

layer_of_ground::~layer_of_ground() {
    delete[] pressure;
    delete[] a_diag;
    delete[] b_diag;
    delete[] c_diag;
    delete[] d_column;
    delete[] alpha_koef;
    delete[] betta_koef;
    std::cout << "object delete" << std::endl;
}

void layer_of_ground::define_boundary_pressure (double p1, double p2) {
    d_column[0] = p1;
    d_column[NX-1] = p2;
    std::cout << "boundary conditions defined" << std::endl;
}

void layer_of_ground::define_pressure_start (double p) {
    for (unsigned int i = 0; i < NX; ++i)
    {
        pressure[i] = p;
    }
    std::cout << "boundary start defined" << std::endl;
}

void layer_of_ground::define_times (double calculation_time, double time_step) {
    this->tau = time_step / 10.0;
    this->calculation_time = calculation_time;
    std::cout << "calculation time found" << std::endl;
}

double layer_of_ground::find_rho (double p) {
    return rho_0 * (1.0 + cf*(p - p_0));
}

double layer_of_ground::find_rho_plus(unsigned int index) {
    if(pressure[index] >= pressure[index+1])
    {
        return find_rho(pressure[index]);
    }
    else
    {
        return find_rho(pressure[index+1]);
    }
}

double layer_of_ground::find_rho_minus(unsigned int index) {
    if(pressure[index] > pressure[index-1])
    {
        return find_rho(pressure[index]);
    }
    else
    {
        return find_rho(pressure[index-1]);
    }
}

void layer_of_ground::find_diagonals() {
    for (unsigned int i = 1; i < NX-1; ++i)
    {
        c_diag[i] = k * find_rho_minus(i) / (mu * this->h * this->h);
        b_diag[i] = k * find_rho_plus(i) / (mu * this->h * this->h);
        a_diag[i] = - c_diag[i] - b_diag[i] - phi*cf*rho_0/(this->tau);
        d_column[i] = - pressure[i]*phi*cf*rho_0/(this->tau);
    }

    alpha_koef[1] = 0.0;
    betta_koef[1] = d_column[0];
    for (unsigned int i = 1; i < NX; ++i)
    {
        alpha_koef[i+1] = -1.0 * b_diag[i] / (c_diag[i]*alpha_koef[i] + a_diag[i]);
        betta_koef[i+1] = (d_column[i] - c_diag[i]*betta_koef[i]) / (c_diag[i]*alpha_koef[i] + a_diag[i]);
    }
}

void layer_of_ground::save_results() {
    std::ofstream outfile;
    outfile.open("layer_of_ground.txt", std::ios_base::app);

    for (unsigned int i = 0; i < this->NX - 1; ++i)
    {
        outfile << pressure[i] * 1e-6 << " ";   //MPa
    }
    outfile << pressure[this->NX - 1] * 1e-6 << std::endl;

    outfile.close();
    std::cout << "results in time " << this->calculation_time << " was saved" << std::endl;
}

void layer_of_ground::start_solving() {
    std::cout << "calculation is starting" << std::endl;
    double t = 0.0;
    while (t < this->calculation_time)
    {
        t += this->tau;
        this->find_diagonals();

        pressure[this->NX - 1] = (d_column[this->NX - 1] - c_diag[this->NX - 1]*betta_koef[this->NX - 1]) / (c_diag[this->NX - 1]*alpha_koef[this->NX - 1] + a_diag[this->NX - 1]);
        for (unsigned i = this->NX - 2; i > 0; i--)
        {
            pressure[i] = pressure[i+1]*alpha_koef[i+1] + betta_koef[i+1];
        }
        pressure[0] = pressure[1]*alpha_koef[1] + betta_koef[1];
    }
    std::cout << "calculations completed" << std::endl;
}

int main() {
    layer_of_ground tube(L, Z, B, N);
    tube.define_boundary_pressure(p_inj, p_prod);
    tube.define_pressure_start(start_pressure);
    tube.define_times(2.0*3600.0, 3600.0);  //2 hours
    tube.start_solving();
    tube.save_results();
    tube.define_times(4.0*3600.0, 3600.0);  //6 hours
    tube.start_solving();
    tube.save_results();
    tube.define_times(6.0*3600.0, 3600.0);  //12 hours
    tube.start_solving();
    tube.save_results();
    tube.define_times(12.0*3600.0, 3600.0); //24 hours
    tube.start_solving();
    tube.save_results();
    tube.define_times(12.0*3600.0, 3600.0); //36 hours
    tube.start_solving();
    tube.save_results();
    tube.define_times(12.0*3600.0, 3600.0); //48 hours
    tube.start_solving();
    tube.save_results();
    tube.define_times(24.0*3600.0, 3600.0); //72 hours
    tube.start_solving();
    tube.save_results();
    tube.define_times(168.0*3600.0, 3600.0);//240 hours
    tube.start_solving();
    tube.save_results();
return 0;
}
