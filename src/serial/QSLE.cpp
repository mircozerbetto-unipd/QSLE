/***********************************************************************************
 * QSLE v1.0 - Program to simulate the Quantum-Stochastic Liouville Equation       *
 * Copyright (C) 2024  Riccardo Cortivo, Mirco Zerbetto                            *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or any later version.                                           *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 ***********************************************************************************
 * Authors: Riccardo Cortivo, Mirco Zerbetto                                       *
 * Dipartimento di Scienze Chimiche - Universita' di Padova - Italy                *
 * E-mail: mirco.zerbetto@unipd.it                                                 *
 ***********************************************************************************/
#include "QSLE.h"

#define HARTREE 4.359744722207185e-18      // Hartree constant in J
#define EV 1.60218e-19                     // Electronvolt constant in J
#define PLANCK 6.62607015e-34              // Planck constant in J/s
#define BOLTZMANN 1.380649e-23             // Boltzmann constant in J/K
#define DALTON 1.66053906660e-27           // dalton in Kg
#define SPACE_UNIT 1.0e-10                 // The unit of measure for the space is the angstrom
#define TIME_UNIT 1.0e-15                  // The unit of measure for the time is the picosecond

const double boltzmann = BOLTZMANN /SPACE_UNIT /SPACE_UNIT/ DALTON* TIME_UNIT* TIME_UNIT;
const double energy_conversion = EV /SPACE_UNIT /SPACE_UNIT/ DALTON* TIME_UNIT* TIME_UNIT;

using namespace std;
using namespace arma;

void bad_termination(){
    cout << endl << "### Terminating QSLE because of errors ###" << endl << endl;
    exit(1);
}

void get_2C_data(const string& filename, vector<double> &vecX, vector<double> &vecY){
    ifstream fin(filename);
    if (!fin.good()){
        cout << endl << "### ERROR: File " << filename << " not found ###" <<  endl;
        bad_termination();
    }
    double x,y ;
    for (string line; getline(fin, line);){
        istringstream iss(line);
        if (iss >> x >> y){
            vecX.push_back(x);
            vecY.push_back(y);
        }
    }
    if (vecX.size()<3){
        cout << endl << "### ERROR: Bad formatting of data in file " << filename << " ###" << endl;
        bad_termination();
    }
}

void get_3C_data(const string& filename, vector<double> &vecX, vector<double> &vecY,  vector<double> &vecZ){
    ifstream fin(filename);
    if (!fin.good()){
        cout << endl << " ### ERROR: File " << filename << " not found ###" <<  endl;
        bad_termination();
    }
    double x,y,z ;
    for (string line; getline(fin, line);){
        istringstream iss(line);
        if (iss >> x >> y >> z){
            vecX.push_back(x);
            vecY.push_back(y);
            vecZ.push_back(z);
        }
    }
    if (vecX.size()<3){
        cout << endl << "### ERROR: Bad formatting of data in file " << filename << " ###" << endl;
        bad_termination();
    }
}

double norm2D_spec(const vec& X, const vec& Y, vec Z){
    long long int Ysize =  Y.n_elem;
    long long int Xsize =  X.n_elem;
    long long int Zsize =  Z.n_elem;
    vec integral_y1(Ysize);
    vec integral_x1(Xsize);
    vec integral_y2(Ysize);
    vec integral_x2(Xsize);
    for (long long int i=0; i < Xsize; i++){
        for (long long int ii=0; ii < Ysize; ii++){
            integral_y1(ii) = Z( i*Ysize + ii);
            integral_y2(ii) = Z( i*Ysize + ii + Xsize*Ysize);
        }
        mat int_trap = trapz( Y, integral_y1);
        integral_x1(i)= int_trap(0,0);
        int_trap = trapz( Y, integral_y2);
        integral_x2(i)= int_trap(0,0);
    }
    double norm= as_scalar( trapz( X, integral_x1) + trapz( X, integral_x2) ) ;
    return norm;
}


//First derivative of the Cartesian position ca respect the Z-matrix coordinate q
vec vectorArmadilloConverter(vectorOfDoubles * v){
    vec der = zeros(3); // implementation from molecule.cpp
    for (int i = 0; i < 3; i++){
        der(i) = v->at(i);
    }
    return der;
}

// Function to obtain the inertia tensor in scaled units
mat inertia(molecule *m){
    m->calculateInertiaTensor();
    auto matrix = m->getInertiaTensor();
    mat I = zeros(3,3);
    I(0,0) = (&matrix)->getRow1().at(0);
    I(1,1) = (&matrix)->getRow2().at(1);
    I(2,2) = (&matrix)->getRow3().at(2);
    // I *= DALTON * ANGSTROM * ANGSTROM;      
    return I;
}

// Function to obtain the gauge potential (adimensional)
vec gauge(molecule *m){
    mat I = inertia(m);
    int Natoms = m->getNAtoms();
    //cout<<"This is the number of active atoms "<<Natoms<<endl;
    vec gauge = zeros(3);
    for (int i = 0; i < Natoms; i++){
        // here the atom index starts from 1! See molecule.cpp
        // 4 stands for the dihedral angle
        auto der = m->getFirstDerivativeInMF(i+1, 3, 4);
        auto pos = m->getAtomPositionInMF(i+1);
        vec c_alpha_prime = vectorArmadilloConverter(&der);
        vec c_alpha = vectorArmadilloConverter(&pos);
        gauge += m->getMassOfAtom(i) * cross(c_alpha, c_alpha_prime);
        //cout << " gauge is " << gauge << endl;
    }
    gauge = inv(I) * gauge ;
    return gauge;
}

// Implementation of the metric tensor in scaled units 
double metric(molecule * m){
    int Natoms = m->getNAtoms();
    vec A = gauge(m);
    mat I = inertia(m);
    double g = 0.;
    for (int i = 0; i < Natoms; i++){
        auto der = m->getFirstDerivativeInMF( i+1, 3, 4);
        vec c_alpha_prime = vectorArmadilloConverter(&der);
        g += m->getMassOfAtom(i)  * norm(c_alpha_prime) * norm(c_alpha_prime);
    }
    g-= as_scalar(A.t() * I * A);
    return g;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}



/***************/
/* Constructor */
/***************/
QSLE::QSLE(const string& input_name)
{
    cout << setprecision(16);
    restart << setprecision(16);
    dite_input = input_name;
    base_name = "QSLE";
    restart.str("");
    PES_file = pdb_file = d01_file = reff = C = viscosity = "";
    distribution = "boltzmann" ;
    transition = "none";
    method = "FD";
    close_q = 15;
    close_p = 15;
    ts = 0.005;
    total_steps = 200000000;
    restart_step = 20000;
    normalization_step = 0;
    angle_points = 300;
    mom_points = 301;
    max_Ek = 10.;
    temperature = 300.0;
    csi = 0.;
    tau_dec = 0.;
    mom_max = 0.;
    initial_time = 0.;
    start_angle = start_mom = 0. ;
    duration_dite2 = duration_evo = duration_m = duration_FP = 0;

    ifstream fin(input_name);
    if (!fin.good()){
        cout << endl << "### ERROR: Input file " << input_name << " not found ###" <<  endl;
        bad_termination();
    }

    cout << "=======================================================================" << endl;
    cout << "Reading from input file " << input_name << endl;
    cout << "=======================================================================" <<  endl << endl;

    string x,y, y2, y3, y4;
    for (string line; getline(fin, line);){
        istringstream iss(line);
        if (iss >> x >> y ){
            transform(x.begin(), x.end(), x.begin(), [](unsigned char c) {
                return tolower(c);
            });
            if (x == "pes"){
                PES_file = y;
            }
            else if (x == "d01"){
                d01_file = y;
            }
            else if (x == "pdb"){
                pdb_file = y;
            }
            else if (x == "prefix"){
                base_name = y;
            }
            else if (x == "restart") {
                restart.str(y);
                initial_time = stod(y);
            }
            else if (x == "total_steps") {
                total_steps = stoll(y);
            }
            else if (x == "restart_step") {
                restart_step = stoll(y);
            }
            else if (x == "normalization_step") {
                normalization_step = stoll(y);
            }
            else if (x=="timestep"){
                ts = stod(y);
            }
            else if (x=="temperature"){
                temperature = stod(y);
            }
            else if (x=="ek_max"){
                max_Ek = stod(y);
            }
            else if (x=="mom_max"){
                mom_max = stod(y);
            }
            else if (x=="angle_points"){
                angle_points = stoll(y);
            }
            else if (x=="momentum_points"){
                mom_points = stoll(y);
                if (mom_points % 2 != 1){
                    cout << "### WARNING: MOMENTUM_POINTS must be an odd number >> adding 1 to the number given in input" << endl;
                    mom_points ++;
                }
            }
            else if (x=="friction"){
                csi = stod(y);
            }
            else if (x=="decoherence_time"){
                tau_dec = stod(y);
            }
            else if (x == "initial_distribution") {
                distribution = y;
                transform(y.begin(), y.end(), y.begin(), [](unsigned char c) {
                    return tolower(c);
                });
                if (y == "dirac"){
                    distribution = y;
                    if (iss >> y2 >> y3){
                        start_angle = stod(y2);
                        start_mom = stod(y3);
                    }
                    else{
                        cout << endl << "### ERROR: No starting point given for INITIAL_DISTRIBUTION DIRAC ###" << endl;
                        bad_termination();
                    }
                }
                else if (y == "boltzmann"){
                    distribution = y;
                }
            }
            else if (x == "derivation_method") {
                method = y;
                transform(method.begin(), method.end(), method.begin(), [](unsigned char c) {
                    return toupper(c);
                });
                if (method!= "LRBF" && method != "FD"){
                    cout << endl << "### ERROR: String associated to keyword DERIVATION_METHOD not recognised ###" << endl;
                    bad_termination();
                }
            }
            else if (x=="angle_domain"){
                close_q = stoi(y);
                if (close_q % 2 != 1 || close_q < 4){
                    cout << endl << "### ERROR: ANGLE_DOMAIN must be an odd number > 4" << endl;
                    bad_termination();
                }
            }
            else if (x=="momentum_domain"){
                close_p = stoi(y);
                if (close_p % 2 != 1 || close_p < 4){
                    cout << endl << "### ERROR: ANGLE_DOMAIN must be an odd number > 4" << endl;
                    bad_termination();
                }
            }
            else if (x == "force_transition") {
                transition = y;
                transform(transition.begin(), transition.end(), transition.begin(), [](unsigned char c) {
                    return tolower(c);
                });
                if (transition!= "up" && transition != "none" && transition!= "down"){
                    cout << endl << "### ERROR: String associated to keyword FORCE_TRANSITION not recognised ###" << endl;
                    bad_termination();
                }
            }   
            else if (x=="refatoms"){
                if (iss >> y2 >> y3 >> y4){
                    dihedral.push_back(stoi(y));
                    dihedral.push_back(stoi(y2));
                    dihedral.push_back(stoi(y3));
                    dihedral.push_back(stoi(y4));
                }
                else{
                    cout << endl << "### ERROR: Indexes associated to keyword REFATOMS not recognised ###!" << endl;
                    bad_termination();
                }
            }
            else if (x == "reff"){
                reff = y;
            }
            else if (x == "c"){
                C = y;
            }
            else if (x == "viscosity"){
                viscosity = y;
            }
        }
    }
    if (dihedral.size() == 0){
        cout << endl << "### ERROR: Dihedral setting not provided in the input file ###" << endl;
        bad_termination();
    }

    number_of_points = angle_points * mom_points;
    if (normalization_step == 0){
        normalization_step = restart_step;
    }


    cout << "***************************" << endl;
    cout << "* Parsing QSLE input file *" << endl;
    cout << "***************************" << endl;
    cout << "- Prefix given to output files       : " << base_name << endl;
    if (restart.str().empty()){
    }
    else{
        cout << "- Restart time                       : " << restart.str() << endl;
    }
    cout << "- Input PDB file                     : " << pdb_file << endl;
    cout << "- Input PES file                     : " << PES_file << endl;
    cout << "- Input d01 file                     : " << d01_file << endl;
    cout << "- ID of reference atoms              : " << dihedral[0] << " " << dihedral[1] << " " << dihedral[2] << " " << dihedral[3] << endl;
    cout << "- Number of points for angle axis    : " << angle_points << endl;
    cout << "- Number of points for momentum axis : " << mom_points << endl;
    if (method == "LRBF"){
        cout << "- Numerical method for derivatives   : Local Radial Basis Function" << endl;
        cout << "- Local domain for angle             : " << close_q << endl;
        cout << "- Local domain for momentum          : " << close_p << endl;
    }
    else{
        cout << "- Numerical method for derivatives   : Finite Difference" << endl;
    }
    cout << "- Timestep in femtoseconds           : " << ts << endl; 
    cout << "- Total number of steps              : " << total_steps << endl; 
    cout << "- Restart step number                : " << restart_step << endl;
    cout << "- Normalization step number          : " << normalization_step << endl;
    cout << "- Temperature                        : " << temperature << endl;
    if (mom_max == 0){
        cout << "- Maximum kinetic energy in kbT      : " << max_Ek << endl;
    }
    else{
        cout << "- Maximum momentum in scaled units   : " << mom_max << endl;
    }
    if (csi != 0.){
        cout << "- Friction in SI units               : " << csi << endl;
        csi = csi / DALTON / SPACE_UNIT / SPACE_UNIT * TIME_UNIT;
    }
    if (tau_dec != 0.){
        cout << "- Decoherence time in femtoseconds   : " << tau_dec << endl;
    }
    cout << "- Initial distribution               : " << distribution << endl;
    if (distribution == "dirac"){
        cout << "- Coordinates of central point       : (" << start_angle << "," << start_mom << ")" << endl;
    }
    cout << "- Force initial transition           : " << transition << endl;
    if (!reff.empty()){
        cout << "- Effective radius (DiTe2 setting)   : " << reff << endl;
    }
    if (!C.empty()){
        cout << "- C (DiTe2 setting)                  : " << C << endl;
    }
    if (!viscosity.empty()){
        cout << "- Viscosity (DiTe2 setting)          : " << viscosity << endl;
    }

    cout << "*********************" << endl;
    cout << "* End of QSLE input *" << endl;
    cout << "*********************" << endl << endl;

    GammaM = zeros(number_of_points * 2, number_of_points * 2);
    prob_dens = zeros(number_of_points*2);

    set_grid();

}



/*****************/
/* Deconstructor */
/*****************/
QSLE::~QSLE()
= default;



/********************************************************/
/* Setting armadillo matrices to restart the simulation */
/********************************************************/
void QSLE::recover_data()
{
    cout << "**********************************" << endl;
    cout << "* Reading data from previous run *" << endl;
    cout << "**********************************" << endl;
    GammaM.load(base_name + "_matrix");
    cout << "- " << base_name << "_matrix" << endl;
    if (initial_time != 0.) {
        prob_dens.load(base_name + "_time_" + restart.str() + ".dat", raw_ascii);
        cout << "- " << base_name << "_time_" << restart.str() << ".dat" << endl << endl;
    }
    else{
        cout << endl;
        set_potential_vectors();
        set_initial_condition();
    }

}



/********************************************************/
/* Setting armadillo matrices to restart the simulation */
/********************************************************/
bool QSLE::restart_check()
{
    bool check;
    if (restart.str().empty() ){
        check = false;
    }
    else{
        check =true;
    }
    return check;
}



/********************************************/
/* Setting the number of points of the grid */
/********************************************/
void QSLE::set_grid()
{
	if (initial_time == 0){
		set_g0();
		double g0_max = *max_element(g0.begin(), g0.end());
		if (mom_max == 0) {
			mom_max = sqrt(max_Ek * 2. * g0_max * boltzmann * temperature);
		} else {
			max_Ek = mom_max * mom_max / 2. / g0_max / boltzmann / temperature;
		}	
	}
    angle_grid = vector<double>(angle_points);
    mom_grid  = vector<double>(mom_points);
    dx = 2.0 *M_PI /angle_points;
    dy = 2.0 *mom_max/(mom_points-1);
    for (long long int i=0; i< angle_points; i++){
        angle_grid[i] = -M_PI + i* dx;
    }
    for (long long int i=0; i< mom_points; i++){
        mom_grid[i] = -mom_max + i* dy;
    }
}




/********************************/
/* Setting contravariant tensor */
/********************************/
void QSLE::set_g0()
{
    cout << "*************************************" << endl;
    cout << "* Computing generalized mass vector *" << endl;
    cout << "*************************************" << endl << endl;

    config configuration;
    std::string str;
    char* qsleHome;
    qsleHome = getenv ("QSLEINFO");
    if (qsleHome == NULL)
    {
        cout << endl << "### ERROR: the QSLEINFO envinronment variable is not set ###" << endl;
        bad_termination();
    }
    str.assign(qsleHome); str.append("/vdw.dat");
    configuration.setVdwFile(str);
    str.assign(qsleHome); str.append("/atoms_masses.dat");
    configuration.setMassFile(str);
    molecule mol(&configuration);
    mol.setIOFileName(pdb_file);
    mol.getFileName();
    mol.loadMolecule();
    mol.expressXYZinMF();
    mol.calculateDerivativesOfCM();
    mol.setMainDihedralAngle(dihedral[0],dihedral[1],dihedral[2],dihedral[3]);
    mol.buildZMatrix();

    ifstream fin("zmatrix.zmt");
    double a,b,c,d,e,f,g , L4, A4;
    for (string line; getline(fin, line);){
        istringstream iss(line);
        iss >> a >> b >> c >> d >> e >> f >> g;
        if (a==4 && b==3){
            L4= c;
            A4= e;
        }
    }
    supp_g0 = vector<double>(angle_points);
    g0  = vector<double>(angle_points);
    dg0  = vector<double>(angle_points);
    for (long long int i=0; i< angle_points; i++){
        supp_g0[i]= -M_PI + i * 2.0* M_PI/angle_points;
        mol.changeZMatrix(4,3,L4,2,A4 *M_PI /180.0, 1, -M_PI + i * 2.0* M_PI/angle_points);
        mol.expressXYZinMF();
        g0[i] = metric(&mol);
    }
    tk::spline g0_fun(supp_g0, g0);
    for (long long int i=0; i< angle_points; i++){
        dg0[i] = g0_fun.deriv(1, -M_PI + i * 2.0* M_PI/angle_points);
    }

    remove("zmatrix.zmt");
}



/************************************/
/* Setting the vectors for the PESs */
/************************************/
void QSLE::set_potential_vectors()
{
    cout << "=======================================================================" << endl;
    cout << "Reading PES from file " << PES_file << endl;
    cout << "=======================================================================" <<  endl << endl;
    get_3C_data(PES_file, supp_PES, GS, ES);
    double min_PES = *min_element(GS.begin(), GS.end());
    MS = deltaE  = vector<double>(supp_PES.size());
    //FARE CICLO PER TROVARE MS E DELTAE
    for (long long int i=0; i< supp_PES.size(); i++) {
        GS[i] = (GS[i] - min_PES ) * energy_conversion;
        ES[i] = (ES[i] - min_PES ) * energy_conversion;
        MS[i] = (GS[i] + ES[i]) / 2.0;
        deltaE[i] = (ES[i] - GS[i]);
    }
    tk::spline GS_fun(supp_PES, GS);
    tk::spline ES_fun(supp_PES, ES);
    FG = FE = deltaE_01 = vector<double>(angle_points);
    for (long long int i=0; i< angle_points; i++){
        deltaE_01[i] = GS_fun(angle_grid[i]) - ES_fun(angle_grid[i]);
        FG[i]= - GS_fun.deriv(1, angle_grid[i]);
        FE[i]= - ES_fun.deriv(1, angle_grid[i]);
    }
}



/********************************/
/* Setting friction coefficient */
/********************************/
void QSLE::set_friction()
{
    if (csi==0){
        auto start_for = chrono::high_resolution_clock::now();

        cout << "******************************************************************" << endl;
        cout << "* Using DiTe2 to compute friction coefficient for dihedral angle *" << endl;
        cout << "******************************************************************" << endl;
          
        cout << "   ___ _            _     ___  _ _____    ___ " << endl;
        cout << "  / __| |_ __ _ _ _| |_  |   \\(_)_   _|__|_  )" << endl;
        cout << "  \\__ \\  _/ _` | '_|  _| | |) | | | |/ -_)/ / " << endl;
        cout << "  |___/\\__\\__,_|_|  \\__| |___/|_| |_|\\___/___|" << endl << endl;
        
        // Parse the input file
        inputParser ip;
        ip.parseInputFile(dite_input);
        DiTe2 diff(&ip);
        diff.calculateDiffusionTensor(base_name +"_DiTe2_friction.dat", base_name +"_DiTe2_diffusion.dat");
        csi= 1.e-20 * diff.FrictionDihedral;

        cout << "   ___         _    ___  _ _____    ___ " << endl;
        cout << "  | __|_ _  __| |  |   \\(_)_   _|__|_  )" << endl;
        cout << "  | _|| ' \\/ _` |  | |) | | | |/ -_)/ / " << endl;
        cout << "  |___|_||_\\__,_|  |___/|_| |_|\\___/___|" << endl << endl;

        remove("zmatrix.zmt");

        cout << "=======================================================================" << endl;
        cout << "Computed friction coefficient in SI units: " << csi << endl;
        cout << "DiTe2 results are stored in:" << endl;
	    cout << "- " << base_name << "_DiTe2_friction.dat" <<  endl;
	    cout << "- " << base_name << "_DiTe2_diffusion.dat" << endl ;
        cout << "=======================================================================" <<  endl << endl;

        cout << "=======================================================================" << endl;
	    cout << "For more information about DiTe2 program:" << endl;
	    cout << "- https://doi.org/10.1002/jcc.25742"  << endl;
	    cout << "- https://wwwdisc.chimica.unipd.it/mirco.zerbetto/?p=172"  << endl;
	    cout << "=======================================================================" << endl << endl;

        auto end_for = chrono::high_resolution_clock::now();
        auto duration_for = chrono::duration_cast<chrono::milliseconds>(end_for - start_for);
        duration_dite2 = duration_for.count();

        csi = csi / DALTON / SPACE_UNIT / SPACE_UNIT * TIME_UNIT;

    }


}



/***********************************/
/* Setting the vectors for the d01 */
/***********************************/
void QSLE::set_nonadiabatic_coeff()
{
    cout << "=======================================================================" << endl;
    cout << "Reading d01 from file " << d01_file << endl;
    cout << "=======================================================================" <<  endl << endl;
    get_2C_data(d01_file, supp_d01, d01);
}



/*************************************************/
/* Computing the x^2 coefficient of the parabola */
/*************************************************/
/*
double QSLE::parabolic_fitting()
{
    // choice of the angle range for the parabolic fitting, based on the place where d01 has its highest absolute value
    double max_angle, diff_angle_tmp;
    int max_d01_abs_element;
    int CI_grid_element;
    double max_d01 = *max_element(d01.begin(), d01.end());
    double min_d01 = *min_element(d01.begin(), d01.end());
    if ( abs(max_d01) > abs(min_d01)){
        max_d01_abs_element = max_element(d01.begin(), d01.end()) - d01.begin();
    }
    else{
        max_d01_abs_element = min_element(d01.begin(), d01.end()) - d01.begin();
    }
    max_angle = supp_d01[max_d01_abs_element];

    diff_angle_tmp= 2* M_PI;
    CI_grid_element=0;
    for (int i=0; i< supp_PES.size(); i++){
        if (abs(supp_PES[i]- max_angle) < diff_angle_tmp){
            diff_angle_tmp = abs(supp_PES[i]- max_angle);
            CI_grid_element= i;
        }
    }
    
    tk::spline MS_fun(supp_PES, MS);

    return MS_fun.deriv(2, supp_PES[CI_grid_element]);  // the resulting value is the coefficient of x^2 (l'apertura della parabola) that is useful to compute the decoherence time
}
*/


/*********************************************************************/
/* Computing the decoherence time t_dec based on the parabola method */
/*********************************************************************/
/*
void QSLE::compute_parabolic_decoherence_time()
{
    int i;
    double k, g0_CI, l1, l2, a1, a2;
    // Find the x^2 coefficient of the parabola
    k = parabolic_fitting();
    // Find g0 near the conical intersection
    double angle_CI, d01_tmp, diff_angle_tmp, sqrt_arg;
    for (i= 0; i< d01.size(); i++)
    {
        if (d01_tmp < abs(d01[i])){
            d01_tmp = abs(d01[i]);
            angle_CI = supp_d01[i];
        }
    }
    diff_angle_tmp= 2* M_PI;
    for (i=0; i< g0.size(); i++)
    {
        if( abs(supp_g0[i]-angle_CI) < diff_angle_tmp){
            diff_angle_tmp = abs(supp_g0[i]-angle_CI);
            g0_CI = g0[i];
        }
    }

    // Computing the decoherence time for different values of k
    if (k <= 0.0){
        tau_dec= g0_CI/csi;
    }
    else{
        sqrt_arg = csi*csi/g0_CI/g0_CI - 4.0*k/g0_CI;
        if (sqrt_arg <= 0.0){
            tau_dec = g0_CI/csi/2.0;
        }
        else{
            l1= 0.5*csi/g0_CI + 0.5* sqrt(sqrt_arg);
            l2= 0.5*csi/g0_CI - 0.5* sqrt(sqrt_arg);
            a1= l1*l1/(-k/g0_CI+l1*l1);
            a2= -k/g0_CI/(-k/g0_CI+l1*l1);

            // Bisection method to find at which time the correlation function is equal to 1/e
            double a, b, c, corr_a, corr_c, tolerance;
            a= 0.0;
            b= 3.0/l1;
            tolerance = 1.0*pow(10.0, -6.0);
            corr_a = 1.0;
            while ((b - a) >= tolerance) {
                // Find the midpoint
                c = (a + b) / 2;
                corr_c = a1 * exp(- c *l1) + a2 * exp(-c * l2) - exp(-1.0);
                // Check if the root is found
                if ( corr_c == 0.0) {
                    break;
                }
                    // Update the interval
                else if (corr_c * corr_a < 0) {
                    b = c;
                } else {
                    a = c;
                    corr_a = corr_c;
                }
            }

            tau_dec= c;

        }
    }
}
*/


/****************************************/
/* Computing the decoherence time t_dec */
/****************************************/
void QSLE::compute_decoherence_time()
{
    long long int i;
    double g0_CI;
    // Find g0 near the conical intersection
    double angle_CI, d01_tmp, diff_angle_tmp, sqrt_arg;
    for (i= 0; i< d01.size(); i++)
    {
        if (d01_tmp < abs(d01[i])){
            d01_tmp = abs(d01[i]);
            angle_CI = supp_d01[i];
        }
    }
    diff_angle_tmp= 2* M_PI;
    for (i=0; i< g0.size(); i++)
    {
        if( abs(supp_g0[i]-angle_CI) < diff_angle_tmp){
            diff_angle_tmp = abs(supp_g0[i]-angle_CI);
            g0_CI = g0[i];
        }
    }

    tau_dec= g0_CI/csi;

    cout << "=======================================================================" << endl;
    cout << "Computed decoherence time in femtoseconds: " << tau_dec << endl;
    cout << "=======================================================================" <<  endl << endl;

}



/*****************************************/
/* Computing the transition coefficients */
/*****************************************/
void QSLE::compute_transition_coefficients()
{
    auto start_for = chrono::high_resolution_clock::now();
    if (tau_dec == 0.){
        //compute_parabolic_decoherence_time();
        compute_decoherence_time();
    }

    cout << "*****************************************" << endl;
    cout << "* Computing transition rates (m01, m10) *" << endl;
    cout << "*****************************************" << endl;

    m01= m10 = vector<double>(number_of_points);
    tk::spline MS_fun(supp_PES, MS);
    tk::spline deltaE_fun(supp_PES, deltaE);
    tk::spline d01_fun(supp_d01, d01);
    tk::spline g0_fun(supp_g0, g0);

    double timestep = tau_dec/10000.0;
    double htag = PLANCK / 2.0 / M_PI / SPACE_UNIT /SPACE_UNIT /DALTON * TIME_UNIT;

    long long int count = 0;

    cout << "+ Building jump rates = 0%" << endl;

    #pragma omp parallel for
    for (long long int i=0; i < number_of_points; i++){

        #pragma omp atomic
        count++;

        if (count%(number_of_points/10)==0){
            cout << "+ Building jump rates = " << 100*count/number_of_points << "%" << endl;
        }

        double angle_init, mom_init;
        angle_init= angle_grid[i/ mom_points];
        mom_init = mom_grid[i % mom_points];

        if (angle_init < supp_d01[0] || angle_init> supp_d01.back() || mom_init ==0. || d01_fun(angle_init) == 0. ){
            m01[i] = m10[i] = 0.;
        }
        else {
            double angle_old, angle_tmp, mom_tmp, v_tmp, W_tmp, omega0_old, mom_up2, mom_down2, time, D0, D1_tmp, M_tmp, M_old, m01_tmp, m10_tmp;
            mom_up2 = mom_init*mom_init - 2.0* deltaE_fun(angle_init) *g0_fun(angle_init);
            if (mom_up2 >= 0.){
                mom_tmp = sgn(mom_init)*sqrt(mom_init*mom_init - deltaE_fun(angle_init) *g0_fun(angle_init));
                time=0.;
                m01_tmp=0.;
                angle_old =angle_init;
                D0= d01_fun(angle_init)* mom_init /g0_fun(angle_init);
                D1_tmp= d01_fun(angle_init)* mom_tmp /g0_fun(angle_init);
                M_old= 2.0* D0* D1_tmp;
                W_tmp=0.;
                omega0_old = deltaE_fun(angle_init);
                if (abs(supp_d01[0]) == abs(supp_d01.back()) && abs(supp_d01[0]) >= 3.14){
					while (time < 7.0* tau_dec){
						time += timestep;
						angle_tmp = angle_old - timestep * mom_tmp/g0_fun(angle_old) + 0.5 * timestep*timestep* MS_fun.deriv(1, angle_old) /g0_fun(angle_old);
						if (angle_tmp < -M_PI){
							angle_tmp +=  2.0*M_PI;
						}
						if (angle_tmp > M_PI){
							angle_tmp -= 2.0*M_PI;
						}
						if (timestep*omega0_old > M_PI_4){
							break;
						}
						v_tmp = mom_tmp/g0_fun(angle_old) + 0.5 * timestep * (MS_fun.deriv(1, angle_old) /g0_fun(angle_old) + MS_fun.deriv(1, angle_tmp) /g0_fun(angle_tmp));
						mom_tmp = v_tmp *g0_fun(angle_tmp);
						D1_tmp= d01_fun(angle_tmp)* mom_tmp /g0_fun(angle_tmp);
						W_tmp= W_tmp + (omega0_old +deltaE_fun(angle_tmp)/htag) * timestep/2.0;
						M_tmp = 2.0 * D0 * D1_tmp * exp(- time/ tau_dec) * cos(W_tmp);
						m01_tmp = m01_tmp + (M_old + M_tmp)*timestep/2.0;
						angle_old = angle_tmp;
						omega0_old = deltaE_fun(angle_old)/htag;
						M_old= M_tmp;
					}
				}
				else{
					while (time < 7.0* tau_dec){
						time += timestep;
						angle_tmp = angle_old - timestep * mom_tmp/g0_fun(angle_old) + 0.5 * timestep*timestep* MS_fun.deriv(1, angle_old) /g0_fun(angle_old);
						if (angle_tmp < -M_PI){
							angle_tmp +=  2.0*M_PI;
						}
						if (angle_tmp > M_PI){
							angle_tmp -= 2.0*M_PI;
						}
						if (angle_tmp < supp_d01[0] || angle_tmp > supp_d01.back()){
							break;
						}
						if (timestep*omega0_old > M_PI_4){
							break;
						}
						v_tmp = mom_tmp/g0_fun(angle_old) + 0.5 * timestep * (MS_fun.deriv(1, angle_old) /g0_fun(angle_old) + MS_fun.deriv(1, angle_tmp) /g0_fun(angle_tmp));
						mom_tmp = v_tmp *g0_fun(angle_tmp);
						D1_tmp= d01_fun(angle_tmp)* mom_tmp /g0_fun(angle_tmp);
						W_tmp= W_tmp + (omega0_old +deltaE_fun(angle_tmp)/htag) * timestep/2.0;
						M_tmp = 2.0 * D0 * D1_tmp * exp(- time/ tau_dec) * cos(W_tmp);
						m01_tmp = m01_tmp + (M_old + M_tmp)*timestep/2.0;
						angle_old = angle_tmp;
						omega0_old = deltaE_fun(angle_old)/htag;
						M_old= M_tmp;
					}
				}
                if (m01_tmp < 10e-5){  //10e-5 is based on the fact that it takes almost a ps to relax a 1% of population ( p=p0*e^(-m*t) considering in the differential equation only the transition coefficient )
                    m01_tmp =0;
                }
                m01[i] = m01_tmp;
            }
            else{
                m01[i] = 0.;
            }

            mom_down2 = mom_init*mom_init + 2.0* deltaE_fun(angle_init) *g0_fun(angle_init);
            mom_tmp = sgn(mom_init)*sqrt(mom_init*mom_init +  deltaE_fun(angle_init) *g0_fun(angle_init));
            time=0;
            m10_tmp=0;
            angle_old =angle_init;
            D0= d01_fun(angle_init)* mom_init /g0_fun(angle_init);
            D1_tmp= d01_fun(angle_init)* mom_tmp /g0_fun(angle_init);
            M_old= 2.0* D0* D1_tmp;
            W_tmp=0.;
            omega0_old = deltaE_fun(angle_init);
            /*
            if (i== number_of_points/2 + mom_points/2 ){
                tempo.push_back(angle_init);
                testM.push_back(mom_init);
                tempo.push_back(0);
                testM.push_back(M_old);
            }
            */
            while (time < 7.0* tau_dec){
                time += timestep;
                angle_tmp = angle_old - timestep * mom_tmp/g0_fun(angle_old) + 0.5 * timestep*timestep* MS_fun.deriv(1, angle_old) /g0_fun(angle_old);
                if (angle_tmp < -M_PI){
                    angle_tmp +=  2.0*M_PI;
                }
                if (angle_tmp > M_PI){
                    angle_tmp -= 2.0*M_PI;
                }
                if (angle_tmp < supp_d01[0] || angle_tmp > supp_d01.back()){
                    break;
                }
                v_tmp = mom_tmp/g0_fun(angle_old) + 0.5 * timestep * (MS_fun.deriv(1, angle_old) /g0_fun(angle_old) + MS_fun.deriv(1, angle_tmp) /g0_fun(angle_tmp));
                mom_tmp = v_tmp *g0_fun(angle_tmp);
                D1_tmp= d01_fun(angle_tmp)* mom_tmp /g0_fun(angle_tmp);
                W_tmp= W_tmp + (omega0_old +deltaE_fun(angle_tmp)/htag) * timestep/2.0;
                M_tmp = 2.0 * D0 * D1_tmp * exp(- time/ tau_dec) * cos(W_tmp);
                m10_tmp = m10_tmp + (M_old + M_tmp)*timestep/2.0;
                angle_old = angle_tmp;
                omega0_old = deltaE_fun(angle_old)/htag;
                M_old= M_tmp;
                /*
                if (i== number_of_points/2 + mom_points/2 ){
                    tempo.push_back(time);
                    testM.push_back(M_old);
                }
                */
            }
            if (m10_tmp < 10e-5){
                m10_tmp =0;
            }
            m10[i] = m10_tmp;

        }
    }

    auto end_for = chrono::high_resolution_clock::now();
    auto duration_for = chrono::duration_cast<chrono::milliseconds>(end_for - start_for);
    duration_m = duration_for.count();
}



/********************************/
/* Computing the c coefficients */
/********************************/
void QSLE::find_c()
{
    cout << "**********************************************" << endl;
    cout << "* Computing shape parameters for LRBF scheme *" << endl;
    cout << "**********************************************" << endl;
                
    Phi_q = zeros(close_q, close_q);
    Phi_p = zeros(close_p, close_p);

    vec eigvalq, eigvalp ;
    mat eigvecq, eigvecp;
    double c = 0.;                         
    double rate = 0.;
    while (rate< pow(10,13) || rate> pow(10,15)){
        c += 0.00001;
        for (int ii = 0; ii < close_q ; ii++) {
            for (int jj= ii; jj< close_q; jj++) {
                Phi_q(jj, ii)= exp(- c* (ii-jj) * (ii-jj) ) ;
                Phi_q(ii, jj)= Phi_q(jj, ii);
            }
        }
        eig_sym(eigvalq, eigvecq, Phi_q, "std");
        rate= eigvalq.max()/eigvalq.min();
    }
    cq = c /dx/dx;
    cout << "- cq = " << cq << endl;
           
    if (close_q != close_p){
        c = 0.;                         
        rate = 0.;
        while (rate< pow(10,13) || rate> pow(10,15)){
            c += 0.00001;
            for (int ii = 0; ii < close_p ; ii++) {
                for (int jj= ii; jj< close_p; jj++) {
                    Phi_p(jj, ii)= exp(- c* (ii-jj) * (ii-jj) ) ;
                    Phi_p(ii, jj)= Phi_p(jj, ii);
                }
            }
            eig_sym(eigvalp, eigvecp, Phi_p, "std");
            rate= eigvalp.max()/eigvalp.min();
        }
    }
    else{
        for (int ii = 0; ii < close_p ; ii++) {
            for (int jj= ii; jj< close_p; jj++) {
                Phi_p(jj, ii)= exp(- c* (ii-jj) * (ii-jj) ) ;
                Phi_p(ii, jj)= Phi_p(jj, ii);
            }
        }
    }     
    cp = c /dy/dy;
    cout << "- cp = " << cp << endl << endl;  

}



/*****************************/
/* Computing the QSLE matrix */
/*****************************/
void QSLE::compute_QSLE_matrix()
{
    auto start_for = chrono::high_resolution_clock::now();

    if (method == "LRBF"){
        find_c();
    }

    cout << "**********************************************" << endl;
    cout << "* Computing generalized Fokker-Planck matrix *" << endl;
    cout << "**********************************************" << endl;

    if (method == "LRBF"){
        compute_LRBF_matrix();
    }
    else{
        compute_FD_matrix();
    }

    GammaM.save(base_name + "_matrix");
    cout << endl << "=======================================================================" << endl;
    cout << "Fokker-Planck matrix saved in " << base_name << "_matrix"  << endl;
    cout << "=======================================================================" << endl << endl;

    auto end_for = chrono::high_resolution_clock::now();
    auto duration_for = chrono::duration_cast<chrono::milliseconds>(end_for - start_for);
    duration_FP = duration_for.count();
}




/****************************************/
/* Computing the QSLE matrix using LRBF */
/****************************************/
void QSLE::compute_LRBF_matrix()
{
    cout << "+ Building Fokker-Planck matrix = 0%" << endl;
    long long int count =0;

    #pragma omp parallel for
    for (long long int i=0; i < number_of_points; i++){

        #pragma omp atomic
        count++;

        if (count%(number_of_points/10)==0){
            cout << "+ Building Fokker-Planck matrix = " << 100*count/number_of_points << "%" << endl;
        }

        long long int index_angle, index_mom;
        index_angle = i / mom_points;
        index_mom = i % mom_points;
        double q_i = angle_grid[index_angle];
        double p_i = mom_grid[index_mom];

        long long int j;
        double distance;                                                 // distance between two points of the grid
        long long int dist_grid;
        double esponential;

        vec wq;
        vec functq(close_q);
        vec LRBF_qdomain(close_q);
        for (int ii = 0; ii < close_q ; ii++) {
            dist_grid = ii -close_q/2;
            if (index_angle + dist_grid < 0 ){
                LRBF_qdomain(ii) = index_mom + mom_points* (angle_points + index_angle + dist_grid);
            }
            else if (index_angle + dist_grid > angle_points -1){
                LRBF_qdomain(ii) = index_mom + mom_points* (index_angle + dist_grid - angle_points);
            }
            else{
                LRBF_qdomain(ii) = i + mom_points* dist_grid;
            }
            distance =  dist_grid * dx;
            esponential = exp(- cq* distance * distance);
            functq(ii)=  esponential * p_i / g0[index_angle] * 2.0 *cq * distance ;
        }
        wq= solve(Phi_q, functq);
        for (int ii = 0; ii < close_q ; ii++) {
            j= LRBF_qdomain(ii);
            GammaM(i,j) -= wq(ii);
            GammaM(i+ number_of_points,j+ number_of_points) -= wq(ii);
        }

        if (i % mom_points != 0 && i % mom_points != mom_points-1) {
            vec wpG;
            vec wpE;
            vec functpG(close_p);
            vec functpE(close_p);
            vec LRBF_pdomain(close_p);
            long long int start_index;
            if (index_mom -close_p/2 < 0 ){
                start_index = 0;
            }
            else if (index_mom +close_p/2 > mom_points -1){
                start_index = mom_points - close_p;
            }
            else{
                start_index = index_mom -close_p/2;
            }
            for (int ii = 0; ii < close_p ; ii++) {
                dist_grid = start_index + ii - index_mom;
                distance =  dist_grid * dy;
                LRBF_pdomain(ii) = start_index + ii + index_angle * mom_points;
                esponential = exp(- cp* distance * distance );
                functpG(ii)=  esponential * ( 2.0 * cp * distance * ( FG[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi * ( p_i/g0[index_angle] + boltzmann* temperature * 2.*cp*distance ) ) + 2.*cp*csi*boltzmann*temperature ) ;
                functpE(ii)=  esponential * ( 2.0 * cp * distance * ( FE[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi * ( p_i/g0[index_angle] + boltzmann* temperature * 2.*cp*distance ) ) + 2.*cp*csi*boltzmann*temperature ) ;
            }
            wpG= solve(Phi_p, functpG);
            wpE= solve(Phi_p, functpE);
            for (int ii = 0; ii < close_p ; ii++) {
                j= LRBF_pdomain(ii);
                GammaM(i,j) -= wpG(ii);
                GammaM(i+number_of_points,j+number_of_points) -= wpE(ii);
            }
            GammaM(i,i) += csi/g0[index_angle];
            GammaM(i+number_of_points, i+number_of_points) += csi/g0[index_angle];

        } else if (i % mom_points == 0) {
            GammaM(i, i + 1) = -(FG[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i + dy)) / 2.0 / dy -
                                 boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + 1 + number_of_points) = -
                    (FE[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i + dy)) / 2.0 / dy -
                    boltzmann * temperature * csi / dy / dy;
            GammaM(i, i) = -boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + number_of_points) = -boltzmann * temperature * csi / dy / dy;
        } else if (i % mom_points == mom_points-1) {
            GammaM(i, i - 1) = (FG[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i - dy)) / 2.0 / dy -
                                boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i - 1 + number_of_points) = 
                    (FE[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i - dy)) / 2.0 / dy -
                    boltzmann * temperature * csi / dy / dy;
            GammaM(i, i) = -boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + number_of_points) = -boltzmann * temperature * csi / dy / dy;
        }

        double radice;                         // useful to calculate the point after the moment jump
        double Pbeta;                          // useful to calculate the point after the moment jump
        long long int supp;
        long long int iteration;
        if (m01[i] != 0) {
            radice = p_i * p_i + 2 * deltaE_01[index_angle] * g0[index_angle];
            iteration = close_p/2;
            Pbeta = p_i / abs(p_i) * sqrt(radice);
            if (-mom_max + iteration*dy > Pbeta) {
            }else if (mom_max - iteration * dy < Pbeta) {
                iteration = mom_points - 1 - close_p/2;
            } else {
                for (long long int ii = close_p/2; ii < mom_points - close_p/2; ii++) {
                    if (-mom_max + ii * dy - dy / 2.0 <= Pbeta && Pbeta <= -mom_max + ii * dy + dy / 2.0) {
                        break;
                    }
                    iteration = iteration + 1;
                }
            }
            supp = i - i % mom_points + iteration;
            GammaM(i, supp + number_of_points) += m01[i];
            GammaM(i, i) -= m01[i];
        }
        if (m10[i] != 0) {
            radice = p_i * p_i - 2 * deltaE_01[index_angle] * g0[index_angle];
            iteration = close_p/2;
            Pbeta = p_i / abs(p_i) * sqrt(radice);
            if (-mom_max + iteration*dy > Pbeta) {
            } else if (mom_max - iteration * dy < Pbeta) {
                iteration = mom_points - 1 - close_p/2;
            } else {
                for (long long int ii = close_p/2; ii < mom_points - close_p/2; ii++) {
                    if (-mom_max + ii * dy - dy / 2.0 <= Pbeta && Pbeta <= -mom_max + ii * dy + dy / 2.0) {
                        break;
                    }
                    iteration = iteration + 1;
                }
            }
            supp = i - i % mom_points + iteration;
            GammaM(i + number_of_points, supp) += m10[i];
            GammaM(i + number_of_points, i + number_of_points) -= m10[i];
        }
    }
}



/************************************************************/
/* Computing the QSLE matrix using Finite Difference method */
/************************************************************/
void QSLE::compute_FD_matrix()
{
    cout << "+ Building Fokker-Planck matrix = 0%" << endl;
    long long int count =0;

    #pragma omp parallel for
    for (long long int i=0; i < number_of_points; i++){

        #pragma omp atomic
        count++;

        if (count%(number_of_points/10)==0){
            cout << "+ Building Fokker-Planck matrix = " << 100*count/number_of_points << "%" << endl;
        }
        double q_i = angle_grid[i/ mom_points];
        double p_i = mom_grid[i % mom_points];
        long long int index_angle, index_dx, index_sx;
        index_angle = i / mom_points;
        index_sx = i - mom_points;
        index_dx = i + mom_points;
        if (index_sx < 0) {
            index_sx = (angle_points - 1) * mom_points + i;
        } else if (index_dx >= number_of_points) {
            index_dx = i % mom_points;
        }
        GammaM(i, index_dx) = GammaM(i + number_of_points, index_dx + number_of_points) = -p_i / g0[index_angle] / 2.0 / dx;
        GammaM(i, index_sx) = GammaM(i + number_of_points, index_sx + number_of_points) = p_i / g0[index_angle] / 2.0 / dx;
        if (i % mom_points != 0 && i % mom_points != mom_points-1) {
            GammaM(i, i + 1) = -(FG[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i + dy)) / 2.0 / dy +
                               boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + 1 + number_of_points) =
                    -(FE[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i + dy)) / 2.0 / dy +
                    boltzmann * temperature * csi / dy / dy;
            GammaM(i, i - 1) = (FG[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i - dy)) / 2.0 / dy +
                               boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i - 1 + number_of_points) =
                    (FE[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i - dy)) / 2.0 / dy +
                    boltzmann * temperature * csi / dy / dy;
            GammaM(i, i) = -2.0 * boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + number_of_points) = -2.0 * boltzmann * temperature * csi / dy / dy;
        } else if (i % mom_points == 0) {
            GammaM(i, i + 1) = -(FG[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i + dy)) / 2.0 / dy -
                                 boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + 1 + number_of_points) = -
                    (FE[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i + dy)) / 2.0 / dy -
                    boltzmann * temperature * csi / dy / dy;
            GammaM(i, i) = -boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + number_of_points) = -boltzmann * temperature * csi / dy / dy;
        } else if (i % mom_points == mom_points-1) {
            GammaM(i, i - 1) = (FG[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i - dy)) / 2.0 / dy -
                                boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i - 1 + number_of_points) = 
                    (FE[index_angle] + dg0[index_angle]*p_i*p_i/2.0/g0[index_angle]/g0[index_angle] - csi / g0[index_angle] * (p_i - dy)) / 2.0 / dy -
                    boltzmann * temperature * csi / dy / dy;
            GammaM(i, i) = -boltzmann * temperature * csi / dy / dy;
            GammaM(i + number_of_points, i + number_of_points) = -boltzmann * temperature * csi / dy / dy;
        }

        double radice;                         // useful to calculate the point after the moment jump
        double Pbeta;                          // useful to calculate the point after the moment jump
        long long int iteration;
        long long int supp;
        if (m01[i] != 0) {
            radice = p_i * p_i + 2 * deltaE_01[index_angle] * g0[index_angle];
            iteration = 1;
            Pbeta = p_i / abs(p_i) * sqrt(radice);
            if (-mom_max + dy > Pbeta) {
            } else if (mom_max - dy < Pbeta) {
                iteration = mom_points - 2;
            } else {
                for (long long int ii = 1; ii < mom_points-1; ii++) {
                    if (-mom_max + ii * dy - dy / 2.0 <= Pbeta && Pbeta <= -mom_max + ii * dy + dy / 2.0) {
                        break;
                    }
                    iteration = iteration + 1;
                }
            }
            supp = i - i % mom_points + iteration;
            GammaM(i, supp + number_of_points) += m01[i];
            GammaM(i, i) -= m01[i];
        }
        if (m10[i] != 0) {
            radice = p_i * p_i - 2 * deltaE_01[index_angle] * g0[index_angle];
            iteration = 1;
            Pbeta = p_i / abs(p_i) * sqrt(radice);
            if (-mom_max +dy > Pbeta) {
            } else if (mom_max -dy < Pbeta) {
                iteration = mom_points - 2;
            } else {
                for (long long int ii = 1; ii < mom_points-1; ii++) {
                    if (-mom_max + ii * dy - dy / 2.0 <= Pbeta && Pbeta <= -mom_max + ii * dy + dy / 2.0) {
                        break;
                    }
                    iteration = iteration + 1;
                }
            }
            supp = i - i % mom_points + iteration;
            GammaM(i + number_of_points, supp) += m10[i];
            GammaM(i + number_of_points, i + number_of_points) -= m10[i];
        }
    }
}


/*********************************/
/* Setting the initial condition */
/*********************************/
void QSLE::set_initial_condition()
{
    cout << "*********************************" << endl;
    cout << "* Setting the initial condition *" << endl;
    cout << "*********************************" << endl;

    if (distribution == "boltzmann") {
        tk::spline GS_fun(supp_PES, GS);
        tk::spline ES_fun(supp_PES, ES);
        for (long long int i = 0; i < number_of_points; i++) {
            prob_dens(i) = exp(-(mom_grid[i % mom_points] * mom_grid[i % mom_points] / 2.0 / g0[i / mom_points] +
                             GS_fun(angle_grid[i/mom_points])) / boltzmann / temperature);
            prob_dens(i+ number_of_points) = exp(-(mom_grid[i % mom_points] * mom_grid[i % mom_points] / 2.0 / g0[i / mom_points] +
                             GS_fun(angle_grid[i/mom_points])) / boltzmann / temperature);
        }
        cout << "- Boltzmann configuration applied" << endl;
    }
    else if (distribution == "dirac"){
        for (long long int i = 0; i < number_of_points; i++) {
            double distance = abs(angle_grid[i / mom_points] - start_angle);
            if (distance > M_PI){
                distance = 2.*M_PI - distance;
            }
            prob_dens(i) = exp(- (mom_grid[i % mom_points] - start_mom)*(mom_grid[i % mom_points] - start_mom)/dy/dy - distance*distance/dx/dx);
        }
        cout << "- Dirac configuration (" << start_angle << "," << start_mom <<") applied " << endl;
    }
    else{
        cout << "- Reading from " << distribution << endl;
        ifstream fin(distribution);
        if (!fin.good()) {
            cout << endl << "### ERROR: Configuration associated to keyword INITIAL_DISTRIBUTION has not been recognised (wrong filename or keyword) ###" << endl;
            bad_termination();
        }
        double value;
        long long int iteration = 0;
        for (string line; getline(fin, line);){
            istringstream iss(line);
            if (iss >> value ){
                if (iteration == number_of_points*2){
                    cout << endl << "### ERROR: Mismatch between the number of points given in input and the number of points in the input file " << distribution <<  " ###" << endl;
                    bad_termination();
                }
                prob_dens(iteration) = value;
                iteration ++;
            }
        }
        if (iteration < number_of_points*2){
            cout << "### ERROR: Mismatch between the number of points given in input and the number of points in the input file " << distribution <<  " ###" << endl;
            bad_termination();
        }
    }
    if (transition == "up"){
        for (long long int i=0; i < number_of_points; i++){
            prob_dens(i+ number_of_points) = prob_dens(i+ number_of_points) + prob_dens(i);
            prob_dens(i) =0.;
        }
        cout << "- 0 --> 1 vertical transition applied " << endl;
    }
    else if (transition == "down"){
        for (long long int i=0; i < number_of_points; i++){
            prob_dens(i) = prob_dens(i+ number_of_points) + prob_dens(i);
            prob_dens(i+ number_of_points) =0.;
        }
        cout << "- 1 --> 0 vertical transition applied " << endl;
    }
    else{
        cout << "- No vertical transition applied " << endl;
    }

    cout << endl;

    normalization = norm2D_spec(angle_grid, mom_grid, prob_dens);
    prob_dens = prob_dens/ normalization;

    if (initial_time==0){
        restart.str("");
        restart << initial_time;
        print_restart();
    }
}


/******************/
/* Time evolution */
/******************/
void QSLE::time_evolution()
{
    cout << "***************************" << endl;
    cout << "* Starting time evolution *" << endl;
    cout << "***************************" << endl;
    cout << "+ Time: " << initial_time << " fs" << endl;
    

    double time;
    vec k1( number_of_points*2 );
    vec k2( number_of_points*2);
    vec k3( number_of_points*2);
    vec k4( number_of_points*2);

    integral_Y1 =zeros(mom_points);
    integral_X1=zeros(angle_points);
    integral_Y2=zeros(mom_points);
    integral_X2=zeros(angle_points);
    X = angle_grid;
    Y = mom_grid;
    

    auto start_for = chrono::high_resolution_clock::now();
    for (long long int t=1; t< total_steps+1; t++){
        //Propagation using Runge Kutta
        k1= GammaM * prob_dens;
        k2= GammaM * (prob_dens + 0.5*ts*k1 );
        k3= GammaM * (prob_dens + 0.5*ts*k2 );
        k4= GammaM * (prob_dens + ts*k3 );
        prob_dens = prob_dens +  ts/6.0 * ( k1 + 2.0*k2 + 2.0*k3 +k4 ) ;

        //Propagation using forward Euler
        //prob_dens += ts*GammaM*prob_dens;
        
        prob_dens.for_each( [](vec::elem_type& val) {
            if (val < 0){
                val=0;
            }
        });
        
        if (t % normalization_step==0) {
            //normalization = norm2D_spec(angle_grid, mom_grid, prob_dens);
            //Normalization using Fubini rule for double integration
            for (long long int i=0; i < angle_points; i++){
                integral_Y1 = prob_dens.subvec(i*mom_points, (i+1)*mom_points -1);
                integral_Y2 = prob_dens.subvec(i*mom_points + number_of_points, angle_points*mom_points + (i+1)*mom_points -1);
                mat int_trap = trapz( Y, integral_Y1);
                integral_X1(i)= int_trap(0,0);
                int_trap = trapz( Y, integral_Y2);
                integral_X2(i)= int_trap(0,0) ;
            }
            mat int_trap = trapz( X, integral_X1);
            mat int_trap2 = trapz( X, integral_X2);
            normalization = as_scalar( int_trap + int_trap2 ) ;
            prob_dens = prob_dens/ normalization;
        }
        
        if (t % restart_step==0){
            time = initial_time + t*ts;
            restart.str("");
            restart << time;
            print_restart();
            cout << "+ Time: " << time << " fs" << endl;
        }
        

    }
    auto end_for = chrono::high_resolution_clock::now();
    auto duration_for = chrono::duration_cast<chrono::milliseconds>(end_for - start_for);
    duration_evo = duration_for.count();

    cout << "*************************" << endl;
    cout << "* End of time evolution *" << endl;
    cout << "*************************" << endl << endl;

}






/************************/
/* Debugging & printing */
/************************/
void QSLE::print_m(){
    string filename_m01 = base_name + "_m01.dat";
    string filename_m10 = base_name + "_m10.dat";

    ofstream fout_m01(filename_m01);
    fout_m01 << setprecision(16);
    for (long long int i=0; i < number_of_points; i++){
        fout_m01 << angle_grid[i/ mom_points] << " " <<  mom_grid[i % mom_points] << " " << m01[i] << endl;
        if (i % mom_points == mom_points-1){
            fout_m01 << endl ;
        }
    }
    fout_m01.close();
    ofstream fout_m10(filename_m10);
    fout_m10 << setprecision(16);
    for (long long int i=0; i < number_of_points; i++){
        fout_m10 << angle_grid[i/ mom_points] << " " <<  mom_grid[i % mom_points] << " " << m10[i] << endl;
        if (i % mom_points == mom_points-1){
            fout_m10 << endl ;
        }
    }
    fout_m10.close();

    cout << "=======================================================================" << endl;
    cout << "Transition rates saved in:" << endl;
    cout << "- " << base_name << "_m01.dat" << endl;
    cout << "- " << base_name << "_m10.dat" << endl;
    cout << "=======================================================================" << endl << endl;;
}

void QSLE::print_restart(){
    prob_dens.save(base_name + "_time_" + restart.str() + ".dat", raw_ascii);
    string filename_R = base_name + "_restart.inp";
    ofstream fout_R(filename_R);
    fout_R << setprecision(16);
    fout_R << "Prefix               " << base_name << endl;
    fout_R << "Restart              " << restart.str() << endl;
    fout_R << "Pes                  " << PES_file << endl;
    fout_R << "d01                  " << d01_file << endl;
    fout_R << "pdb                  " << pdb_file << endl;
    fout_R << "refatoms             " << dihedral[0] << " "  << dihedral[1] << " "  << dihedral[2] << " "  << dihedral[3] << endl;
    fout_R << "Angle_points         " << angle_points << endl;
    fout_R << "Momentum_points      " << mom_points << endl;
    fout_R << "Derivation_method    " << method << endl;
    if (method == "LRBF"){
        fout_R << "Angle_domain         " << close_q << endl;
        fout_R << "Momentum_domain      " << close_p << endl;
    }
    fout_R << "Timestep             " << ts << endl;
    fout_R << "Total_steps          " << total_steps << endl;
    fout_R << "Restart_step         " << restart_step << endl;
    fout_R << "Normalization_step   " << normalization_step << endl;
    fout_R << "Temperature          " << temperature << endl;
    fout_R << "Mom_max              " << mom_max << endl;
    fout_R << "Ek_max               " << max_Ek << endl;
    fout_R << "Friction             " << csi * DALTON * SPACE_UNIT * SPACE_UNIT / TIME_UNIT << endl;
    fout_R << "Decoherence_time     " << tau_dec << endl;
    if (distribution == "boltzmann"){
        fout_R << "Initial_distribution Boltzmann " << endl;
    }
    else if (distribution == "dirac"){
        fout_R << "Initial_distribution Dirac " << start_angle << " " << start_mom << endl;
    }
    else{
        fout_R << "Initial_distribution " << distribution << endl;
    }
    fout_R << "Force_transition     " << transition << endl;
    fout_R << "Reff                 " << reff << endl;
    fout_R << "C                    " << C << endl;
    fout_R << "Viscosity            " << viscosity << endl;
    fout_R.close();
}

void QSLE::print_g0(){
    string filename_g0 = base_name + "_g0.dat";
    ofstream fout_g0(filename_g0);
    fout_g0 << setprecision(16);
    for (long long int i=0; i < angle_points; i++){
        fout_g0 << angle_grid[i] << " " << g0[i] << endl;
    }
    fout_g0.close();

    cout << "=======================================================================" << endl;
    cout << "Generalized mass vector saved in " << base_name << "_g0.dat" << endl;
    cout << "=======================================================================" << endl << endl;;

}

/*
void QSLE::print_dg0(){
    string filename_dg0 = base_name + "_dg0.dat";
    ofstream fout_dg0(filename_dg0);
    fout_dg0 << setprecision(16);
    for (long long int i=0; i < angle_points; i++){
        fout_dg0 << angle_grid[i] << " " << dg0[i] << endl;
    }
    fout_dg0.close();
}

void QSLE::print_testM(){
    string filename_testM = base_name + "_testM.dat";
    ofstream fout_testM(filename_testM);
    fout_testM << setprecision(16);
    for (int i=0; i < testM.size(); i++){
        fout_testM << tempo[i] << " " << testM[i] << endl;
    }
    fout_testM.close();
}

void QSLE::print_forces(){
    string filename_FG = base_name + "_FG.dat";
    string filename_FE = base_name + "_FE.dat";
    ofstream fout_FG(filename_FG);
    fout_FG << setprecision(16);
    for (long long int i=0; i < angle_points; i++){
        fout_FG << angle_grid[i] << " " << FG[i] << endl;
    }
    fout_FG.close();
    ofstream fout_FE(filename_FE);
    fout_FE << setprecision(16);
    for (long long int i=0; i < angle_points; i++){
        fout_FE << angle_grid[i] << " " << FE[i] << endl;
    }
    fout_FE.close();
}
*/
