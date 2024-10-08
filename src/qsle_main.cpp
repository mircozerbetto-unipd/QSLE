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

#include <iostream>
#include "QSLE_parallel.h"
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include "armadillo"

using namespace std;

int main(int argc, char* argv[]) {
    auto start_time = chrono::high_resolution_clock::now();
    srand(time(0));
          
cout << setprecision(16) << endl;

cout << "                           *****************" << endl;
cout << "                           * Q   S   L   E *" << endl;
cout << "                           *****************" << endl << endl;
/*
cout << "      _  _  _  _         _  _  _  _       _                 _  _  _  _  _ " << endl;                        
cout << "    _(_)(_)(_)(_)_     _(_)(_)(_)(_)_    (_)               (_)(_)(_)(_)(_)" << endl;                       
cout << "   (_)          (_)   (_)          (_)   (_)               (_)            " << endl;                    
cout << "   (_)          (_)   (_)_  _  _  _      (_)               (_) _  _       " << endl;                    
cout << "   (_)     _    (_)     (_)(_)(_)(_)_    (_)               (_)(_)(_)      " << endl;                    
cout << "   (_)    (_) _ (_)    _           (_)   (_)               (_)            " << endl;                    
cout << "   (_)_  _  _(_) _    (_)_  _  _  _(_)   (_) _  _  _  _    (_) _  _  _  _ " << endl;                    
cout << "     (_)(_)(_)  (_)     (_)(_)(_)(_)     (_)(_)(_)(_)(_)   (_)(_)(_)(_)(_)" << endl << endl;  


cout << "                                       _          _  _    " << endl;
cout << "                                    _ (_)       _(_)(_)_  " << endl;
cout << "                    _       _      (_)(_)      (_)    (_) " << endl;
cout << "                   (_)_   _(_)        (_)     (_)      (_)" << endl;
cout << "                     (_)_(_)    _     (_)  _   (_)_  _(_) " << endl;
cout << "                       (_)     (_)    (_) (_)    (_)(_)   " << endl << endl;
*/



cout << "    ______    ______   __        ________                            " << endl;
cout << "   /      \\  /      \\ |  \\      |        \\        __        ______   " << endl;
cout << "  |  $$$$$$\\|  $$$$$$\\| $$      | $$$$$$$$      _/  \\      /      \\  " << endl;
cout << "  | $$  | $$| $$___\\$$| $$      | $$__         |   $$     |  $$$$$$\\ " << endl;
cout << "  | $$  | $$ \\$$    \\ | $$      | $$  \\         \\$$$$     | $$  | $$ " << endl;
cout << "  | $$__| $$ _\\$$$$$$\\| $$      | $$$$$          | $$     | $$  | $$ " << endl;
cout << "  | $$  \\ $$|  \\__| $$| $$_____ | $$_____       _| $$_  __| $$__| $$ " << endl;
cout << "   \\$$ $$ $$ \\$$    $$| $$     \\| $$     \\     |   $$ \\|  \\\\$$    $$ " << endl;
cout << "    \\$$$$$$\\  \\$$$$$$  \\$$$$$$$$ \\$$$$$$$$      \\$$$$$$ \\$$ \\$$$$$$  " << endl;
cout << "        \\$$$                                                         " << endl << endl << endl;


cout << "                #######################################" << endl;
cout << "                #                                     #" << endl;
cout << "                #            < T | C | G >            #" << endl;
cout << "                #                                     #" << endl;
cout << "                #   Department of Chemical Sciences   #" << endl;
cout << "                #         University of Padua         #" << endl;
cout << "                #       Via Francesco Marzolo 1       #" << endl;
cout << "                #             35131 Padua             #" << endl;
cout << "                #                Italy                #" << endl;
cout << "                #                                     #" << endl;
cout << "                #######################################" << endl << endl << endl;

cout << "***********" << endl;
cout << "* Authors *" << endl;
cout << "***********" << endl;
cout << "- Riccardo Cortivo   : riccardo.cortivo@phd.unipd.it" << endl;
cout << "- Mirco Zerbetto     : mirco.zerbetto@unipd.it" << endl << endl << endl;


cout << "********************" << endl;
cout << "* Citation Details *" << endl;
cout << "********************" << endl;
cout << "- Jonathan Campeggio, Riccardo Cortivo and Mirco Zerbetto" << endl;
cout << "  A multiscale approach to coupled nuclear and electronic dynamics. I. Quantum-stochastic Liouville equation in natural internal coordinates" << endl;
cout << "  The Journal of Chemical Physics, Vol. 157, pp. 244104, 2023" << endl;
cout << "- Riccardo Cortivo, Jonathan Campeggio and Mirco Zerbetto" << endl;
cout << "  A multiscale approach to coupled nuclear and electronic dynamics. II. Exact and approximated evaluation of nonradiative transition rates" << endl;
cout << "  The Journal of Chemical Physics, Vol. 158, pp. 244105, 2023" << endl;
cout << "- Jonathan Campeggio, Antonino Polimeno and Mirco Zerbetto" << endl;
cout << "  DiTe2: Calculating the diffusion tensor for flexible molecules" << endl;
cout << "  The Journal of Computational Chemistry, Vol. 40, pp. 697, 2018" << endl;
cout << "- Conrad Sanderson and Ryan Curtin" << endl; 
cout << "  Armadillo: a template-based C++ library for linear algebra" << endl; 
cout << "  Journal of Open Source Software, Vol. 1, No. 2, pp. 26, 2016" << endl;  
cout << "- Conrad Sanderson and Ryan Curtin" << endl; 
cout << "  Practical Sparse Matrices in C++ with Hybrid Storage and Template-Based Expression Optimisation" << endl;
cout << "  Mathematical and Computational Applications, Vol. 24, No. 3, 2019" << endl << endl << endl;





    // Reading input file
    if (argc != 2){
        cout << endl << "ERROR >> Missing input filename" << endl;
        cout << endl <<  "Terminating QSLE because of errors" << endl << endl;
        exit(1);
    }

    // Initialization of the vectors
    QSLE QSLE_4BEAD(argv[1]);
    // First run section
    if (QSLE_4BEAD.restart_check() == false ) {
        // Setting the vectors using the input files
        QSLE_4BEAD.set_potential_vectors();
        /*MZ 05.09.2024*/QSLE_4BEAD.set_mu01();
        /*RC 01.10.2024*/QSLE_4BEAD.set_wrad();
        /*RC 01.10.2024*/QSLE_4BEAD.set_wnorad();
        QSLE_4BEAD.set_nonadiabatic_coeff();
        QSLE_4BEAD.set_friction();

        // Print files for debugging
        //QSLE_4BEAD.print_forces();
        //QSLE_4BEAD.print_dg0();
        QSLE_4BEAD.print_g0();

        // Compute and print transition coefficients
        QSLE_4BEAD.compute_transition_coefficients();
        cout << endl ;
        QSLE_4BEAD.print_m();

        // Compute and save Fokker-Planck matrix
        QSLE_4BEAD.compute_QSLE_matrix();
        QSLE_4BEAD.set_initial_condition();

    }
    // Restart section
    else{
        QSLE_4BEAD.recover_data();
    }
    // Time evolution
    QSLE_4BEAD.time_evolution();


    auto end_time = chrono::high_resolution_clock::now();
    auto duration_tot = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    long long int duration = duration_tot.count();

    long long int days, hours, minutes, seconds, msec;
    days = duration / 86400000;
    hours = (duration - days*86400000) / 3600000;
    minutes = (duration - days*86400000 - hours*3600000) / 60000;
    seconds = (duration - days*86400000 - hours*3600000 - minutes * 60000) / 1000;
    msec = duration - days*86400000 - hours*3600000 - minutes * 60000 - seconds *1000;

    cout << "***********************************" << endl;
    cout << "* Timings for individual modules: *" << endl;
    cout << "***********************************" << endl;
    cout << "- Total run time         ...   " << duration/1000.0 << " sec " << endl;
    cout << "- DiTe2                  ...   " << QSLE_4BEAD.duration_dite2/1000.0 << " sec --> " << round((QSLE_4BEAD.duration_dite2/1000.0)/(duration/1000.0)*10000.0)/100.0 << "%" << endl;
    cout << "- Transition rates       ...   " << QSLE_4BEAD.duration_m/1000.0 <<" sec --> " << round((QSLE_4BEAD.duration_m/1000.0)/(duration/1000.0)*10000.0)/100.0 << "%" << endl;
    cout << "- Fokker-Planck matrix   ...   " << QSLE_4BEAD.duration_FP/1000.0 << " sec --> " << round((QSLE_4BEAD.duration_FP/1000.0)/(duration/1000.0)*10000.0)/100.0 << "%" << endl;
    cout << "- Time evolution         ...   " << QSLE_4BEAD.duration_evo/1000.0 << " sec --> " << round((QSLE_4BEAD.duration_evo/1000.0)/(duration/1000.0)*10000.0)/100.0 << "%" << endl<< endl;


    cout << "=======================================================================" << endl << endl;
    cout << "TOTAL RUN TIME: " << days << " days " <<  hours << " hours " <<  minutes << " minutes " << seconds << " seconds "  << msec <<  " msec " << endl << endl;
    cout << "QSLE TERMINATED NORMALLY" << endl << endl;

    int random_value = rand() % 7;
    if (random_value == 0) {
        cout << "Yesterday is history, tomorrow is a mystery, but today is a gift. That\nis why it is called the present. –Master Oogway" << endl;
    }
    else if (random_value == 1){
        cout << "You have to believe in yourself. That’s the secret. – Po " << endl;
    }
    else if (random_value == 2){
        cout << "The strongest of us sometimes have the hardest time fighting what’s on\nthe inside. – Po " << endl;
    }
    else if (random_value == 3){
        cout << "The only true limit is the one you set for yourself. – Shifu" << endl;
    }
    else if (random_value == 4){
        cout << "If you only do what you can do, you’ll never be better than what you\nare. – Shifu" << endl;
    }
    else if (random_value == 5){
        cout << "The measure of a real champion is not whether they can triumph, but\nwhether they can overcome defeat. – Shifu" << endl;
    }
    else if (random_value == 6){
        cout << "Your mind is like water, my friend. When it is agitated, it becomes"<< endl;
        cout << "difficult to see. But if you allow it to settle, then the answer\nbecomes clear. – Master Oogway" << endl;
    }
    cout << endl;


    return 0;
}
