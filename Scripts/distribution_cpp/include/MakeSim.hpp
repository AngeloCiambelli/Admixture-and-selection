#ifndef DEF_MAKESIM_HPP
#define DEF_MAKESIM_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>
#include <Python.h>
#include <random>
#include <omp.h>

#include "matplotlibcpp.h" // Display with matplotlib
#include "Vecteur.hpp"
#include "Admixture.hpp"
#include "Simulation.hpp"


namespace plt = matplotlibcpp;
using namespace std;

//===========================================================================
//                          Description
//===========================================================================
//
// Class MakeSim : Make N simulation of addmixture with random mating and selection 
// given a distribution of contribution vector (contribution) and a distribution of selection matrix (Phi) .
// to randomly draw parameters for the N simulations. 
//
// Save the parameters in a table under the simulation id. 
// Then save under the same simulation id in another table the result summary statistics.
//
//===========================================================================
//                          Define MakeSim class
//===========================================================================

class MakeSim
{
    public:
    string scenario;

    //MakeSim(string scen, int N);
    MakeSim(Vecteur<string> scen_list, int N, int nStep, int nThreads);

};

float f(float x, float x_min, float x_max, float y_min, float y_max, float a)
{
    float val = a * (y_max-y_min) * (1 - (x-x_min)/(x_max-x_min)) / (a + (x-x_min)/(x_max-x_min)) + y_min;
    return val;
}

float g(float x, float x_min, float x_max, float y_min, float y_max, float a)
{
    float val = a * (y_max-y_min) * ((x-x_min)/(x_max-x_min)) / (a + 1 - (x-x_min)/(x_max-x_min)) + y_min;
    return val;
}

MakeSim::MakeSim(Vecteur<string> scenario_list, int N, int nStep, int nThreads)
{
    string fileName = "./output/" + scenario +"/";
    int nGen = 20;
    int nPop = 2;

    #pragma omp parallel for num_threads(nThreads)
    for (const auto& scenario : scenario_list)
    {
        //#pragma omp parallel for num_threads(nThreads)
        for (int i = 1; i <= N; i++) 
        {
            //Contribution vector
            Vecteur<Vecteur<float>> contribution(1,{0.,0.,0.});                       
            for (int i=0; i<nGen; i++){contribution.push_back({0.,0.,0.});}

            //Get simulation coefficients in scenario/simu_i/simu_i.txt
            float ngen, Ne, s1_i, s2_i;
            fstream file;
            string param_file = "./output/parameters/"+scenario+"/simu_"+to_string(i)+"/simu_"+to_string(i)+".par";
            file.open(param_file, ios::in);
            string line;
            getline(file, line);
            while(file >> ngen >> Ne >> s1_i >> s2_i)
            {   
                if (int(ngen)<=nGen)
                {
                    contribution[int(ngen)][0] = s1_i;
                    contribution[int(ngen)][1] = s2_i;
                    contribution[int(ngen)][2] = 1 - s1_i - s2_i;
                }
            }
            file.close();

            //Get simualtion id
            string id = scenario+"_simu_"+to_string(i);

            //Create simulation object
            Simulation sim;

            //Draw selection parameters and run the corresponding simulation
            if (scenario == "AfrDE-EurDE-SelIN") 
            {
                //Set random device
                std::random_device rd{};
                std::mt19937 gen{rd()}; //Using generation 32-bit Mersenne Twister by Matsumoto and Nishimura, 1998 (one of the best)
                normal_distribution<float> dist_phi(0,1);
                float ln_Phi11_20 = dist_phi(gen); float ln_Phi22_20 = dist_phi(gen);
                float ln_Phi11_1 = dist_phi(gen); float ln_Phi22_1 = dist_phi(gen);

                //Redraw parameters if they are not compatible with scenario
                while ((abs(3.f*ln_Phi11_1) > abs(ln_Phi11_20)) || (ln_Phi11_20 * ln_Phi11_1 < 0))
                {
                    ln_Phi11_1 = dist_phi(gen);
                }
                while ((abs(3.f*ln_Phi22_1) > abs(ln_Phi22_20)) || (ln_Phi22_20 * ln_Phi22_1 < 0))
                {
                    ln_Phi22_1 = dist_phi(gen);
                }

                //Compute Phi parameters (coeff at gen 1 and 20)
                float Phi11_1 = exp(ln_Phi11_1); float Phi22_1 = exp(ln_Phi22_1);
                float Phi11_20 = exp(ln_Phi11_20); float Phi22_20 = exp(ln_Phi22_20);

                //Draw stepness of increase
                uniform_real_distribution<float> dist_u(0.f,0.5f);
                float u1 = dist_u(gen); float u2 = dist_u(gen);
                float a1 = pow(u1,2)/(1-2*u1); float a2 = pow(u2,2.f)/(1.f-2.f*u2);

                //Make the matrix Phi 1 and 20 and increasing parameter
                Vecteur<Vecteur<Vecteur<float>>> Phi(nGen, Vecteur<Vecteur<float>> ({Vecteur<float>({1.,1.}),Vecteur<float>({1.,1.})})); 
                for (int i=0; i<nGen; i++)
                {
                    float Phi_11_i = g(i+1, 1, 20, Phi11_1, Phi11_20, a1);
                    float Phi_22_i = g(i+1, 1, 20, Phi22_1, Phi22_20, a2);

                    Phi[i][0][0] = Phi_11_i;
                    Phi[i][1][1] = Phi_22_i;
                }

                //Save the phi parameters
                string param_phi_file = "./output/our_parameters/"+scenario+"/simu_"+to_string(i)+"/simu_"+to_string(i)+"_phi.par";
                ofstream phi_file;
                phi_file.open(param_phi_file);
                phi_file << "Phi11;" << "Phi12;" << "Phi21;" << "Phi22\n";
                for (int i=0; i<nGen; i++)
                {
                    phi_file << Phi[i][0][0] << ";" << Phi[i][0][1] << ";"  << Phi[i][1][0] << ";" << Phi[i][1][1] << "\n";
                }
                phi_file.close();

                //Run simulation
                sim = Simulation (contribution, Phi, nGen, fileName, nStep);
            }
            else if (scenario == "AfrDE-EurDE-SelDE") 
            {
                //Set random device
                std::random_device rd{};
                std::mt19937 gen{rd()}; //Using generation 32-bit Mersenne Twister by Matsumoto and Nishimura, 1998 (one of the best)
                normal_distribution<float> dist_phi(0,1);
                float ln_Phi11_20 = dist_phi(gen); float ln_Phi22_20 = dist_phi(gen);
                float ln_Phi11_1 = dist_phi(gen); float ln_Phi22_1 = dist_phi(gen);

                //Redraw parameters if they are not compatible with scenario
                while ((abs(3.f*ln_Phi11_20) > abs(ln_Phi11_1)) || (ln_Phi11_20 * ln_Phi11_1 < 0))
                {
                    ln_Phi11_20 = dist_phi(gen);
                }
                while ((abs(3.f*ln_Phi22_20) > abs(ln_Phi22_1)) || (ln_Phi22_20 * ln_Phi22_1 < 0) )
                {
                    ln_Phi22_20 = dist_phi(gen);
                }
                
                //Draw Phi parameters (coeff at gen 1 and 20)
                float Phi11_1 = exp(ln_Phi11_1); float Phi22_1 = exp(ln_Phi22_1);
                float Phi11_20 = exp(ln_Phi11_20); float Phi22_20 = exp(ln_Phi22_20);

                //Draw stepness of increasing parameter
                uniform_real_distribution<float> dist_u(0.f,0.5f);
                float u1 = dist_u(gen); float u2 = dist_u(gen);
                float a1 = pow(u1,2)/(1-2*u1); float a2 = pow(u2,2.f)/(1.f-2.f*u2);

                //Make the matrix Phi 1 and 20 and decreasing parameters
                Vecteur<Vecteur<Vecteur<float>>> Phi(nGen, Vecteur<Vecteur<float>> ({Vecteur<float>({1.,1.}),Vecteur<float>({1.,1.})})); 
                for (int i=0; i<nGen; i++)
                {
                    float Phi_11_i = f(i+1, 1, 20, Phi11_20, Phi11_1, a1);
                    float Phi_22_i = f(i+1, 1, 20, Phi22_20, Phi22_1, a2);

                    Phi[i][0][0] = Phi_11_i;
                    Phi[i][1][1] = Phi_22_i;
                }

                //Save the phi parameters
                string param_phi_file = "./output/our_parameters/"+scenario+"/simu_"+to_string(i)+"/simu_"+to_string(i)+"_phi.par";
                ofstream phi_file;
                phi_file.open(param_phi_file);
                phi_file << "Phi11;" << "Phi12;" << "Phi21;" << "Phi22\n";
                for (int i=0; i<nGen; i++)
                {
                    phi_file << Phi[i][0][0] << ";" << Phi[i][0][1] << ";"  << Phi[i][1][0] << ";" << Phi[i][1][1] << "\n";
                }
                phi_file.close();

                //Run simulation
                sim = Simulation (contribution, Phi, nGen, fileName, nStep);
            }
            else if (scenario == "AfrDE-EurDE-SelPC")
            {
                std::random_device rd{};
                std::mt19937 gen{rd()}; 

                //Draw selection change times
                uniform_int_distribution<int> dist_t(1,nGen);
                int t1 = dist_t(gen); int t2 = dist_t(gen);
                while (t1>=t2)
                {
                    t1 = dist_t(gen); t2 = dist_t(gen);
                }

                normal_distribution<float> dist_phi(0,1);
                float Phi11_0 = exp(dist_phi(gen)); float Phi22_0 = exp(dist_phi(gen));
                float Phi11_1 = exp(dist_phi(gen)); float Phi22_1 = exp(dist_phi(gen));
                float Phi11_2 = exp(dist_phi(gen)); float Phi22_2 = exp(dist_phi(gen));

                //Make Phi list
                Vecteur<Vecteur<Vecteur<float>>> Phi(nGen, Vecteur<Vecteur<float>> ({Vecteur<float>({1.,1.}),Vecteur<float>({1.,1.})})); 
                for (int i=0; i<nGen; i++)
                {
                    if (i+1 < t1)
                    {
                        Phi[i][0][0] = Phi11_0; Phi[i][1][1] = Phi22_0;
                    }
                    else if (i+1 >= t1 && i+1 < t2)
                    {
                        Phi[i][0][0] = Phi11_1; Phi[i][1][1] = Phi22_1;
                    }
                    else if (i+1 >= t2)
                    {
                        Phi[i][0][0] = Phi11_2; Phi[i][1][1] = Phi22_2;
                    }
                }

                //Save the phi parameters
                string param_phi_file = "./output/our_parameters/"+scenario+"/simu_"+to_string(i)+"/simu_"+to_string(i)+"_phi.par";
                ofstream phi_file;
                phi_file.open(param_phi_file);
                phi_file << "Phi11;" << "Phi12;" << "Phi21;" << "Phi22\n";
                for (int i=0; i<nGen; i++)
                {
                    phi_file << Phi[i][0][0] << ";" << Phi[i][0][1] << ";"  << Phi[i][1][0] << ";" << Phi[i][1][1] << "\n";
                }
                phi_file.close();

                //Run simulation
                sim = Simulation (contribution, Phi, nGen, fileName, nStep);
            }   
            else
            {
                //Run the simulations if no selection
                sim = Simulation (contribution, nGen, fileName, nStep);
            }

            //Check simulations run well
            cout << scenario << " : end sim" + to_string(i) << endl;

            //Create the directory if it don't exist
            string path = "./output/fortes_lima_nStep="+to_string(nStep)+"/"+scenario;
            //string path = "./output/Selection_nStep="+to_string(nStep)+"/"+scenario;
            filesystem::create_directories(path);

            //Convert to H1 to H to use other functions of Admixture
            Admixture converted_pool = sim.final_admixture_proportion.convert_fraction_H1_to_fraction();
            map<string, float> sumStat = converted_pool.summaryStat(10);
            ofstream newfile;
            newfile.open("./output/fortes_lima_nStep="+to_string(nStep)+"/"+scenario+"/"+scenario+"sumStat.txt", ios::app);
            //newfile.open("./output/Selection_nStep="+to_string(nStep)+"/"+scenario+"/"+scenario+"sumStat.txt", ios::app);
            newfile << id << ";" << sumStat["mean"] << ";" 
                    << sumStat["variance"] << ";" 
                    << sumStat["skewness"] << ";" 
                    << sumStat["kurtosis"] << ";" 
                    << sumStat["mode"] << ";" 
                    << sumStat["min"] << ";" 
                    << sumStat["max"] << ";" 
                    << sumStat["Q1"] << ";" 
                    << sumStat["Q2"] << ";" 
                    << sumStat["Q3"] << ";" 
                    << sumStat["Q4"] << ";" << sumStat["Q5"] 
                    << ";" << sumStat["Q6"] << ";" 
                    << sumStat["Q7"] << ";" 
                    << sumStat["Q8"] << ";" 
                    << sumStat["Q9"] << "\n";
            newfile.close();
        }

        //Check simulations run well
        cout << scenario << " : end sim" << endl;
    }
}

#endif
